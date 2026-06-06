#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from typing import Dict, Tuple

import joblib
import pandas as pd
from scipy.stats import loguniform, randint, uniform
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, f1_score, matthews_corrcoef
from sklearn.model_selection import RandomizedSearchCV, train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

try:
    from xgboost import XGBClassifier
    HAS_XGBOOST = True
except Exception:
    HAS_XGBOOST = False


@dataclass
class RunResult:
    model_name: str
    best_params: Dict
    cv_mcc: float
    test_mcc: float
    test_accuracy: float
    test_f1: float


def load_dataset(path: str) -> Tuple[pd.DataFrame, pd.Series]:
    if path.endswith(".csv"):
        df = pd.read_csv(path)
    else:
        df = pd.read_csv(path, sep="\t")

    if "label" not in df.columns:
        raise ValueError("Input data must contain a 'label' column.")

    feature_cols = [c for c in df.columns if c not in {"label", "sequence"}]
    x = df[feature_cols].copy()
    y = df["label"].astype(int)
    return x, y


def split_data(x, y, random_state):
    x_trainval, x_test, y_trainval, y_test = train_test_split(
        x, y, test_size=0.20, stratify=y, random_state=random_state
    )

    x_train, x_val, y_train, y_val = train_test_split(
        x_trainval,
        y_trainval,
        test_size=0.25,
        stratify=y_trainval,
        random_state=random_state,
    )

    return {
        "train": (x_train, y_train),
        "val": (x_val, y_val),
        "test": (x_test, y_test),
        "trainval": (x_trainval, y_trainval),
    }


def build_search_spaces(scale_pos_weight: float, random_state: int):
    searches = {}

    searches["svm_linear"] = {
        "estimator": Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            ("clf", SVC(kernel="linear", class_weight="balanced", probability=True)),
        ]),
        "param_distributions": {
            "clf__C": loguniform(1e-3, 1e3),
        },
    }

    searches["svm_rbf"] = {
        "estimator": Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            ("clf", SVC(kernel="rbf", class_weight="balanced", probability=True)),
        ]),
        "param_distributions": {
            "clf__C": loguniform(1e-3, 1e3),
            "clf__gamma": loguniform(1e-4, 1),
        },
    }

    searches["logistic_regression"] = {
        "estimator": Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            ("clf", LogisticRegression(
                class_weight="balanced",
                max_iter=4000,
                solver="liblinear",
                random_state=random_state,
            )),
        ]),
        "param_distributions": {
            "clf__C": loguniform(1e-3, 1e3),
        },
    }

    searches["random_forest"] = {
        "estimator": Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("clf", RandomForestClassifier(
                class_weight="balanced",
                random_state=random_state,
                n_jobs=-1,
            )),
        ]),
        "param_distributions": {
            "clf__n_estimators": randint(200, 1201),
            "clf__max_depth": [None] + list(range(3, 31)),
            "clf__min_samples_leaf": randint(1, 11),
            "clf__min_samples_split": randint(2, 21),
        },
    }

    searches["knn"] = {
        "estimator": Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            ("clf", KNeighborsClassifier()),
        ]),
        "param_distributions": {
            "clf__n_neighbors": randint(3, 51),
            "clf__weights": ["uniform", "distance"],
            "clf__p": [1, 2],
        },
    }

    if HAS_XGBOOST:
        searches["xgboost"] = {
            "estimator": Pipeline([
                ("imputer", SimpleImputer(strategy="median")),
                ("clf", XGBClassifier(
                    random_state=random_state,
                    objective="binary:logistic",
                    eval_metric="logloss",
                    n_jobs=-1,
                    scale_pos_weight=scale_pos_weight,
                )),
            ]),
            "param_distributions": {
                "clf__n_estimators": randint(200, 1201),
                "clf__max_depth": randint(3, 13),
                "clf__min_child_weight": randint(1, 11),
                "clf__learning_rate": uniform(0.01, 0.29),
                "clf__subsample": uniform(0.6, 0.4),
                "clf__colsample_bytree": uniform(0.6, 0.4),
            },
        }

    return searches


def tune_and_evaluate(
    model_name,
    estimator,
    param_distributions,
    x_trainval,
    y_trainval,
    x_test,
    y_test,
    n_iter,
    random_state,
):
    search = RandomizedSearchCV(
        estimator=estimator,
        param_distributions=param_distributions,
        n_iter=n_iter,
        scoring="matthews_corrcoef",
        cv=5,
        n_jobs=-1,
        random_state=random_state,
        refit=True,
        verbose=0,
    )
    search.fit(x_trainval, y_trainval)

    y_pred = search.best_estimator_.predict(x_test)

    result = RunResult(
        model_name=model_name,
        best_params=search.best_params_,
        cv_mcc=float(search.best_score_),
        test_mcc=float(matthews_corrcoef(y_test, y_pred)),
        test_accuracy=float(accuracy_score(y_test, y_pred)),
        test_f1=float(f1_score(y_test, y_pred, zero_division=0)),
    )

    return result, search.best_estimator_


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to feature table (.tsv/.csv)")
    parser.add_argument("--n-iter", type=int, default=30)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--output-json", default="model_comparison_results.json")
    parser.add_argument("--output-model", default="best_model.pkl")
    parser.add_argument(
        "--save-model-name",
        default=None,
        help="Optional: save this specific model instead of the top-ranked model (e.g. xgboost)",
    )
    args = parser.parse_args()

    x, y = load_dataset(args.input)
    splits = split_data(x, y, random_state=args.random_state)

    x_train, y_train = splits["train"]
    x_val, y_val = splits["val"]
    x_test, y_test = splits["test"]
    x_trainval, y_trainval = splits["trainval"]

    print("Split sizes:")
    print(f"  Train:      {len(x_train)}")
    print(f"  Validation: {len(x_val)}")
    print(f"  Test:       {len(x_test)}")

    positives = int(y_trainval.sum())
    negatives = int((1 - y_trainval).sum())
    scale_pos_weight = negatives / max(positives, 1)

    searches = build_search_spaces(
        scale_pos_weight=scale_pos_weight,
        random_state=args.random_state,
    )

    results = []
    estimators = {}

    for name, spec in searches.items():
        print(f"Tuning {name}...")
        run_result, best_estimator = tune_and_evaluate(
            model_name=name,
            estimator=spec["estimator"],
            param_distributions=spec["param_distributions"],
            x_trainval=x_trainval,
            y_trainval=y_trainval,
            x_test=x_test,
            y_test=y_test,
            n_iter=args.n_iter,
            random_state=args.random_state,
        )
        results.append(run_result)
        estimators[name] = best_estimator

    if not HAS_XGBOOST:
        print("xgboost is not installed: XGBoost model was skipped.")

    results_sorted = sorted(
        results,
        key=lambda r: (r.test_mcc, -abs(r.cv_mcc - r.test_mcc)),
        reverse=True,
    )

    print("\nModel comparison (sorted by test MCC):")
    for r in results_sorted:
        gap = abs(r.cv_mcc - r.test_mcc)
        print(
            f"- {r.model_name:20s} | CV MCC: {r.cv_mcc:.4f} | "
            f"Test MCC: {r.test_mcc:.4f} | Test Acc: {r.test_accuracy:.4f} | "
            f"Test F1: {r.test_f1:.4f} | |CV-Test| gap: {gap:.4f}"
        )

    best = results_sorted[0]
    print(f"\nBest model: {best.model_name}")
    print(f"Best params: {best.best_params}")

    if args.save_model_name is not None:
        if args.save_model_name not in estimators:
            raise ValueError(
                f"Requested model '{args.save_model_name}' not found. "
                f"Available models: {list(estimators.keys())}"
            )
        model_to_save_name = args.save_model_name
    else:
        model_to_save_name = best.model_name

    model_to_save = estimators[model_to_save_name]
    joblib.dump(model_to_save, args.output_model)
    print(f"Saved fitted model '{model_to_save_name}' to: {args.output_model}")

    output_payload = {
        "input": args.input,
        "n_samples": len(x),
        "split_sizes": {
            "train": len(x_train),
            "validation": len(x_val),
            "test": len(x_test),
        },
        "results": [
            {
                "model": r.model_name,
                "best_params": r.best_params,
                "cv_mcc": r.cv_mcc,
                "test_mcc": r.test_mcc,
                "test_accuracy": r.test_accuracy,
                "test_f1": r.test_f1,
                "generalization_gap_abs": abs(r.cv_mcc - r.test_mcc),
            }
            for r in results_sorted
        ],
        "best_model": {
            "model": best.model_name,
            "best_params": best.best_params,
        },
    }

    with open(args.output_json, "w", encoding="utf-8") as f:
        json.dump(output_payload, f, indent=2)

    print(f"Saved results to: {args.output_json}")


if __name__ == "__main__":
    main()