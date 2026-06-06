#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import os
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, Tuple, List

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

from Bio.SeqUtils.ProtParam import ProteinAnalysis
from propy.CTD import CalculateCTD

try:
    from xgboost import XGBClassifier
    HAS_XGBOOST = True
except Exception:
    HAS_XGBOOST = False


# =============================
# EXTERNAL DATASETS
# =============================

external_sequences_s3 = [
    ("SGGHQTAVPKISKQGLGGDFEEIPSDEIIE", 1),   # avathrin
    ("SDEAVRAIPKMYSTAPPGDFETIPDDAIEEREMKAR", 1),
    ("SDEAVRAIPKMYSTAPPGDFETIPDDAIEER", 1),
    ("SDEAVRAIPKMYSTAPPGDFETIPDDAIEE", 1),
    ("SDEAVRAIPKMYSTAPPGDFEEIPDDAIEE", 1),        # ultravariegin
    ("SDQGDVAIPKMYSTAPPGDFEEIPDDAIEE", 1),        # UV003
    ("SDEAVRAEPKMHKTAPPGDFEEIPDDAIEE", 1),        # UV004
    ("SDEAVRAIPKMYSTAPPGDFEEIPEEYLDDES", 1),      # UV005
    ("SDEAVRAIPKMYSTAPPGDFEEIPDDEIEE", 1),        # UV012
    ("SDEAVRAIPKMYSQAPPGDFEEIPDDAIEE", 1),        # UV013
    ("MYSTAPPGDFEEIPDDAIEE", 1),                  # UV011
]

external_sequences_s5 = [
    ("LVYTDCTESGQNLCLCEGSNVCGQGNKCILGSDGEKNQCVTGEGTPKPQSHNDGDFEEIPEEYLQ", 1),  # lepirudin
    ("VVYTDCTESGQNLCLCEGSNVCGQGNKCILGSDGEKNQCVTGEGTPKPQSHNDGDFEEIPEEYLQ", 1),  # desirudin
    ("FPRPGGGGNGDFEEIPEEYL", 1),  # bivalirudin
]

amino_acids = list("ACDEFGHIKLMNPQRSTVWY")
dipeptides = [a + b for a in amino_acids for b in amino_acids]


@dataclass
class RunResult:
    model_name: str
    best_params: Dict
    cv_mcc: float
    test_mcc: float
    test_accuracy: float
    test_f1: float
    generalization_gap_abs: float


def dataset_tag_from_path(path: str) -> str:
    return Path(path).stem


def load_dataset(path: str) -> Tuple[pd.DataFrame, pd.Series, List[str]]:
    if path.endswith(".csv"):
        df = pd.read_csv(path)
    else:
        df = pd.read_csv(path, sep="\t")

    if "label" not in df.columns:
        raise ValueError("Input data must contain a 'label' column.")

    feature_cols = [c for c in df.columns if c not in {"label", "sequence"}]
    x = df[feature_cols].copy()
    y = df["label"].astype(int)
    return x, y, feature_cols


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
        "param_distributions": {"clf__C": loguniform(1e-3, 1e3)},
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
        "param_distributions": {"clf__C": loguniform(1e-3, 1e3)},
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

    cv_mcc = float(search.best_score_)
    test_mcc = float(matthews_corrcoef(y_test, y_pred))
    test_acc = float(accuracy_score(y_test, y_pred))
    test_f1 = float(f1_score(y_test, y_pred, zero_division=0))

    result = RunResult(
        model_name=model_name,
        best_params=search.best_params_,
        cv_mcc=cv_mcc,
        test_mcc=test_mcc,
        test_accuracy=test_acc,
        test_f1=test_f1,
        generalization_gap_abs=abs(cv_mcc - test_mcc),
    )

    return result, search.best_estimator_


def featurize_external(sequence_label_pairs, feature_cols: List[str]):
    rows = []

    use_ctd = any("_" in c for c in feature_cols)
    use_dpc = any(len(c) == 2 and c in dipeptides for c in feature_cols)

    for seq, label in sequence_label_pairs:
        seq = seq.strip().upper()
        if len(seq) == 0 or "X" in seq:
            continue

        try:
            analysis = ProteinAnalysis(seq)
            length = len(seq)

            feature_dict = {
                "sequence": seq,
                "length": length,
                "mw": analysis.molecular_weight(),
                "aromaticity": analysis.aromaticity(),
                "pI": analysis.isoelectric_point(),
                "instability": analysis.instability_index(),
                "label": label,
            }

            for aa in amino_acids:
                feature_dict[aa] = (seq.count(aa) / length) * 100

            if use_ctd:
                try:
                    ctd = CalculateCTD(seq)
                    feature_dict.update(ctd)
                except Exception:
                    print(f"CTD failed for external sequence: {seq}")

            if use_dpc:
                dpc = dict.fromkeys(dipeptides, 0)
                total_dipeptides = len(seq) - 1
                if total_dipeptides > 0:
                    for i in range(total_dipeptides):
                        dp = seq[i:i+2]
                        if dp in dpc:
                            dpc[dp] += 1
                    for dp in dpc:
                        dpc[dp] = (dpc[dp] / total_dipeptides) * 100
                feature_dict.update(dpc)

            rows.append(feature_dict)

        except Exception as e:
            print(f"Error with external sequence: {seq}")
            print(e)

    return pd.DataFrame(rows)


def save_external_predictions(dataset_name, sequence_label_pairs, feature_cols, estimators, output_dir):
    df_external = featurize_external(sequence_label_pairs, feature_cols)
    X_external = df_external[feature_cols].copy()

    for model_name, model in estimators.items():
        probs = model.predict_proba(X_external)[:, 1]
        preds = model.predict(X_external)

        df_out = df_external.copy()
        df_out["predicted_label"] = preds
        df_out["prediction_probability"] = probs
        df_out = df_out.sort_values("prediction_probability", ascending=False)

        out_path = os.path.join(output_dir, f"{model_name}_{dataset_name}_predictions.tsv")
        df_out.to_csv(out_path, sep="\t", index=False)
        print(f"Saved {dataset_name} predictions: {out_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="One training feature table")
    parser.add_argument("--n-iter", type=int, default=30)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--models-root", default="saved_models")
    parser.add_argument("--predictions-root", default="external_predictions")
    parser.add_argument("--results-root", default="results_by_dataset")
    args = parser.parse_args()

    dataset_tag = dataset_tag_from_path(args.input)

    models_dir = os.path.join(args.models_root, dataset_tag)
    preds_dir = os.path.join(args.predictions_root, dataset_tag)
    results_dir = os.path.join(args.results_root, dataset_tag)

    os.makedirs(models_dir, exist_ok=True)
    os.makedirs(preds_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)

    x, y, feature_cols = load_dataset(args.input)
    splits = split_data(x, y, random_state=args.random_state)

    x_train, y_train = splits["train"]
    x_val, y_val = splits["val"]
    x_test, y_test = splits["test"]
    x_trainval, y_trainval = splits["trainval"]

    print(f"\n=== DATASET: {dataset_tag} ===")
    print("Split sizes:")
    print(f"  Train:      {len(x_train)}")
    print(f"  Validation: {len(x_val)}")
    print(f"  Test:       {len(x_test)}")

    positives = int(y_trainval.sum())
    negatives = int((1 - y_trainval).sum())
    scale_pos_weight = negatives / max(positives, 1)

    searches = build_search_spaces(scale_pos_weight=scale_pos_weight, random_state=args.random_state)

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

        model_path = os.path.join(models_dir, f"{name}.pkl")
        joblib.dump(best_estimator, model_path)
        print(f"Saved model: {model_path}")

    results_sorted = sorted(
        results,
        key=lambda r: (r.test_mcc, -r.generalization_gap_abs),
        reverse=True,
    )

    print("\nInternal S1 model comparison:")
    for r in results_sorted:
        print(
            f"- {r.model_name:20s} | CV MCC: {r.cv_mcc:.4f} | "
            f"Test MCC: {r.test_mcc:.4f} | Acc: {r.test_accuracy:.4f} | "
            f"F1: {r.test_f1:.4f} | Gap: {r.generalization_gap_abs:.4f}"
        )

    results_json = os.path.join(results_dir, "internal_results.json")
    with open(results_json, "w", encoding="utf-8") as f:
        json.dump(
            {
                "dataset": dataset_tag,
                "input": args.input,
                "n_samples": len(x),
                "split_sizes": {
                    "train": len(x_train),
                    "validation": len(x_val),
                    "test": len(x_test),
                },
                "results": [asdict(r) for r in results_sorted],
            },
            f,
            indent=2,
        )
    print(f"Saved internal results: {results_json}")

    save_external_predictions("s3", external_sequences_s3, feature_cols, estimators, preds_dir)
    save_external_predictions("s5", external_sequences_s5, feature_cols, estimators, preds_dir)


if __name__ == "__main__":
    main()
    