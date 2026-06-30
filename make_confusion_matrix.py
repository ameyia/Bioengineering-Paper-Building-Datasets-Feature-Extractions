#!/usr/bin/env python3

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import ConfusionMatrixDisplay, confusion_matrix, matthews_corrcoef
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

try:
    from xgboost import XGBClassifier
    HAS_XGB = True
except ImportError:
    HAS_XGB = False


def load_table(path):
    sep = "," if str(path).endswith(".csv") else "\t"
    df = pd.read_csv(path, sep=sep)
    y = df["label"].astype(int)
    drop_cols = ["label"]
    if "sequence" in df.columns:
        drop_cols.append("sequence")
    X = df.drop(columns=drop_cols)
    X = X.apply(pd.to_numeric, errors="coerce")
    X = X.fillna(X.median(numeric_only=True))
    return X, y


def get_models(random_state):
    models = {
        "random_forest": RandomForestClassifier(
            n_estimators=500,
            min_samples_leaf=2,
            class_weight="balanced",
            random_state=random_state,
            n_jobs=-1,
        ),
        "logistic_regression": Pipeline(
            [
                ("scaler", StandardScaler()),
                (
                    "model",
                    LogisticRegression(
                        max_iter=5000,
                        class_weight="balanced",
                        random_state=random_state,
                    ),
                ),
            ]
        ),
        "svm_rbf": Pipeline(
            [
                ("scaler", StandardScaler()),
                (
                    "model",
                    SVC(
                        kernel="rbf",
                        C=10,
                        gamma="scale",
                        probability=True,
                        class_weight="balanced",
                        random_state=random_state,
                    ),
                ),
            ]
        ),
        "knn": Pipeline(
            [
                ("scaler", StandardScaler()),
                ("model", KNeighborsClassifier(n_neighbors=5)),
            ]
        ),
    }

    if HAS_XGB:
        models["xgboost"] = XGBClassifier(
            n_estimators=300,
            max_depth=3,
            learning_rate=0.03,
            subsample=0.8,
            colsample_bytree=0.8,
            eval_metric="logloss",
            random_state=random_state,
            n_jobs=-1,
        )

    return models


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--table", default="pcp_aac_ctd_dpc.tsv")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--model",
        default=None,
        help="Train one model by name, e.g. svm_rbf. If omitted, choose the best model by validation MCC.",
    )
    parser.add_argument("--out", default="confusion_matrix_pcp_aac_ctd_dpc_seed0.png")
    args = parser.parse_args()

    X, y = load_table(args.table)

    X_train, X_temp, y_train, y_temp = train_test_split(
        X,
        y,
        test_size=0.40,
        stratify=y,
        random_state=args.seed,
    )

    X_val, X_test, y_val, y_test = train_test_split(
        X_temp,
        y_temp,
        test_size=0.50,
        stratify=y_temp,
        random_state=args.seed,
    )

    all_models = get_models(args.seed)
    if args.model:
        if args.model not in all_models:
            raise ValueError(f"Unknown model '{args.model}'. Choose from: {sorted(all_models)}")
        models = {args.model: all_models[args.model]}
    else:
        models = all_models

    val_scores = {}

    for model_name, model in models.items():
        model.fit(X_train, y_train)
        val_pred = model.predict(X_val)
        val_scores[model_name] = matthews_corrcoef(y_val, val_pred)

    best_model_name = max(val_scores, key=val_scores.get)
    best_model = models[best_model_name]

    y_pred = best_model.predict(X_test)
    cm = confusion_matrix(y_test, y_pred, labels=[0, 1])

    disp = ConfusionMatrixDisplay(
        confusion_matrix=cm,
        display_labels=["Negative", "Antithrombotic"],
    )
    disp.plot(values_format="d", cmap="Blues")
    plt.title(f"Confusion Matrix: {Path(args.table).stem}, seed {args.seed}, {best_model_name}")
    plt.tight_layout()
    plt.savefig(args.out, dpi=300, bbox_inches="tight")

    print(f"Best model: {best_model_name}")
    print(f"Validation MCC: {val_scores[best_model_name]:.3f}")
    print(f"Confusion matrix [[TN, FP], [FN, TP]]:")
    print(cm)
    print(f"Saved: {args.out}")


if __name__ == "__main__":
    main()
