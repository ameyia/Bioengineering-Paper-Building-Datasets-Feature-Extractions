#!/usr/bin/env python3

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    ConfusionMatrixDisplay,
    accuracy_score,
    confusion_matrix,
    f1_score,
    matthews_corrcoef,
    precision_score,
    recall_score,
)
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

    if "label" not in df.columns:
        raise ValueError(f"{path} does not contain a label column")

    y = df["label"].astype(int)
    drop_cols = ["label"]
    if "sequence" in df.columns:
        drop_cols.append("sequence")

    X = df.drop(columns=drop_cols)
    X = X.apply(pd.to_numeric, errors="coerce")
    X = X.replace([np.inf, -np.inf], np.nan)
    X = X.fillna(X.median(numeric_only=True))
    return X, y


def get_models(random_state):
    models = {
        "svm_linear": Pipeline(
            [
                ("scaler", StandardScaler()),
                (
                    "model",
                    SVC(
                        kernel="linear",
                        probability=True,
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
        "random_forest": RandomForestClassifier(
            n_estimators=500,
            min_samples_leaf=2,
            class_weight="balanced",
            random_state=random_state,
            n_jobs=-1,
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


def split_data(X, y, seed):
    X_train, X_temp, y_train, y_temp = train_test_split(
        X,
        y,
        test_size=0.40,
        stratify=y,
        random_state=seed,
    )

    X_val, X_test, y_val, y_test = train_test_split(
        X_temp,
        y_temp,
        test_size=0.50,
        stratify=y_temp,
        random_state=seed,
    )

    return X_train, X_val, X_test, y_train, y_val, y_test


def metric_row(dataset, seed, model_name, y_val, val_pred, y_test, test_pred):
    tn, fp, fn, tp = confusion_matrix(y_test, test_pred, labels=[0, 1]).ravel()

    return {
        "dataset": dataset,
        "seed": seed,
        "model": model_name,
        "val_mcc": matthews_corrcoef(y_val, val_pred),
        "test_mcc": matthews_corrcoef(y_test, test_pred),
        "test_accuracy": accuracy_score(y_test, test_pred),
        "test_f1": f1_score(y_test, test_pred, zero_division=0),
        "test_precision": precision_score(y_test, test_pred, zero_division=0),
        "test_recall": recall_score(y_test, test_pred, zero_division=0),
        "tn": int(tn),
        "fp": int(fp),
        "fn": int(fn),
        "tp": int(tp),
    }


def save_confusion_matrix(cm, title, out_path):
    disp = ConfusionMatrixDisplay(
        confusion_matrix=cm,
        display_labels=["Negative", "Antithrombotic"],
    )
    disp.plot(values_format="d", cmap="Blues")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--tables",
        nargs="+",
        default=["data/processed/baseline/pcp_aac.tsv", "data/processed/baseline/pcp_aac_ctd.tsv", "data/processed/baseline/pcp_aac_ctd_dpc.tsv"],
        help="Feature tables to evaluate",
    )
    parser.add_argument("--seeds", type=int, default=10, help="Number of seeds, starting at 0")
    parser.add_argument("--outdir", default="results/figures/confusion_matrices")
    parser.add_argument(
        "--plot-mode",
        choices=["none", "best", "all"],
        default="best",
        help="none: only CSVs; best: plot validation-selected model per dataset/seed; all: plot every model",
    )
    args = parser.parse_args()

    outdir = Path(args.outdir)
    plot_dir = outdir / "plots"
    outdir.mkdir(parents=True, exist_ok=True)
    if args.plot_mode != "none":
        plot_dir.mkdir(parents=True, exist_ok=True)

    rows = []

    for table in args.tables:
        dataset = Path(table).stem
        print(f"\n=== {dataset} ===")
        X, y = load_table(table)

        for seed in range(args.seeds):
            print(f"  seed {seed}")
            X_train, X_val, X_test, y_train, y_val, y_test = split_data(X, y, seed)

            seed_results = []
            seed_predictions = {}
            models = get_models(seed)

            for model_name, model in models.items():
                print(f"    {model_name}")
                model.fit(X_train, y_train)
                val_pred = model.predict(X_val)
                test_pred = model.predict(X_test)

                row = metric_row(dataset, seed, model_name, y_val, val_pred, y_test, test_pred)
                seed_results.append(row)
                seed_predictions[model_name] = test_pred

            best_model = max(seed_results, key=lambda row: row["val_mcc"])["model"]

            for row in seed_results:
                row["selected_by_val_mcc"] = row["model"] == best_model
                rows.append(row)

                should_plot = args.plot_mode == "all" or (
                    args.plot_mode == "best" and row["model"] == best_model
                )
                if should_plot:
                    cm = np.array([[row["tn"], row["fp"]], [row["fn"], row["tp"]]])
                    out_path = plot_dir / f"{dataset}_seed{seed}_{row['model']}.png"
                    title = f"{dataset}, seed {seed}, {row['model']}"
                    save_confusion_matrix(cm, title, out_path)

    results_df = pd.DataFrame(rows)
    results_path = outdir / "confusion_matrix_results_by_model_seed.csv"
    results_df.to_csv(results_path, index=False)

    summary_df = (
        results_df.groupby(["dataset", "model"], as_index=False)
        .agg(
            mean_test_mcc=("test_mcc", "mean"),
            std_test_mcc=("test_mcc", "std"),
            mean_test_accuracy=("test_accuracy", "mean"),
            mean_test_f1=("test_f1", "mean"),
            mean_test_precision=("test_precision", "mean"),
            mean_test_recall=("test_recall", "mean"),
            mean_tn=("tn", "mean"),
            mean_fp=("fp", "mean"),
            mean_fn=("fn", "mean"),
            mean_tp=("tp", "mean"),
            times_selected_by_val_mcc=("selected_by_val_mcc", "sum"),
        )
        .sort_values(["dataset", "mean_test_mcc"], ascending=[True, False])
    )

    summary_path = outdir / "confusion_matrix_summary_by_dataset_model.csv"
    summary_df.to_csv(summary_path, index=False)

    selected_df = results_df[results_df["selected_by_val_mcc"]].copy()
    selected_summary_df = (
        selected_df.groupby("dataset", as_index=False)
        .agg(
            mean_test_mcc=("test_mcc", "mean"),
            std_test_mcc=("test_mcc", "std"),
            mean_test_accuracy=("test_accuracy", "mean"),
            mean_test_f1=("test_f1", "mean"),
            mean_tn=("tn", "mean"),
            mean_fp=("fp", "mean"),
            mean_fn=("fn", "mean"),
            mean_tp=("tp", "mean"),
        )
        .sort_values("mean_test_mcc", ascending=False)
    )

    selected_summary_path = outdir / "confusion_matrix_summary_validation_selected.csv"
    selected_summary_df.to_csv(selected_summary_path, index=False)

    print("\nSaved:")
    print(f"  {results_path}")
    print(f"  {summary_path}")
    print(f"  {selected_summary_path}")
    if args.plot_mode != "none":
        print(f"  plots in {plot_dir}")


if __name__ == "__main__":
    main()
