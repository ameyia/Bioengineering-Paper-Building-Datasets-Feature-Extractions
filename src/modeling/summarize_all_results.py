#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def safe_read_json(path: Path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def summarize_internal(results_root: Path) -> pd.DataFrame:
    rows = []

    if not results_root.exists():
        return pd.DataFrame()

    for dataset_dir in sorted(results_root.iterdir()):
        if not dataset_dir.is_dir():
            continue

        json_path = dataset_dir / "internal_results.json"
        if not json_path.exists():
            continue

        payload = safe_read_json(json_path)
        dataset = payload.get("dataset", dataset_dir.name)

        for r in payload.get("results", []):
            model_name = r.get("model", r.get("model_name"))
            if model_name is None:
                continue

            row = {
                "dataset": dataset,
                "model": model_name,
                "cv_mcc": r.get("cv_mcc"),
                "test_mcc": r.get("test_mcc"),
                "test_accuracy": r.get("test_accuracy"),
                "test_f1": r.get("test_f1"),
                "generalization_gap_abs": r.get("generalization_gap_abs"),
            }
            rows.append(row)

    return pd.DataFrame(rows)


def summarize_external(predictions_root: Path) -> pd.DataFrame:
    rows = []

    if not predictions_root.exists():
        return pd.DataFrame()

    for dataset_dir in sorted(predictions_root.iterdir()):
        if not dataset_dir.is_dir():
            continue

        dataset = dataset_dir.name

        for tsv_path in sorted(dataset_dir.glob("*_s3_predictions.tsv")):
            model = tsv_path.name.replace("_s3_predictions.tsv", "")
            df = pd.read_csv(tsv_path, sep="\t")

            if "predicted_label" not in df.columns:
                continue

            n_total = len(df)
            n_positive = int(df["predicted_label"].sum())
            n_negative = int((df["predicted_label"] == 0).sum())

            row = {
                "dataset": dataset,
                "model": model,
                "s3_total": n_total,
                "s3_predicted_positive": n_positive,
                "s3_predicted_negative": n_negative,
                "s3_positive_fraction": n_positive / n_total if n_total > 0 else None,
                "s3_mean_probability": df["prediction_probability"].mean() if "prediction_probability" in df.columns else None,
                "s3_min_probability": df["prediction_probability"].min() if "prediction_probability" in df.columns else None,
                "s3_max_probability": df["prediction_probability"].max() if "prediction_probability" in df.columns else None,
            }

            if "sequence" in df.columns:
                neg_df = df[df["predicted_label"] == 0]
                row["s3_negative_sequences"] = "; ".join(neg_df["sequence"].astype(str).tolist())
            else:
                row["s3_negative_sequences"] = ""

            rows.append(row)

    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-root", default="results/baseline/results_by_dataset")
    parser.add_argument("--predictions-root", default="data/external/s3_predictions")
    parser.add_argument("--output-dir", default="results/baseline/model_summaries")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    internal_df = summarize_internal(Path(args.results_root))
    external_df = summarize_external(Path(args.predictions_root))

    if internal_df.empty:
        print("No internal results found.")
    else:
        internal_df = internal_df.sort_values(
            by=["dataset", "test_mcc", "generalization_gap_abs"],
            ascending=[True, False, True],
        )
        internal_df.to_csv(output_dir / "internal_summary.tsv", sep="\t", index=False)
        print(f"Saved: {output_dir / 'internal_summary.tsv'}")

    if external_df.empty:
        print("No external S3 prediction files found.")
    else:
        external_df = external_df.sort_values(
            by=["dataset", "s3_predicted_positive", "s3_mean_probability"],
            ascending=[True, False, False],
        )
        external_df.to_csv(output_dir / "external_s3_summary.tsv", sep="\t", index=False)
        print(f"Saved: {output_dir / 'external_s3_summary.tsv'}")

    if not internal_df.empty and not external_df.empty:
        merged = pd.merge(
            internal_df,
            external_df,
            on=["dataset", "model"],
            how="outer",
        )

        merged["overall_rank_score"] = (
            merged["test_mcc"].fillna(0)
            + merged["s3_positive_fraction"].fillna(0)
            - merged["generalization_gap_abs"].fillna(0)
        )

        merged = merged.sort_values(
            by=["dataset", "overall_rank_score", "test_mcc", "s3_positive_fraction"],
            ascending=[True, False, False, False],
        )

        merged.to_csv(output_dir / "combined_model_dataset_summary.tsv", sep="\t", index=False)
        print(f"Saved: {output_dir / 'combined_model_dataset_summary.tsv'}")

        best_rows = []
        for dataset, grp in merged.groupby("dataset"):
            best_rows.append(grp.iloc[0])

        best_df = pd.DataFrame(best_rows)
        best_df.to_csv(output_dir / "best_by_dataset.tsv", sep="\t", index=False)
        print(f"Saved: {output_dir / 'best_by_dataset.tsv'}")

        best_overall = merged.iloc[0:1]
        best_overall.to_csv(output_dir / "best_overall.tsv", sep="\t", index=False)
        print(f"Saved: {output_dir / 'best_overall.tsv'}")

        print("\nTop combinations:")
        cols = [
            "dataset", "model", "test_mcc", "generalization_gap_abs",
            "s3_predicted_positive", "s3_total", "s3_mean_probability",
            "overall_rank_score"
        ]
        print(merged[cols].head(10).to_string(index=False))


if __name__ == "__main__":
    main()