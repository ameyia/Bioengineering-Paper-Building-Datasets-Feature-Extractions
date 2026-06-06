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

            rows.append({
                "dataset": dataset,
                "model": model_name,
                "cv_mcc": r.get("cv_mcc"),
                "test_mcc": r.get("test_mcc"),
                "test_accuracy": r.get("test_accuracy"),
                "test_f1": r.get("test_f1"),
                "generalization_gap_abs": r.get("generalization_gap_abs"),
            })

    return pd.DataFrame(rows)


def summarize_external_for_dataset(predictions_root: Path, dataset_suffix: str) -> pd.DataFrame:
    rows = []

    if not predictions_root.exists():
        return pd.DataFrame()

    for dataset_dir in sorted(predictions_root.iterdir()):
        if not dataset_dir.is_dir():
            continue

        training_dataset = dataset_dir.name

        for tsv_path in sorted(dataset_dir.glob(f"*_{dataset_suffix}_predictions.tsv")):
            model = tsv_path.name.replace(f"_{dataset_suffix}_predictions.tsv", "")
            df = pd.read_csv(tsv_path, sep="\t")

            if "predicted_label" not in df.columns:
                continue

            n_total = len(df)
            n_positive = int(df["predicted_label"].sum())
            n_negative = int((df["predicted_label"] == 0).sum())

            row = {
                "dataset": training_dataset,
                "model": model,
                f"{dataset_suffix}_total": n_total,
                f"{dataset_suffix}_predicted_positive": n_positive,
                f"{dataset_suffix}_predicted_negative": n_negative,
                f"{dataset_suffix}_positive_fraction": n_positive / n_total if n_total > 0 else None,
                f"{dataset_suffix}_mean_probability": df["prediction_probability"].mean() if "prediction_probability" in df.columns else None,
                f"{dataset_suffix}_min_probability": df["prediction_probability"].min() if "prediction_probability" in df.columns else None,
                f"{dataset_suffix}_max_probability": df["prediction_probability"].max() if "prediction_probability" in df.columns else None,
            }

            if "sequence" in df.columns:
                neg_df = df[df["predicted_label"] == 0]
                row[f"{dataset_suffix}_negative_sequences"] = "; ".join(neg_df["sequence"].astype(str).tolist())
            else:
                row[f"{dataset_suffix}_negative_sequences"] = ""

            rows.append(row)

    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-root", default="results_by_dataset")
    parser.add_argument("--predictions-root", default="external_predictions")
    parser.add_argument("--output-dir", default="summary_outputs")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    internal_df = summarize_internal(Path(args.results_root))
    s3_df = summarize_external_for_dataset(Path(args.predictions_root), "s3")
    s5_df = summarize_external_for_dataset(Path(args.predictions_root), "s5")

    if not internal_df.empty:
        internal_df.to_csv(output_dir / "internal_summary.tsv", sep="\t", index=False)
        print(f"Saved: {output_dir / 'internal_summary.tsv'}")

    if not s3_df.empty:
        s3_df.to_csv(output_dir / "external_s3_summary.tsv", sep="\t", index=False)
        print(f"Saved: {output_dir / 'external_s3_summary.tsv'}")

    if not s5_df.empty:
        s5_df.to_csv(output_dir / "external_s5_summary.tsv", sep="\t", index=False)
        print(f"Saved: {output_dir / 'external_s5_summary.tsv'}")

    merged = internal_df.copy()

    if not s3_df.empty:
        merged = pd.merge(merged, s3_df, on=["dataset", "model"], how="outer")

    if not s5_df.empty:
        merged = pd.merge(merged, s5_df, on=["dataset", "model"], how="outer")

    if not merged.empty:
        merged["overall_rank_score"] = (
            merged["test_mcc"].fillna(0)
            + merged.get("s3_positive_fraction", pd.Series(0, index=merged.index)).fillna(0)
            + merged.get("s5_positive_fraction", pd.Series(0, index=merged.index)).fillna(0)
            - merged["generalization_gap_abs"].fillna(0)
        )

        merged = merged.sort_values(
            by=["dataset", "overall_rank_score", "test_mcc"],
            ascending=[True, False, False],
        )

        merged.to_csv(output_dir / "combined_model_dataset_summary.tsv", sep="\t", index=False)
        print(f"Saved: {output_dir / 'combined_model_dataset_summary.tsv'}")

        best_rows = []
        for dataset, grp in merged.groupby("dataset"):
            best_rows.append(grp.iloc[0])

        best_df = pd.DataFrame(best_rows)
        best_df.to_csv(output_dir / "best_by_dataset.tsv", sep="\t", index=False)
        print(f"Saved: {output_dir / 'best_by_dataset.tsv'}")

        merged.iloc[0:1].to_csv(output_dir / "best_overall.tsv", sep="\t", index=False)
        print(f"Saved: {output_dir / 'best_overall.tsv'}")

        print("\nTop combinations:")
        cols = [
            "dataset", "model", "test_mcc", "generalization_gap_abs",
            "s3_predicted_positive", "s3_total",
            "s5_predicted_positive", "s5_total",
            "overall_rank_score"
        ]
        cols = [c for c in cols if c in merged.columns]
        print(merged[cols].head(10).to_string(index=False))


if __name__ == "__main__":
    main()