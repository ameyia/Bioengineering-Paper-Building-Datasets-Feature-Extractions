#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def load_average_tables(shap_root: Path) -> pd.DataFrame:
    rows = []

    if not shap_root.exists():
        return pd.DataFrame()

    for dataset_dir in sorted(shap_root.iterdir()):
        if not dataset_dir.is_dir():
            continue

        dataset = dataset_dir.name

        for model_dir in sorted(dataset_dir.iterdir()):
            if not model_dir.is_dir():
                continue

            model = model_dir.name
            avg_path = model_dir / "top_features_average.tsv"

            if not avg_path.exists():
                continue

            df = pd.read_csv(avg_path, sep="\t")
            df["dataset"] = dataset
            df["model"] = model
            rows.append(df)

    if not rows:
        return pd.DataFrame()

    return pd.concat(rows, ignore_index=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--shap-root",
        default="shap_outputs_high_quality",
        help="Root directory containing SHAP outputs"
    )
    parser.add_argument(
        "--output-dir",
        default="shap_summary_outputs",
        help="Directory to save SHAP summary files"
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=20,
        help="Top N features per dataset/model to keep in filtered summaries"
    )
    args = parser.parse_args()

    shap_root = Path(args.shap_root)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df = load_average_tables(shap_root)

    if df.empty:
        print("No SHAP average tables found.")
        return

    # Full merged table
    df.to_csv(output_dir / "all_shap_features.tsv", sep="\t", index=False)
    print(f"Saved: {output_dir / 'all_shap_features.tsv'}")

    # Rank within each dataset-model
    df_ranked = df.copy()
    df_ranked["rank_within_model"] = (
        df_ranked.groupby(["dataset", "model"])["mean_abs_shap_mean"]
        .rank(method="first", ascending=False)
        .astype(int)
    )
    df_ranked.to_csv(output_dir / "all_shap_features_ranked.tsv", sep="\t", index=False)
    print(f"Saved: {output_dir / 'all_shap_features_ranked.tsv'}")

    # Top N per dataset-model
    top_df = df_ranked[df_ranked["rank_within_model"] <= args.top_n].copy()
    top_df.to_csv(output_dir / f"top_{args.top_n}_features_per_model.tsv", sep="\t", index=False)
    print(f"Saved: {output_dir / f'top_{args.top_n}_features_per_model.tsv'}")

    # Best features by dataset/model
    best_rows = []
    for (dataset, model), grp in df_ranked.groupby(["dataset", "model"]):
        best_rows.append(grp.sort_values("mean_abs_shap_mean", ascending=False).iloc[0])

    best_df = pd.DataFrame(best_rows)
    best_df.to_csv(output_dir / "best_feature_per_model.tsv", sep="\t", index=False)
    print(f"Saved: {output_dir / 'best_feature_per_model.tsv'}")

    # Feature frequency across all dataset-model combos
    freq_df = (
        top_df.groupby("feature")
        .agg(
            count_in_top_n=("feature", "count"),
            mean_importance=("mean_abs_shap_mean", "mean"),
            std_importance=("mean_abs_shap_mean", "std"),
        )
        .reset_index()
        .sort_values(["count_in_top_n", "mean_importance"], ascending=[False, False])
    )
    freq_df.to_csv(output_dir / "feature_frequency_across_models.tsv", sep="\t", index=False)
    print(f"Saved: {output_dir / 'feature_frequency_across_models.tsv'}")

    # Cross-dataset/model pivot for top features
    pivot_df = top_df.pivot_table(
        index="feature",
        columns=["dataset", "model"],
        values="mean_abs_shap_mean",
        aggfunc="mean"
    )
    pivot_df.to_csv(output_dir / "feature_importance_pivot.tsv", sep="\t")
    print(f"Saved: {output_dir / 'feature_importance_pivot.tsv'}")

    # Top features by dataset only (averaged across models)
    dataset_avg = (
        top_df.groupby(["dataset", "feature"])
        .agg(
            mean_importance=("mean_abs_shap_mean", "mean"),
            mean_std=("mean_abs_shap_std", "mean"),
            appearances=("feature", "count"),
        )
        .reset_index()
    )

    dataset_top_rows = []
    for dataset, grp in dataset_avg.groupby("dataset"):
        grp = grp.sort_values(["mean_importance", "appearances"], ascending=[False, False])
        grp["dataset_rank"] = range(1, len(grp) + 1)
        dataset_top_rows.append(grp)

    dataset_top_df = pd.concat(dataset_top_rows, ignore_index=True)
    dataset_top_df.to_csv(output_dir / "top_features_by_dataset.tsv", sep="\t", index=False)
    print(f"Saved: {output_dir / 'top_features_by_dataset.tsv'}")

    # Top features by model only (averaged across datasets)
    model_avg = (
        top_df.groupby(["model", "feature"])
        .agg(
            mean_importance=("mean_abs_shap_mean", "mean"),
            mean_std=("mean_abs_shap_std", "mean"),
            appearances=("feature", "count"),
        )
        .reset_index()
    )

    model_top_rows = []
    for model, grp in model_avg.groupby("model"):
        grp = grp.sort_values(["mean_importance", "appearances"], ascending=[False, False])
        grp["model_rank"] = range(1, len(grp) + 1)
        model_top_rows.append(grp)

    model_top_df = pd.concat(model_top_rows, ignore_index=True)
    model_top_df.to_csv(output_dir / "top_features_by_model.tsv", sep="\t", index=False)
    print(f"Saved: {output_dir / 'top_features_by_model.tsv'}")

    # Print a quick preview
    print("\nMost recurring important features across models:")
    print(freq_df.head(15).to_string(index=False))


if __name__ == "__main__":
    main()