#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import shap


def load_feature_table(path: str):
    if path.endswith(".csv"):
        df = pd.read_csv(path)
    else:
        df = pd.read_csv(path, sep="\t")

    feature_cols = [c for c in df.columns if c not in ["sequence", "label"]]
    X = df[feature_cols].copy()
    y = df["label"].copy() if "label" in df.columns else None
    return df, X, y, feature_cols


def get_sample(X: pd.DataFrame, n: int, random_state: int):
    if len(X) <= n:
        return X.copy()
    return X.sample(n=n, random_state=random_state)


def get_pipeline_parts(model):
    if hasattr(model, "named_steps"):
        steps = model.named_steps
        estimator = steps.get("clf", model)

        def transform_X(X):
            X_out = X.copy()
            if "imputer" in steps:
                X_out = steps["imputer"].transform(X_out)
            if "scaler" in steps:
                X_out = steps["scaler"].transform(X_out)
            return X_out

        return transform_X, estimator
    else:
        return lambda X: X, model


def choose_predict_fn(estimator):
    if hasattr(estimator, "predict_proba"):
        return lambda X: estimator.predict_proba(X)[:, 1]
    elif hasattr(estimator, "decision_function"):
        return lambda X: estimator.decision_function(X)
    else:
        return lambda X: estimator.predict(X)


def build_explainer(estimator, X_background_df, model_name: str, max_evals: int):
    model_name_lower = model_name.lower()

    # Faster/more exact for tree models
    if "xgboost" in model_name_lower or "random_forest" in model_name_lower:
        return shap.Explainer(estimator, X_background_df)

    # Generic fallback for non-tree models
    predict_fn = choose_predict_fn(estimator)
    return shap.Explainer(predict_fn, X_background_df, max_evals=max_evals)


def normalize_shap_values(shap_values):
    """
    Convert SHAP output into a 2D matrix: (n_samples, n_features)
    """
    values = shap_values.values

    # Common case: already 2D
    if values.ndim == 2:
        return values

    # Binary classifier or multiclass style outputs
    if values.ndim == 3:
        # Usually last dim is class index; use positive class if available
        if values.shape[2] >= 2:
            return values[:, :, 1]
        else:
            return values[:, :, 0]

    raise ValueError(f"Unsupported SHAP value shape: {values.shape}")


def make_2d_explanation(shap_values):
    """
    Returns a SHAP Explanation object usable for beeswarm/bar plots in binary classification
    """
    values = normalize_shap_values(shap_values)

    data = shap_values.data
    if isinstance(data, pd.DataFrame):
        data_arr = data.values
        feature_names = list(data.columns)
    else:
        data_arr = np.array(data)
        feature_names = list(shap_values.feature_names)

    base_values = shap_values.base_values
    base_values = np.array(base_values)

    if base_values.ndim == 2:
        if base_values.shape[1] >= 2:
            base_values = base_values[:, 1]
        else:
            base_values = base_values[:, 0]

    return shap.Explanation(
        values=values,
        base_values=base_values,
        data=data_arr,
        feature_names=feature_names,
    )


def mean_abs_shap_table(shap_values_2d, feature_names):
    values = shap_values_2d.values
    mean_abs = np.abs(values).mean(axis=0)

    df_top = pd.DataFrame({
        "feature": feature_names,
        "mean_abs_shap": mean_abs
    }).sort_values("mean_abs_shap", ascending=False)

    return df_top


def save_beeswarm_plot(shap_values_2d, out_path, max_display=20):
    plt.figure()
    shap.plots.beeswarm(shap_values_2d, max_display=max_display, show=False)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()


def save_bar_plot(shap_values_2d, out_path, max_display=20):
    plt.figure()
    shap.plots.bar(shap_values_2d, max_display=max_display, show=False)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--datasets",
        nargs="+",
        default=["pcp_aac.tsv", "pcp_aac_ctd.tsv", "pcp_aac_ctd_dpc.tsv"],
        help="Feature tables to analyze"
    )
    parser.add_argument(
        "--models-root",
        default="saved_models",
        help="Root directory containing per-dataset model subfolders"
    )
    parser.add_argument(
        "--output-root",
        default="shap_outputs_high_quality",
        help="Root directory to save SHAP outputs"
    )
    parser.add_argument(
        "--background-size",
        type=int,
        default=120,
        help="Number of samples for SHAP background set"
    )
    parser.add_argument(
        "--explain-size",
        type=int,
        default=80,
        help="Number of samples to explain"
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=25,
        help="Number of top features to save"
    )
    parser.add_argument(
        "--repeats",
        type=int,
        default=2,
        help="Number of repeated SHAP runs"
    )
    parser.add_argument(
        "--base-seed",
        type=int,
        default=42,
        help="Base random seed"
    )
    args = parser.parse_args()

    models_root = Path(args.models_root)
    output_root = Path(args.output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    summary_rows = []

    for dataset_path in args.datasets:
        dataset_tag = Path(dataset_path).stem
        dataset_models_dir = models_root / dataset_tag
        dataset_output_dir = output_root / dataset_tag
        dataset_output_dir.mkdir(parents=True, exist_ok=True)

        if not dataset_models_dir.exists():
            print(f"Skipping {dataset_tag}: no model folder found at {dataset_models_dir}")
            continue

        print(f"\n=== High-quality SHAP for dataset: {dataset_tag} ===")

        df, X_raw, y, feature_names = load_feature_table(dataset_path)
        min_required_max_evals = 2 * len(feature_names) + 1

        for model_file in sorted(dataset_models_dir.glob("*.pkl")):
            model_name = model_file.stem
            print(f"  -> {model_name}")

            model_out_dir = dataset_output_dir / model_name
            model_out_dir.mkdir(parents=True, exist_ok=True)

            try:
                model = joblib.load(model_file)
                transform_X, estimator = get_pipeline_parts(model)

                repeat_tables = []

                for repeat_idx in range(args.repeats):
                    seed = args.base_seed + repeat_idx
                    print(f"     repeat {repeat_idx+1}/{args.repeats} (seed={seed})")

                    X_background_raw = get_sample(X_raw, args.background_size, seed)
                    X_explain_raw = get_sample(X_raw, args.explain_size, seed)

                    X_background = transform_X(X_background_raw)
                    X_explain = transform_X(X_explain_raw)

                    if not isinstance(X_background, pd.DataFrame):
                        X_background_df = pd.DataFrame(X_background, columns=feature_names)
                    else:
                        X_background_df = X_background.copy()

                    if not isinstance(X_explain, pd.DataFrame):
                        X_explain_df = pd.DataFrame(X_explain, columns=feature_names)
                    else:
                        X_explain_df = X_explain.copy()

                    explainer = build_explainer(
                        estimator=estimator,
                        X_background_df=X_background_df,
                        model_name=model_name,
                        max_evals=min_required_max_evals
                    )

                    shap_values = explainer(X_explain_df)
                    shap_values_2d = make_2d_explanation(shap_values)

                    df_top = mean_abs_shap_table(shap_values_2d, feature_names)
                    df_top["repeat"] = repeat_idx + 1
                    repeat_tables.append(df_top)

                    per_repeat_path = model_out_dir / f"top_features_repeat_{repeat_idx+1}.tsv"
                    df_top.head(args.top_n).to_csv(per_repeat_path, sep="\t", index=False)

                    # Save plots only on first repeat
                    if repeat_idx == 0:
                        beeswarm_path = model_out_dir / "shap_beeswarm.png"
                        bar_path = model_out_dir / "shap_bar.png"

                        try:
                            save_beeswarm_plot(shap_values_2d, beeswarm_path, max_display=args.top_n)
                        except Exception as e:
                            print(f"       beeswarm plot failed: {e}")

                        try:
                            save_bar_plot(shap_values_2d, bar_path, max_display=args.top_n)
                        except Exception as e:
                            print(f"       bar plot failed: {e}")

                combined = pd.concat(repeat_tables, ignore_index=True)

                avg_df = (
                    combined.groupby("feature", as_index=False)["mean_abs_shap"]
                    .agg(["mean", "std"])
                    .reset_index()
                )
                avg_df.columns = ["feature", "mean_abs_shap_mean", "mean_abs_shap_std"]
                avg_df = avg_df.sort_values("mean_abs_shap_mean", ascending=False)

                avg_path = model_out_dir / "top_features_average.tsv"
                avg_df.head(args.top_n).to_csv(avg_path, sep="\t", index=False)

                for rank, (_, row) in enumerate(avg_df.head(args.top_n).iterrows(), start=1):
                    summary_rows.append({
                        "dataset": dataset_tag,
                        "model": model_name,
                        "rank": rank,
                        "feature": row["feature"],
                        "mean_abs_shap_mean": row["mean_abs_shap_mean"],
                        "mean_abs_shap_std": row["mean_abs_shap_std"],
                    })

            except Exception as e:
                print(f"     SHAP failed for {dataset_tag} / {model_name}: {e}")
                continue

    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        summary_path = output_root / "shap_top_features_summary.tsv"
        summary_df.to_csv(summary_path, sep="\t", index=False)
        print(f"\nSaved summary: {summary_path}")
    else:
        print("\nNo SHAP summaries were successfully created.")


if __name__ == "__main__":
    main()