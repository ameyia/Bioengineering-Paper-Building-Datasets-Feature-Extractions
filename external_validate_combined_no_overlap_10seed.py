#!/usr/bin/env python3

from __future__ import annotations

import argparse
import math
import os
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import numpy as np
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from propy.CTD import CalculateCTD
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, f1_score, matthews_corrcoef, precision_score, recall_score
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


OVERLAP_THRESHOLD = 0.9845  # 64/65 is reported as 98.5% when rounded.
MIN_OVERLAP_LENGTH = 8
AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
DIPEPTIDES = [a + b for a in AMINO_ACIDS for b in AMINO_ACIDS]


@dataclass(frozen=True)
class ExternalPeptide:
    name: str
    source: str
    sequence: str
    label: int = 1


EXTERNAL_CANDIDATES = [
    ExternalPeptide("avathrin", "s3", "SGGHQTAVPKISKQGLGGDFEEIPSDEIIE"),
    ExternalPeptide("s3_variant_1", "s3", "SDEAVRAIPKMYSTAPPGDFETIPDDAIEEREMKAR"),
    ExternalPeptide("s3_variant_2", "s3", "SDEAVRAIPKMYSTAPPGDFETIPDDAIEER"),
    ExternalPeptide("s3_variant_3", "s3", "SDEAVRAIPKMYSTAPPGDFETIPDDAIEE"),
    ExternalPeptide("ultravariegin", "s3", "SDEAVRAIPKMYSTAPPGDFEEIPDDAIEE"),
    ExternalPeptide("uv003", "s3", "SDQGDVAIPKMYSTAPPGDFEEIPDDAIEE"),
    ExternalPeptide("uv004", "s3", "SDEAVRAEPKMHKTAPPGDFEEIPDDAIEE"),
    ExternalPeptide("uv005", "s3", "SDEAVRAIPKMYSTAPPGDFEEIPEEYLDDES"),
    ExternalPeptide("uv012", "s3", "SDEAVRAIPKMYSTAPPGDFEEIPDDEIEE"),
    ExternalPeptide("uv013", "s3", "SDEAVRAIPKMYSQAPPGDFEEIPDDAIEE"),
    ExternalPeptide("uv011", "s3", "MYSTAPPGDFEEIPDDAIEE"),
    ExternalPeptide("lepirudin", "s5", "LVYTDCTESGQNLCLCEGSNVCGQGNKCILGSDGEKNQCVTGEGTPKPQSHNDGDFEEIPEEYLQ"),
    ExternalPeptide("desirudin", "s5", "VVYTDCTESGQNLCLCEGSNVCGQGNKCILGSDGEKNQCVTGEGTPKPQSHNDGDFEEIPEEYLQ"),
    ExternalPeptide("bivalirudin", "s5", "FPRPGGGGNGDFEEIPEEYL"),
    ExternalPeptide("new_spfrvvfvkp", "new", "SPFRVVFVKP"),
    ExternalPeptide("new_gcs", "new", "GCSGKGARCAPSKCCSGLSCGRHGGNMYKSCEWNWKTG"),
    ExternalPeptide("new_adc", "new", "ADCGGKTCSGGQVCSDGVCVCTKLRCRLLCRNGFLKDENGCEYPCTCA"),
    ExternalPeptide("new_vis", "new", "VISRTQSNVQAAWGQVGGHAADYSAVAIER"),
]


def dataset_tag(path: str) -> str:
    return Path(path).stem


def load_table(path: str):
    sep = "," if path.endswith(".csv") else "\t"
    df = pd.read_csv(path, sep=sep)
    if "label" not in df.columns:
        raise ValueError(f"{path} does not contain a label column")

    sequences = df["sequence"].astype(str).str.upper().tolist() if "sequence" in df.columns else []
    y = df["label"].astype(int)
    drop_cols = ["label"]
    if "sequence" in df.columns:
        drop_cols.append("sequence")
    x = df.drop(columns=drop_cols)
    x = x.apply(pd.to_numeric, errors="coerce")
    x = x.replace([np.inf, -np.inf], np.nan)
    x = x.fillna(x.median(numeric_only=True))
    return df, x, y, sequences, list(x.columns)


def split_data(x, y, seed):
    x_train, x_temp, y_train, y_temp = train_test_split(
        x, y, test_size=0.40, stratify=y, random_state=seed
    )
    x_val, x_test, y_val, y_test = train_test_split(
        x_temp, y_temp, test_size=0.50, stratify=y_temp, random_state=seed
    )
    return x_train, x_val, x_test, y_train, y_val, y_test


def get_models(seed):
    models = {
        "svm_linear": Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            ("model", SVC(kernel="linear", probability=True, class_weight="balanced", random_state=seed)),
        ]),
        "svm_rbf": Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            ("model", SVC(kernel="rbf", C=10, gamma="scale", probability=True, class_weight="balanced", random_state=seed)),
        ]),
        "logistic_regression": Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            ("model", LogisticRegression(max_iter=5000, class_weight="balanced", random_state=seed)),
        ]),
        "random_forest": Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("model", RandomForestClassifier(
                n_estimators=500,
                min_samples_leaf=2,
                class_weight="balanced",
                random_state=seed,
                n_jobs=-1,
            )),
        ]),
        "knn": Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            ("model", KNeighborsClassifier(n_neighbors=5)),
        ]),
    }

    if HAS_XGB:
        models["xgboost"] = Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("model", XGBClassifier(
                n_estimators=300,
                max_depth=3,
                learning_rate=0.03,
                subsample=0.8,
                colsample_bytree=0.8,
                eval_metric="logloss",
                random_state=seed,
                n_jobs=-1,
            )),
        ])

    return models


def longest_common_contiguous(a: str, b: str) -> str:
    max_len = min(len(a), len(b))
    for length in range(max_len, 0, -1):
        for start in range(len(a) - length + 1):
            sub = a[start:start + length]
            if sub in b:
                return sub
    return ""


def overlap_fraction(a: str, b: str, min_overlap_length: int, denominator: str) -> tuple[float, str]:
    if not a or not b:
        return 0.0, ""
    sub = longest_common_contiguous(a, b)
    if len(sub) < min_overlap_length:
        return 0.0, sub
    if denominator == "max":
        denom = max(len(a), len(b))
    elif denominator == "min":
        denom = min(len(a), len(b))
    else:
        raise ValueError("denominator must be 'min' or 'max'")
    return len(sub) / denom, sub


def clean_external_set(training_sequences: list[str], threshold: float, min_overlap_length: int):
    kept = []
    manifest_rows = []
    seen_sequences = {}

    for peptide in EXTERNAL_CANDIDATES:
        sequence = peptide.sequence.upper()
        status = "kept"
        reason = ""
        overlap_with = ""
        overlap_sequence = ""
        overlap_pct = 0.0

        if sequence in seen_sequences:
            status = "removed"
            reason = "duplicate_external_sequence"
            overlap_with = seen_sequences[sequence]
            overlap_sequence = sequence
            overlap_pct = 100.0

        if status == "kept":
            best_training_overlap = (0.0, "", "")
            for train_seq in training_sequences:
                frac, sub = overlap_fraction(sequence, train_seq, min_overlap_length, "max")
                if frac > best_training_overlap[0]:
                    best_training_overlap = (frac, sub, train_seq)
            if best_training_overlap[0] >= threshold:
                status = "removed"
                reason = "overlaps_training_set_at_or_above_threshold"
                overlap_with = best_training_overlap[2]
                overlap_sequence = best_training_overlap[1]
                overlap_pct = best_training_overlap[0] * 100

        if status == "kept":
            for prior in kept:
                frac, sub = overlap_fraction(sequence, prior.sequence, min_overlap_length, "min")
                if frac >= threshold:
                    status = "removed"
                    reason = "overlaps_kept_external_peptide_at_or_above_threshold"
                    overlap_with = prior.name
                    overlap_sequence = sub
                    overlap_pct = frac * 100
                    break

        manifest_rows.append({
            "name": peptide.name,
            "source": peptide.source,
            "sequence": sequence,
            "label": peptide.label,
            "status": status,
            "reason": reason,
            "overlap_with": overlap_with,
            "overlap_sequence": overlap_sequence,
            "overlap_pct": overlap_pct,
        })

        seen_sequences[sequence] = peptide.name
        if status == "kept":
            kept.append(peptide)

    return kept, pd.DataFrame(manifest_rows)


def featurize_external(peptides: list[ExternalPeptide], feature_cols: list[str]) -> pd.DataFrame:
    rows = []
    use_ctd = any("_" in col for col in feature_cols)
    use_dpc = any(len(col) == 2 and col in DIPEPTIDES for col in feature_cols)

    for peptide in peptides:
        seq = peptide.sequence.upper()
        analysis = ProteinAnalysis(seq)
        length = len(seq)
        row = {
            "name": peptide.name,
            "source": peptide.source,
            "sequence": seq,
            "label": peptide.label,
            "length": length,
            "mw": analysis.molecular_weight(),
            "aromaticity": analysis.aromaticity(),
            "pI": analysis.isoelectric_point(),
            "instability": analysis.instability_index(),
        }

        for aa in AMINO_ACIDS:
            row[aa] = (seq.count(aa) / length) * 100

        if use_ctd:
            row.update(CalculateCTD(seq))

        if use_dpc:
            dpc = dict.fromkeys(DIPEPTIDES, 0.0)
            total = len(seq) - 1
            if total > 0:
                for i in range(total):
                    dp = seq[i:i + 2]
                    if dp in dpc:
                        dpc[dp] += 1
                for dp in dpc:
                    dpc[dp] = (dpc[dp] / total) * 100
            row.update(dpc)

        rows.append(row)

    df = pd.DataFrame(rows)
    for col in feature_cols:
        if col not in df.columns:
            df[col] = np.nan
    return df


def ci95(series: pd.Series) -> float:
    values = series.dropna()
    if len(values) <= 1:
        return 0.0
    return 1.96 * values.std(ddof=1) / math.sqrt(len(values))


def summarize_with_ci(df: pd.DataFrame, group_cols: list[str], metric_cols: list[str]) -> pd.DataFrame:
    rows = []
    for keys, group in df.groupby(group_cols, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        row = dict(zip(group_cols, keys))
        row["n_seeds"] = group["seed"].nunique() if "seed" in group.columns else len(group)
        for metric in metric_cols:
            row[f"mean_{metric}"] = group[metric].mean()
            row[f"std_{metric}"] = group[metric].std(ddof=1)
            row[f"ci95_{metric}"] = ci95(group[metric])
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--tables",
        nargs="+",
        default=["pcp_aac.tsv", "pcp_aac_ctd.tsv", "pcp_aac_ctd_dpc.tsv"],
    )
    parser.add_argument("--seeds", type=int, default=10)
    parser.add_argument("--overlap-threshold", type=float, default=OVERLAP_THRESHOLD)
    parser.add_argument("--min-overlap-length", type=int, default=MIN_OVERLAP_LENGTH)
    parser.add_argument("--outdir", default="external_validation_no_overlap")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    first_df, _, _, training_sequences, _ = load_table(args.tables[0])
    kept_external, manifest_df = clean_external_set(
        training_sequences,
        args.overlap_threshold,
        args.min_overlap_length,
    )
    manifest_df.to_csv(outdir / "combined_external_validation_manifest.tsv", sep="\t", index=False)

    kept_df = manifest_df[manifest_df["status"] == "kept"].copy()
    kept_df.to_csv(outdir / "combined_external_validation_set.tsv", sep="\t", index=False)
    print(f"Kept {len(kept_external)} of {len(EXTERNAL_CANDIDATES)} external peptides")
    print(f"Saved: {outdir / 'combined_external_validation_manifest.tsv'}")
    print(f"Saved: {outdir / 'combined_external_validation_set.tsv'}")

    internal_rows = []
    prediction_rows = []

    for table in args.tables:
        tag = dataset_tag(table)
        print(f"\n=== {tag} ===")
        _, x, y, _, feature_cols = load_table(table)
        external_df = featurize_external(kept_external, feature_cols)
        external_x = external_df[feature_cols].apply(pd.to_numeric, errors="coerce")

        for seed in range(args.seeds):
            print(f"  seed {seed}")
            x_train, x_val, x_test, y_train, y_val, y_test = split_data(x, y, seed)

            for model_name, model in get_models(seed).items():
                print(f"    {model_name}")
                model.fit(x_train, y_train)
                val_pred = model.predict(x_val)
                test_pred = model.predict(x_test)

                internal_rows.append({
                    "dataset": tag,
                    "seed": seed,
                    "model": model_name,
                    "val_mcc": matthews_corrcoef(y_val, val_pred),
                    "test_mcc": matthews_corrcoef(y_test, test_pred),
                    "test_accuracy": accuracy_score(y_test, test_pred),
                    "test_f1": f1_score(y_test, test_pred, zero_division=0),
                    "test_precision": precision_score(y_test, test_pred, zero_division=0),
                    "test_recall": recall_score(y_test, test_pred, zero_division=0),
                })

                probs = model.predict_proba(external_x)[:, 1]
                preds = model.predict(external_x)
                for row, prob, pred in zip(external_df.to_dict("records"), probs, preds):
                    prediction_rows.append({
                        "dataset": tag,
                        "seed": seed,
                        "model": model_name,
                        "name": row["name"],
                        "source": row["source"],
                        "sequence": row["sequence"],
                        "label": row["label"],
                        "predicted_label": int(pred),
                        "prediction_probability": float(prob),
                    })

    internal_df = pd.DataFrame(internal_rows)
    predictions_df = pd.DataFrame(prediction_rows)

    internal_df.to_csv(outdir / "internal_metrics_by_seed.tsv", sep="\t", index=False)
    predictions_df.to_csv(outdir / "combined_external_predictions_by_seed.tsv", sep="\t", index=False)

    external_seed_df = (
        predictions_df.groupby(["dataset", "seed", "model"], as_index=False)
        .agg(
            external_total=("sequence", "count"),
            external_predicted_positive=("predicted_label", "sum"),
            external_positive_fraction=("predicted_label", "mean"),
            external_mean_probability=("prediction_probability", "mean"),
            external_min_probability=("prediction_probability", "min"),
            external_max_probability=("prediction_probability", "max"),
        )
    )
    external_seed_df.to_csv(outdir / "combined_external_summary_by_seed.tsv", sep="\t", index=False)

    internal_summary = summarize_with_ci(
        internal_df,
        ["dataset", "model"],
        ["val_mcc", "test_mcc", "test_accuracy", "test_f1", "test_precision", "test_recall"],
    ).sort_values(["dataset", "mean_test_mcc"], ascending=[True, False])
    internal_summary.to_csv(outdir / "internal_metrics_summary_ci.tsv", sep="\t", index=False)

    external_summary = summarize_with_ci(
        external_seed_df,
        ["dataset", "model"],
        ["external_positive_fraction", "external_mean_probability", "external_min_probability", "external_max_probability"],
    ).sort_values(["dataset", "mean_external_positive_fraction", "mean_external_mean_probability"], ascending=[True, False, False])
    external_summary.to_csv(outdir / "combined_external_summary_ci.tsv", sep="\t", index=False)

    peptide_summary = summarize_with_ci(
        predictions_df,
        ["dataset", "model", "name", "source", "sequence"],
        ["predicted_label", "prediction_probability"],
    ).sort_values(["dataset", "model", "name"])
    peptide_summary.to_csv(outdir / "combined_external_peptide_summary_ci.tsv", sep="\t", index=False)

    print("\nSaved:")
    for filename in [
        "internal_metrics_by_seed.tsv",
        "internal_metrics_summary_ci.tsv",
        "combined_external_predictions_by_seed.tsv",
        "combined_external_summary_by_seed.tsv",
        "combined_external_summary_ci.tsv",
        "combined_external_peptide_summary_ci.tsv",
    ]:
        print(f"  {outdir / filename}")


if __name__ == "__main__":
    main()
