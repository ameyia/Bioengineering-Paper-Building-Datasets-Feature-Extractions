#!/usr/bin/env python3

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from propy.CTD import CalculateCTD


AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
DIPEPTIDES = [a + b for a in AMINO_ACIDS for b in AMINO_ACIDS]
AA_SET = set(AMINO_ACIDS)


@dataclass(frozen=True)
class SequenceRecord:
    sequence: str
    label: int
    record_id: str
    source: str
    detail: str = ""
    sample_weight: float | None = None


def is_valid_sequence(sequence: str, min_length: int) -> bool:
    sequence = str(sequence).strip().upper()
    return len(sequence) >= min_length and set(sequence).issubset(AA_SET)


def longest_common_contiguous(a: str, b: str) -> str:
    max_len = min(len(a), len(b))
    for length in range(max_len, 0, -1):
        for start in range(len(a) - length + 1):
            sub = a[start:start + length]
            if sub in b:
                return sub
    return ""


def is_full_containment_overlap(a: str, b: str, min_length: int) -> tuple[bool, str, float]:
    """Return true when one model-ready peptide is fully contained in the other."""
    if len(a) < min_length or len(b) < min_length:
        return False, "", 0.0
    sub = longest_common_contiguous(a, b)
    frac = len(sub) / min(len(a), len(b)) if min(len(a), len(b)) else 0.0
    return frac >= 1.0, sub, frac * 100


def load_s1_positives(path: Path, min_length: int) -> tuple[list[SequenceRecord], list[dict]]:
    records: list[SequenceRecord] = []
    manifest: list[dict] = []
    for record in SeqIO.parse(str(path), "fasta"):
        seq = str(record.seq).strip().upper()
        row = {
            "record_id": record.id,
            "source": "s1_positive",
            "sequence": seq,
            "label": 1,
            "status": "kept",
            "reason": "",
            "overlap_with": "",
            "overlap_sequence": "",
            "overlap_pct": 0.0,
        }
        if not is_valid_sequence(seq, min_length):
            row["status"] = "removed"
            row["reason"] = f"invalid_or_shorter_than_{min_length}"
        else:
            records.append(SequenceRecord(seq, 1, record.id, "s1_positive"))
        manifest.append(row)
    return records, manifest


def load_proxy_positives(path: Path, min_length: int) -> tuple[list[SequenceRecord], list[dict]]:
    df = pd.read_csv(path, sep="\t")
    records: list[SequenceRecord] = []
    manifest: list[dict] = []

    for _, row in df.iterrows():
        seq = str(row["sequence_for_model"]).strip().upper()
        label = int(row["label"])
        record_id = str(row["name"])
        source = f"proxy:{row.get('source_tab', '')}"
        detail = str(row.get("evidence_strength", ""))
        sample_weight = row.get("sample_weight")
        out = {
            "record_id": record_id,
            "source": source,
            "sequence": seq,
            "label": label,
            "status": "kept",
            "reason": "",
            "overlap_with": "",
            "overlap_sequence": "",
            "overlap_pct": 0.0,
            "sample_weight": sample_weight,
            "evidence_strength": detail,
            "notes": row.get("notes", ""),
        }
        if label != 1:
            out["status"] = "removed"
            out["reason"] = "proxy_label_not_positive"
        elif not is_valid_sequence(seq, min_length):
            out["status"] = "removed"
            out["reason"] = f"invalid_or_shorter_than_{min_length}"
        else:
            records.append(SequenceRecord(seq, 1, record_id, source, detail, sample_weight))
        manifest.append(out)
    return records, manifest


def load_negatives(path: Path, min_length: int) -> tuple[list[SequenceRecord], list[dict]]:
    df = pd.read_csv(path, sep="\t", header=None)
    records: list[SequenceRecord] = []
    manifest: list[dict] = []
    for idx, row in df.iterrows():
        record_id = str(row.iloc[0])
        seq = str(row.iloc[1]).strip().upper()
        detail = str(row.iloc[3]) if len(row) > 3 else ""
        out = {
            "record_id": record_id,
            "source": "negative",
            "sequence": seq,
            "label": 0,
            "status": "kept",
            "reason": "",
            "overlap_with": "",
            "overlap_sequence": "",
            "overlap_pct": 0.0,
        }
        if not is_valid_sequence(seq, min_length):
            out["status"] = "removed"
            out["reason"] = f"invalid_or_shorter_than_{min_length}"
        else:
            records.append(SequenceRecord(seq, 0, record_id or f"negative_{idx + 1}", "negative", detail))
        manifest.append(out)
    return records, manifest


def dedupe_and_remove_overlaps(records: list[SequenceRecord], min_length: int) -> tuple[list[SequenceRecord], list[dict]]:
    kept: list[SequenceRecord] = []
    audit: list[dict] = []
    seen_exact: dict[str, SequenceRecord] = {}

    for record in records:
        status = "kept"
        reason = ""
        overlap_with = ""
        overlap_sequence = ""
        overlap_pct = 0.0

        if record.sequence in seen_exact:
            prior = seen_exact[record.sequence]
            status = "removed"
            reason = "duplicate_positive_sequence"
            overlap_with = prior.record_id
            overlap_sequence = record.sequence
            overlap_pct = 100.0

        if status == "kept":
            for prior in kept:
                overlaps, sub, pct = is_full_containment_overlap(record.sequence, prior.sequence, min_length)
                if overlaps:
                    # Preserve the first curated S1/proxy representative encountered.
                    status = "removed"
                    reason = "contained_in_kept_positive_sequence"
                    overlap_with = prior.record_id
                    overlap_sequence = sub
                    overlap_pct = pct
                    break

        audit.append({
            "record_id": record.record_id,
            "source": record.source,
            "sequence": record.sequence,
            "label": record.label,
            "status": status,
            "reason": reason,
            "overlap_with": overlap_with,
            "overlap_sequence": overlap_sequence,
            "overlap_pct": overlap_pct,
            "sample_weight": record.sample_weight,
            "detail": record.detail,
        })

        if status == "kept":
            kept.append(record)
            seen_exact[record.sequence] = record

    return kept, audit


def remove_negative_overlaps(
    negatives: list[SequenceRecord],
    positives: list[SequenceRecord],
    min_length: int,
) -> tuple[list[SequenceRecord], list[dict]]:
    kept: list[SequenceRecord] = []
    audit: list[dict] = []
    seen_exact: set[str] = set()

    for record in negatives:
        status = "kept"
        reason = ""
        overlap_with = ""
        overlap_sequence = ""
        overlap_pct = 0.0

        if record.sequence in seen_exact:
            status = "removed"
            reason = "duplicate_negative_sequence"
            overlap_sequence = record.sequence
            overlap_pct = 100.0

        if status == "kept":
            for positive in positives:
                overlaps, sub, pct = is_full_containment_overlap(record.sequence, positive.sequence, min_length)
                if overlaps:
                    status = "removed"
                    reason = "overlaps_positive_sequence"
                    overlap_with = positive.record_id
                    overlap_sequence = sub
                    overlap_pct = pct
                    break

        audit.append({
            "record_id": record.record_id,
            "source": record.source,
            "sequence": record.sequence,
            "label": record.label,
            "status": status,
            "reason": reason,
            "overlap_with": overlap_with,
            "overlap_sequence": overlap_sequence,
            "overlap_pct": overlap_pct,
            "sample_weight": record.sample_weight,
            "detail": record.detail,
        })

        if status == "kept":
            kept.append(record)
            seen_exact.add(record.sequence)

    return kept, audit


def featurize_sequences(records: list[SequenceRecord]) -> pd.DataFrame:
    features = []

    for record in records:
        seq = record.sequence
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
                "label": record.label,
            }

            for aa in AMINO_ACIDS:
                feature_dict[aa] = (seq.count(aa) / length) * 100

            try:
                feature_dict.update(CalculateCTD(seq))
            except Exception:
                print(f"CTD failed for: {seq}")

            dpc = dict.fromkeys(DIPEPTIDES, 0.0)
            total_dipeptides = len(seq) - 1
            if total_dipeptides > 0:
                for i in range(total_dipeptides):
                    dp = seq[i:i + 2]
                    if dp in dpc:
                        dpc[dp] += 1
                for dp in dpc:
                    dpc[dp] = (dpc[dp] / total_dipeptides) * 100
            feature_dict.update(dpc)
            features.append(feature_dict)
        except Exception as exc:
            print(f"Error with sequence {seq}: {exc}")

    return pd.DataFrame(features)


def write_fasta(records: list[SequenceRecord], path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for idx, record in enumerate(records, 1):
            handle.write(f">{record.record_id}|{record.source}|positive_{idx}\n{record.sequence}\n")


def write_negatives(records: list[SequenceRecord], path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for record in records:
            length_bin = record.detail if record.detail else f"{len(record.sequence)}aa"
            handle.write(f"{record.record_id}\t{record.sequence}\t{len(record.sequence)}\t{length_bin}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s1-positive-fasta", default="TableS1PositiveProteins.fasta")
    parser.add_argument("--negative-tsv", default="final_negatives.tsv")
    parser.add_argument("--proxy-ready-tsv", default="/Users/ameliamccrory/Downloads/proxy_peptides_model_ready.tsv")
    parser.add_argument("--min-length", type=int, default=5)
    parser.add_argument("--output-prefix", default="s1_proxy_no_overlap")
    args = parser.parse_args()

    prefix = args.output_prefix
    s1_records, s1_manifest = load_s1_positives(Path(args.s1_positive_fasta), args.min_length)
    proxy_records, proxy_manifest = load_proxy_positives(Path(args.proxy_ready_tsv), args.min_length)
    negative_records, negative_manifest = load_negatives(Path(args.negative_tsv), args.min_length)

    positives, positive_audit = dedupe_and_remove_overlaps(s1_records + proxy_records, args.min_length)
    negatives, negative_audit = remove_negative_overlaps(negative_records, positives, args.min_length)

    all_records = positives + negatives
    df_full = featurize_sequences(all_records)
    base_columns = [
        "sequence",
        "length",
        "mw",
        "aromaticity",
        "pI",
        "instability",
        "label",
    ] + AMINO_ACIDS
    ctd_columns = [col for col in df_full.columns if "_" in col]

    df_base = df_full[base_columns]
    df_ctd = df_full[base_columns + ctd_columns]
    df_dpc = df_full

    write_fasta(positives, Path(f"{prefix}_positives.fasta"))
    write_negatives(negatives, Path(f"{prefix}_negatives.tsv"))
    df_base.to_csv(f"{prefix}_pcp_aac.tsv", sep="\t", index=False)
    df_ctd.to_csv(f"{prefix}_pcp_aac_ctd.tsv", sep="\t", index=False)
    df_dpc.to_csv(f"{prefix}_pcp_aac_ctd_dpc.tsv", sep="\t", index=False)

    manifest_rows = []
    manifest_rows.extend(s1_manifest)
    manifest_rows.extend(proxy_manifest)
    manifest_rows.extend(negative_manifest)
    pd.DataFrame(manifest_rows).to_csv(f"{prefix}_input_manifest.tsv", sep="\t", index=False)
    pd.DataFrame(positive_audit + negative_audit).to_csv(f"{prefix}_overlap_audit.tsv", sep="\t", index=False)

    kept_proxy = sum(1 for record in positives if record.source.startswith("proxy:"))
    print("Saved proxy-integrated no-overlap dataset")
    print(f"  positives kept: {len(positives)} ({kept_proxy} proxy)")
    print(f"  negatives kept: {len(negatives)}")
    print(f"  total feature rows: {len(df_full)}")
    print(f"  wrote: {prefix}_pcp_aac.tsv")
    print(f"  wrote: {prefix}_pcp_aac_ctd.tsv")
    print(f"  wrote: {prefix}_pcp_aac_ctd_dpc.tsv")


if __name__ == "__main__":
    main()
