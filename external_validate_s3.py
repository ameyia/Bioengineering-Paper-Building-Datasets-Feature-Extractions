#!/usr/bin/env python3

import argparse
import joblib
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from propy.CTD import CalculateCTD

amino_acids = list("ACDEFGHIKLMNPQRSTVWY")

external_sequences_s3 = [
    ("SGGHQTAVPKISKQGLGGDFEEIPSDEIIE", 1),   # avathrin
    ("SDEAVRAIPKMYSTAPPGDFETIPDDAIEEREMKAR", 1),  # DAA34688.1 repeat 1
    ("SDEAVRAIPKMYSTAPPGDFETIPDDAIEER", 1),       # DAA34688.1 repeat 1B
    ("SDEAVRAIPKMYSTAPPGDFETIPDDAIEE", 1),        # DAA34688.1 repeat 1C
    ("SDEAVRAIPKMYSTAPPGDFEEIPDDAIEE", 1),        # ultravariegin
    ("SDQGDVAIPKMYSTAPPGDFEEIPDDAIEE", 1),        # UV003
    ("SDEAVRAEPKMHKTAPPGDFEEIPDDAIEE", 1),        # UV004
    ("SDEAVRAIPKMYSTAPPGDFEEIPEEYLDDES", 1),      # UV005
    ("SDEAVRAIPKMYSTAPPGDFEEIPDDEIEE", 1),        # UV012
    ("SDEAVRAIPKMYSQAPPGDFEEIPDDAIEE", 1),        # UV013
    ("MYSTAPPGDFEEIPDDAIEE", 1),                  # UV011
]


def featurize_external(sequence_label_pairs):
    rows = []

    for seq, label in sequence_label_pairs:
        seq = seq.strip().upper()
        if len(seq) == 0 or "X" in seq:
            continue

        try:
            analysis = ProteinAnalysis(seq)
            length = len(seq)

            aac = {}
            for aa in amino_acids:
                aac[aa] = (seq.count(aa) / length) * 100

            feature_dict = {
                "sequence": seq,
                "length": length,
                "mw": analysis.molecular_weight(),
                "aromaticity": analysis.aromaticity(),
                "pI": analysis.isoelectric_point(),
                "instability": analysis.instability_index(),
                "label": label,
            }

            feature_dict.update(aac)

            try:
                ctd = CalculateCTD(seq)
                feature_dict.update(ctd)
            except Exception:
                print(f"CTD failed for external sequence: {seq}")

            rows.append(feature_dict)

        except Exception as e:
            print(f"Error with external sequence: {seq}")
            print(e)

    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", default="best_model.pkl")
    parser.add_argument("--training-feature-table", default="pcp_aac_ctd.tsv")
    parser.add_argument("--output", default="s3_external_predictions.tsv")
    args = parser.parse_args()

    model = joblib.load(args.model)

    df_external = featurize_external(external_sequences_s3)

    # match training columns exactly
    if args.training_feature_table.endswith(".csv"):
        df_train = pd.read_csv(args.training_feature_table)
    else:
        df_train = pd.read_csv(args.training_feature_table, sep="\t")

    feature_cols = [c for c in df_train.columns if c not in ["sequence", "label"]]
    X_external = df_external[feature_cols]

    probs = model.predict_proba(X_external)[:, 1]
    preds = model.predict(X_external)

    df_external["predicted_label"] = preds
    df_external["prediction_probability"] = probs

    df_external = df_external.sort_values("prediction_probability", ascending=False)
    print(df_external[["sequence", "prediction_probability", "predicted_label"]])

    df_external.to_csv(args.output, sep="\t", index=False)
    print(f"Saved external validation predictions to: {args.output}")


if __name__ == "__main__":
    main()
    