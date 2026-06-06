#!/usr/bin/env python3

from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
import pandas as pd
from propy.CTD import CalculateCTD

amino_acids = list("ACDEFGHIKLMNPQRSTVWY")
dipeptides = [a + b for a in amino_acids for b in amino_acids]


def featurize_sequences(sequence_label_pairs):
    features = []

    for seq, label in sequence_label_pairs:
        seq = str(seq).strip().upper()

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
                print(f"CTD failed for: {seq}")

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
            features.append(feature_dict)

        except Exception as e:
            print(f"Error with sequence: {seq}")
            print(e)

    return pd.DataFrame(features)


# -----------------------------
# LOAD S1 TRAINING DATA
# -----------------------------
positive_sequences = []
for record in SeqIO.parse("TableS1PositiveProteins.fasta", "fasta"):
    seq = str(record.seq).strip().upper()
    positive_sequences.append((seq, 1))

print(f"Loaded {len(positive_sequences)} positive sequences")

neg_df = pd.read_csv("final_negatives.tsv", sep="\t", header=None)
negative_sequences = []
for _, row in neg_df.iterrows():
    seq = str(row[1]).strip().upper()
    negative_sequences.append((seq, 0))

print(f"Loaded {len(negative_sequences)} negative sequences")

all_data = positive_sequences + negative_sequences
print(f"Total sequences: {len(all_data)}")

df_full = featurize_sequences(all_data)

base_columns = [
    "sequence", "length", "mw", "aromaticity",
    "pI", "instability", "label"
] + amino_acids

ctd_columns = [col for col in df_full.columns if "_" in col]

df_base = df_full[base_columns]
df_ctd = df_full[base_columns + ctd_columns]
df_dpc = df_full

df_base.to_csv("pcp_aac.tsv", sep="\t", index=False)
df_ctd.to_csv("pcp_aac_ctd.tsv", sep="\t", index=False)
df_dpc.to_csv("pcp_aac_ctd_dpc.tsv", sep="\t", index=False)

print("\nSaved:")
print("- pcp_aac.tsv")
print("- pcp_aac_ctd.tsv")
print("- pcp_aac_ctd_dpc.tsv")

def count_features(df):
    return len(df.columns) - 2  # remove sequence and label

print("\nFeature counts:")
print("PCP + AAC:", count_features(df_base))
print("PCP + AAC + CTD:", count_features(df_ctd))
print("PCP + AAC + CTD + DPC:", count_features(df_dpc))