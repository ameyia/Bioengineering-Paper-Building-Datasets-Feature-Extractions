from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
import pandas as pd
from propy.CTD import CalculateCTD

# -----------------------------
# DEFINE AMINO ACIDS + DIPEPTIDES
# -----------------------------
amino_acids = list("ACDEFGHIKLMNPQRSTVWY")
dipeptides = [a + b for a in amino_acids for b in amino_acids]

# -----------------------------
# 1. LOAD POSITIVES (FASTA)
# -----------------------------
positive_sequences = []

for record in SeqIO.parse("TableS1PositiveProteins.fasta", "fasta"):
    seq = str(record.seq).strip().upper()
    positive_sequences.append((seq, 1))

print(f"Loaded {len(positive_sequences)} positive sequences")

# -----------------------------
# 2. LOAD NEGATIVES (TSV)
# -----------------------------
neg_df = pd.read_csv("final_negatives.tsv", sep="\t", header=None)

negative_sequences = []

for i, row in neg_df.iterrows():
    seq = str(row[1]).strip().upper()
    negative_sequences.append((seq, 0))

print(f"Loaded {len(negative_sequences)} negative sequences")

# -----------------------------
# 3. COMBINE DATA
# -----------------------------
all_data = positive_sequences + negative_sequences
print(f"Total sequences: {len(all_data)}")

# -----------------------------
# 4. FEATURE EXTRACTION
# -----------------------------
features = []

for seq, label in all_data:
    seq = seq.strip().upper()
    
    # skip bad sequences
    if len(seq) == 0 or "X" in seq:
        continue
    
    try:
        analysis = ProteinAnalysis(seq)
        length = len(seq)
        
        # -------------------------
        # AAC
        # -------------------------
        aac = {}
        for aa in amino_acids:
            aac[aa] = (seq.count(aa) / length) * 100
        
        # -------------------------
        # BASE FEATURES (PCP + AAC)
        # -------------------------
        feature_dict = {
            "sequence": seq,
            "length": length,
            "mw": analysis.molecular_weight(),
            "aromaticity": analysis.aromaticity(),
            "pI": analysis.isoelectric_point(),
            "instability": analysis.instability_index(),
            "label": label
        }
        
        feature_dict.update(aac)

        # -------------------------
        # CTD (SAFE)
        # -------------------------
        try:
            ctd = CalculateCTD(seq)
            feature_dict.update(ctd)
        except Exception:
            print(f"CTD failed for: {seq}")

        # -------------------------
        # DPC (400 features)
        # -------------------------
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

        # store row
        features.append(feature_dict)
    
    except Exception as e:
        print(f"Error with sequence: {seq}")
        print(e)
        continue

# -----------------------------
# 5. CREATE DATAFRAMES
# -----------------------------
df_full = pd.DataFrame(features)

# PCP + AAC
base_columns = [
    "sequence", "length", "mw", "aromaticity",
    "pI", "instability", "label"
] + amino_acids

df_base = df_full[base_columns]

# PCP + AAC + CTD
ctd_columns = [col for col in df_full.columns if "_" in col]
df_ctd = df_full[base_columns + ctd_columns]

# FULL (PCP + AAC + CTD + DPC)
df_dpc = df_full

# -----------------------------
# 6. SAVE FILES
# -----------------------------
df_base.to_csv("pcp_aac.tsv", sep="\t", index=False)
df_ctd.to_csv("pcp_aac_ctd.tsv", sep="\t", index=False)
df_dpc.to_csv("pcp_aac_ctd_dpc.tsv", sep="\t", index=False)

print("\nSaved:")
print("- pcp_aac.tsv (25 features)")
print("- pcp_aac_ctd.tsv (172 features)")
print("- pcp_aac_ctd_dpc.tsv (~572 features)")

print("\nPreview:")
print(df_base.head())

# -----------------------------
# 7. COUNTING FEATURES 
# -----------------------------
# load each dataset
df_aac = pd.read_csv("pcp_aac.tsv", sep="\t")
df_ctd = pd.read_csv("pcp_aac_ctd.tsv", sep="\t")
df_dpc = pd.read_csv("pcp_aac_ctd_dpc.tsv", sep="\t")

# function to count features
def count_features(df):
    return len(df.columns) - 2  # subtract "sequence" and "label"

print("PCP + AAC:", count_features(df_aac))
print("PCP + AAC + CTD:", count_features(df_ctd))
print("PCP + AAC + CTD + DPC:", count_features(df_dpc))