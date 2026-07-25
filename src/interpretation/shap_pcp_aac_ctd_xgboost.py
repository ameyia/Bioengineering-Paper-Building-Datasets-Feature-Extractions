#!/usr/bin/env python3

import joblib
import pandas as pd
import numpy as np
import shap
import matplotlib.pyplot as plt
from pathlib import Path

# -----------------------------
# SETTINGS
# -----------------------------
DATA_PATH = "data/processed/baseline/pcp_aac_ctd.tsv"
MODEL_PATH = "results/models/pcp_aac_ctd/xgboost.pkl"
OUTPUT_DIR = Path("results/interpretability/shap_single_outputs/pcp_aac_ctd_xgboost")
BACKGROUND_SIZE = 100
EXPLAIN_SIZE = 100
TOP_N = 20
RANDOM_STATE = 42

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# -----------------------------
# LOAD DATA
# -----------------------------
df = pd.read_csv(DATA_PATH, sep="\t")
X = df.drop(columns=["sequence", "label"])
feature_names = list(X.columns)

X_background = X.sample(min(BACKGROUND_SIZE, len(X)), random_state=RANDOM_STATE)
X_explain = X.sample(min(EXPLAIN_SIZE, len(X)), random_state=RANDOM_STATE + 1)

# -----------------------------
# LOAD MODEL
# -----------------------------
model = joblib.load(MODEL_PATH)

# transform using pipeline
X_background_proc = model.named_steps["imputer"].transform(X_background)
X_explain_proc = model.named_steps["imputer"].transform(X_explain)

clf = model.named_steps["clf"]

X_background_df = pd.DataFrame(X_background_proc, columns=feature_names)
X_explain_df = pd.DataFrame(X_explain_proc, columns=feature_names)

# -----------------------------
# SHAP
# -----------------------------
explainer = shap.Explainer(clf, X_background_df)
shap_values = explainer(X_explain_df)

values = shap_values.values
if values.ndim == 3:
    values = values[:, :, 1]

explanation = shap.Explanation(
    values=values,
    base_values=np.array(shap_values.base_values),
    data=X_explain_df.values,
    feature_names=feature_names
)

# -----------------------------
# TOP FEATURES TABLE
# -----------------------------
mean_abs = np.abs(values).mean(axis=0)
top_df = pd.DataFrame({
    "feature": feature_names,
    "mean_abs_shap": mean_abs
}).sort_values("mean_abs_shap", ascending=False)

top_df.to_csv(OUTPUT_DIR / "top_features.tsv", sep="\t", index=False)

# -----------------------------
# PLOTS
# -----------------------------
plt.figure()
shap.plots.beeswarm(explanation, max_display=TOP_N, show=False)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "shap_beeswarm.png", dpi=300, bbox_inches="tight")
plt.close()

plt.figure()
shap.plots.bar(explanation, max_display=TOP_N, show=False)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "shap_bar.png", dpi=300, bbox_inches="tight")
plt.close()

print(f"Saved outputs to: {OUTPUT_DIR}")