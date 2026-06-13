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
DATA_PATH = "pcp_aac.tsv"
MODEL_PATH = "saved_models/pcp_aac/knn.pkl"
OUTPUT_DIR = Path("shap_single_outputs/pcp_aac_knn")
BACKGROUND_SIZE = 60
EXPLAIN_SIZE = 40
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

X_background_proc = model.named_steps["imputer"].transform(X_background)
X_background_proc = model.named_steps["scaler"].transform(X_background_proc)

X_explain_proc = model.named_steps["imputer"].transform(X_explain)
X_explain_proc = model.named_steps["scaler"].transform(X_explain_proc)

clf = model.named_steps["clf"]

X_background_df = pd.DataFrame(X_background_proc, columns=feature_names)
X_explain_df = pd.DataFrame(X_explain_proc, columns=feature_names)

# -----------------------------
# SHAP
# -----------------------------
predict_fn = lambda X_: clf.predict_proba(X_)[:, 1]
max_evals = 2 * X_explain_df.shape[1] + 1

explainer = shap.Explainer(predict_fn, X_background_df, max_evals=max_evals)
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