#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import matthews_corrcoef
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


OUTDIR = Path("pca_plots")
DECK_PATH = OUTDIR / "pca_dataset_model_comparison.html"

TABLES = ["pcp_aac.tsv", "pcp_aac_ctd.tsv", "pcp_aac_ctd_dpc.tsv"]
DATASET_LABELS = {
    "pcp_aac": "PCP + AAC",
    "pcp_aac_ctd": "PCP + AAC + CTD",
    "pcp_aac_ctd_dpc": "PCP + AAC + CTD + DPC",
}
DATASET_COLORS = {
    "pcp_aac": "#3f7fbd",
    "pcp_aac_ctd": "#4aa3a1",
    "pcp_aac_ctd_dpc": "#75b266",
}
MODEL_LABELS = {
    "svm_linear": "SVM linear",
    "svm_rbf": "SVM RBF",
    "logistic_regression": "Logistic regression",
    "random_forest": "Random forest",
    "knn": "kNN",
    "xgboost": "XGBoost",
}


def load_table(path):
    df = pd.read_csv(path, sep="\t")
    y = df["label"].astype(int)
    sequences = df["sequence"].astype(str) if "sequence" in df.columns else pd.Series(range(len(df)))
    drop_cols = ["label"]
    if "sequence" in df.columns:
        drop_cols.append("sequence")
    X = df.drop(columns=drop_cols)
    X = X.apply(pd.to_numeric, errors="coerce")
    X = X.replace([np.inf, -np.inf], np.nan)
    X = X.fillna(X.median(numeric_only=True))
    return X, y, sequences


def get_models(random_state):
    models = {
        "svm_linear": Pipeline(
            [
                ("scaler", StandardScaler()),
                (
                    "model",
                    SVC(
                        kernel="linear",
                        probability=True,
                        class_weight="balanced",
                        random_state=random_state,
                    ),
                ),
            ]
        ),
        "svm_rbf": Pipeline(
            [
                ("scaler", StandardScaler()),
                (
                    "model",
                    SVC(
                        kernel="rbf",
                        C=10,
                        gamma="scale",
                        probability=True,
                        class_weight="balanced",
                        random_state=random_state,
                    ),
                ),
            ]
        ),
        "logistic_regression": Pipeline(
            [
                ("scaler", StandardScaler()),
                (
                    "model",
                    LogisticRegression(
                        max_iter=5000,
                        class_weight="balanced",
                        random_state=random_state,
                    ),
                ),
            ]
        ),
        "random_forest": RandomForestClassifier(
            n_estimators=500,
            min_samples_leaf=2,
            class_weight="balanced",
            random_state=random_state,
            n_jobs=-1,
        ),
        "knn": Pipeline(
            [
                ("scaler", StandardScaler()),
                ("model", KNeighborsClassifier(n_neighbors=5)),
            ]
        ),
    }

    if HAS_XGB:
        models["xgboost"] = XGBClassifier(
            n_estimators=300,
            max_depth=3,
            learning_rate=0.03,
            subsample=0.8,
            colsample_bytree=0.8,
            eval_metric="logloss",
            random_state=random_state,
            n_jobs=-1,
        )

    return models


def split_data(X, y, seed):
    X_train, X_temp, y_train, y_temp = train_test_split(
        X,
        y,
        test_size=0.40,
        stratify=y,
        random_state=seed,
    )
    X_val, X_test, y_val, y_test = train_test_split(
        X_temp,
        y_temp,
        test_size=0.50,
        stratify=y_temp,
        random_state=seed,
    )
    return X_train, X_val, X_test, y_train, y_val, y_test


def pca_scores(X):
    scaled = StandardScaler().fit_transform(X)
    pca = PCA(n_components=2, random_state=0)
    scores = pca.fit_transform(scaled)
    return scores, pca.explained_variance_ratio_


def scale_points(points, width, height, pad):
    x = points[:, 0]
    y = points[:, 1]
    xmin, xmax = x.min(), x.max()
    ymin, ymax = y.min(), y.max()
    if xmin == xmax:
        xmin -= 1
        xmax += 1
    if ymin == ymax:
        ymin -= 1
        ymax += 1
    sx = pad + (x - xmin) / (xmax - xmin) * (width - 2 * pad)
    sy = height - pad - (y - ymin) / (ymax - ymin) * (height - 2 * pad)
    return sx, sy


def svg_pca_by_label(dataset, scores, y, variance, width=360, height=310):
    pad = 42
    sx, sy = scale_points(scores, width, height, pad)
    parts = [
        f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="{DATASET_LABELS[dataset]} PCA">',
        '<rect width="100%" height="100%" rx="8" fill="white"/>',
        f'<text x="{width/2}" y="24" text-anchor="middle" font-size="16" font-weight="700" fill="#17324d">{DATASET_LABELS[dataset]}</text>',
        f'<text x="{width/2}" y="{height-8}" text-anchor="middle" font-size="11" fill="#657184">PC1 {variance[0]*100:.1f}% | PC2 {variance[1]*100:.1f}%</text>',
        f'<line x1="{pad}" x2="{width-pad}" y1="{height-pad}" y2="{height-pad}" stroke="#cbd7e5"/>',
        f'<line x1="{pad}" x2="{pad}" y1="{pad}" y2="{height-pad}" stroke="#cbd7e5"/>',
    ]
    for label_value in [0, 1]:
        color = "#4a6f9f" if label_value == 0 else "#d95f5f"
        stroke = "#29435f" if label_value == 0 else "#8d3030"
        radius = 3.2 if label_value == 0 else 4.2
        idx = np.where(y.to_numpy() == label_value)[0]
        for i in idx:
            parts.append(
                f'<circle cx="{sx[i]:.1f}" cy="{sy[i]:.1f}" r="{radius}" fill="{color}" stroke="{stroke}" stroke-width="0.5" opacity="0.76"/>'
            )
    parts.append(f'<circle cx="{width-100}" cy="48" r="4" fill="#4a6f9f"/><text x="{width-90}" y="52" font-size="12" fill="#17324d">negative</text>')
    parts.append(f'<circle cx="{width-100}" cy="68" r="5" fill="#d95f5f"/><text x="{width-90}" y="72" font-size="12" fill="#17324d">positive</text>')
    parts.append("</svg>")
    return "\n".join(parts)


def validation_selected_predictions(X, y, seed=0):
    X_train, X_val, X_test, y_train, y_val, y_test = split_data(X, y, seed)
    models = get_models(seed)
    val_scores = {}
    fitted = {}
    for model_name, model in models.items():
        model.fit(X_train, y_train)
        val_pred = model.predict(X_val)
        val_scores[model_name] = matthews_corrcoef(y_val, val_pred)
        fitted[model_name] = model
    best_model_name = max(val_scores, key=val_scores.get)
    best_model = fitted[best_model_name]
    test_pred = best_model.predict(X_test)
    return best_model_name, X_test.index.to_numpy(), y_test.to_numpy(), test_pred


def svg_pca_test_outcomes(dataset, scores, test_idx, y_test, pred, best_model, width=360, height=310):
    pad = 42
    sx, sy = scale_points(scores, width, height, pad)
    colors = {
        "TN": "#4a6f9f",
        "FP": "#f2a43a",
        "FN": "#d94f4f",
        "TP": "#4ea96b",
    }
    parts = [
        f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="{DATASET_LABELS[dataset]} PCA model outcomes">',
        '<rect width="100%" height="100%" rx="8" fill="white"/>',
        f'<text x="{width/2}" y="24" text-anchor="middle" font-size="15" font-weight="700" fill="#17324d">{DATASET_LABELS[dataset]}</text>',
        f'<text x="{width/2}" y="43" text-anchor="middle" font-size="12" fill="#657184">Seed 0 validation-selected: {MODEL_LABELS[best_model]}</text>',
        f'<line x1="{pad}" x2="{width-pad}" y1="{height-pad}" y2="{height-pad}" stroke="#cbd7e5"/>',
        f'<line x1="{pad}" x2="{pad}" y1="{pad}" y2="{height-pad}" stroke="#cbd7e5"/>',
    ]
    # Show all samples faintly so the held-out test points have context.
    for i in range(scores.shape[0]):
        parts.append(f'<circle cx="{sx[i]:.1f}" cy="{sy[i]:.1f}" r="2.4" fill="#cbd7e5" opacity="0.35"/>')
    for idx, truth, prediction in zip(test_idx, y_test, pred):
        if truth == 0 and prediction == 0:
            key = "TN"
        elif truth == 0 and prediction == 1:
            key = "FP"
        elif truth == 1 and prediction == 0:
            key = "FN"
        else:
            key = "TP"
        parts.append(
            f'<circle cx="{sx[idx]:.1f}" cy="{sy[idx]:.1f}" r="5.2" fill="{colors[key]}" stroke="#24384e" stroke-width="0.8" opacity="0.9"><title>{key}</title></circle>'
        )
    legend = [("TN", 62), ("FP", 82), ("FN", 102), ("TP", 122)]
    for key, y0 in legend:
        parts.append(f'<circle cx="{width-92}" cy="{y0}" r="4.5" fill="{colors[key]}"/><text x="{width-82}" y="{y0+4}" font-size="12" fill="#17324d">{key}</text>')
    parts.append("</svg>")
    return "\n".join(parts)


def svg_model_mcc(summary):
    width, height = 980, 350
    left, top, bottom = 76, 48, 78
    plot_w = width - left - 30
    plot_h = height - top - bottom
    datasets = ["pcp_aac", "pcp_aac_ctd", "pcp_aac_ctd_dpc"]
    models = ["svm_rbf", "xgboost", "random_forest", "logistic_regression", "svm_linear", "knn"]
    x_gap = plot_w / len(models)
    bar_w = x_gap / 4.2
    indexed = summary.set_index(["dataset", "model"])
    parts = [f'<svg viewBox="0 0 {width} {height}">', '<rect width="100%" height="100%" fill="white" rx="8"/>']
    for tick in [0, 0.25, 0.5, 0.75, 1.0]:
        y = top + plot_h - tick * plot_h
        parts.append(f'<line x1="{left}" x2="{width-30}" y1="{y:.1f}" y2="{y:.1f}" stroke="#e3ebf3"/>')
        parts.append(f'<text x="{left-12}" y="{y+4:.1f}" text-anchor="end" font-size="12" fill="#657184">{tick:.2f}</text>')
    for m_idx, model in enumerate(models):
        center = left + x_gap * m_idx + x_gap / 2
        for d_idx, dataset in enumerate(datasets):
            if (dataset, model) not in indexed.index:
                continue
            value = indexed.loc[(dataset, model), "mean_test_mcc"]
            bar_h = value * plot_h
            x = center + (d_idx - 1) * bar_w - bar_w / 2
            y = top + plot_h - bar_h
            parts.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{bar_h:.1f}" fill="{DATASET_COLORS[dataset]}"/>')
        parts.append(f'<text x="{center:.1f}" y="{top+plot_h+28}" text-anchor="middle" font-size="12" fill="#17324d">{MODEL_LABELS[model]}</text>')
    for idx, dataset in enumerate(datasets):
        x = 220 + idx * 210
        parts.append(f'<rect x="{x}" y="16" width="16" height="16" fill="{DATASET_COLORS[dataset]}"/>')
        parts.append(f'<text x="{x+22}" y="29" font-size="13" fill="#17324d">{DATASET_LABELS[dataset]}</text>')
    parts.append("</svg>")
    return "\n".join(parts)


def build_html(dataset_svgs, outcome_svgs, model_svg):
    html = f"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>PCA Dataset and Model Comparison</title>
<style>
body {{ margin: 0; background: #dfe7f1; font-family: Arial, Helvetica, sans-serif; color: #17324d; }}
.slide {{ width: 1280px; min-height: 720px; box-sizing: border-box; margin: 28px auto; padding: 54px 70px 42px; background: #f6f8fb; border: 1px solid #c9d4e1; position: relative; overflow: hidden; box-shadow: 0 8px 24px rgba(20,40,70,.14); page-break-after: always; }}
.topbar {{ position: absolute; top: 0; left: 0; right: 0; height: 16px; background: #6551d8; }}
h1 {{ font-size: 42px; margin: 0 0 8px; }}
h2 {{ font-size: 32px; margin: 0 0 8px; }}
.subtitle {{ color: #657184; font-size: 18px; margin-bottom: 22px; }}
.grid3 {{ display: grid; grid-template-columns: repeat(3, 1fr); gap: 18px; }}
.card {{ background: white; border: 1.5px solid #dbe4ee; border-radius: 8px; padding: 18px 20px; }}
.callout {{ background: #eef7f2; border-left: 6px solid #75b266; padding: 17px 22px; font-size: 21px; line-height: 1.28; margin-top: 22px; }}
svg {{ width: 100%; }}
ul {{ font-size: 22px; line-height: 1.35; padding-left: 24px; }}
li {{ margin-bottom: 10px; }}
.footer {{ position: absolute; left: 70px; right: 70px; bottom: 18px; color: #6d7989; font-size: 12px; border-top: 1px solid #cbd7e5; padding-top: 8px; }}
@media print {{ body {{ background: white; }} .slide {{ margin: 0; box-shadow: none; border: none; page-break-after: always; }} }}
</style>
</head>
<body>

<section class="slide">
<div class="topbar"></div>
<h1>PCA Comparison</h1>
<div class="subtitle">Feature-space visualization for PCP + AAC, PCP + AAC + CTD, and PCP + AAC + CTD + DPC</div>
<div class="grid3">{''.join(f'<div class="card">{svg}</div>' for svg in dataset_svgs)}</div>
<div class="callout">PCA is unsupervised, so it does not prove separability by itself. It shows whether positives and negatives occupy visibly different regions after scaling each feature set.</div>
<div class="footer">Points are all peptides in each feature table; red = positive, blue = negative.</div>
</section>

<section class="slide">
<div class="topbar"></div>
<h2>Model Outcomes on PCA Space</h2>
<div class="subtitle">Seed 0 held-out test predictions from the validation-selected model for each dataset</div>
<div class="grid3">{''.join(f'<div class="card">{svg}</div>' for svg in outcome_svgs)}</div>
<div class="footer">Gray points are all peptides for context; colored points are held-out test predictions.</div>
</section>

<section class="slide">
<div class="topbar"></div>
<h2>Model Comparison</h2>
<div class="subtitle">Average test MCC from the confusion-matrix audit</div>
<div class="card">{model_svg}</div>
<div class="callout">Use the PCA plots as visual context and the MCC/confusion-matrix results as the quantitative evidence. The richest DPC feature set and SVM RBF remain the strongest overall combination.</div>
<div class="footer">MCC values are averaged across seeds.</div>
</section>

</body>
</html>"""
    DECK_PATH.write_text(html)


def main():
    OUTDIR.mkdir(exist_ok=True)
    dataset_svgs = []
    outcome_svgs = []
    for table in TABLES:
        dataset = Path(table).stem
        X, y, _ = load_table(table)
        scores, variance = pca_scores(X)
        dataset_svgs.append(svg_pca_by_label(dataset, scores, y, variance))
        best_model, test_idx, y_test, pred = validation_selected_predictions(X, y, seed=0)
        outcome_svgs.append(svg_pca_test_outcomes(dataset, scores, test_idx, y_test, pred, best_model))

    summary = pd.read_csv("confusion_matrices/confusion_matrix_summary_by_dataset_model.csv")
    model_svg = svg_model_mcc(summary)
    build_html(dataset_svgs, outcome_svgs, model_svg)
    print(f"Saved PCA comparison deck: {DECK_PATH}")


if __name__ == "__main__":
    main()
