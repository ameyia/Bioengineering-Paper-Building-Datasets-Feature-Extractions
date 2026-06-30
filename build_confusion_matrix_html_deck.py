#!/usr/bin/env python3

from pathlib import Path

import pandas as pd


OUTDIR = Path("confusion_matrices")
DECK_PATH = OUTDIR / "confusion_matrix_slide_deck.html"

DATASET_LABELS = {
    "pcp_aac": "PCP + AAC",
    "pcp_aac_ctd": "PCP + AAC + CTD",
    "pcp_aac_ctd_dpc": "PCP + AAC + CTD + DPC",
}

MODEL_LABELS = {
    "svm_linear": "SVM linear",
    "svm_rbf": "SVM RBF",
    "logistic_regression": "Logistic regression",
    "random_forest": "Random forest",
    "knn": "kNN",
    "xgboost": "XGBoost",
}

DATASET_COLORS = {
    "pcp_aac": "#3f7fbd",
    "pcp_aac_ctd": "#4aa3a1",
    "pcp_aac_ctd_dpc": "#75b266",
}

MODEL_COLORS = {
    "svm_rbf": "#355c9c",
    "xgboost": "#6fbf73",
    "random_forest": "#9b6fc7",
    "logistic_regression": "#ef8a47",
    "svm_linear": "#58a6a6",
    "knn": "#d95f5f",
}


def fmt(value, digits=3):
    return f"{value:.{digits}f}"


def svg_bar_chart_selected(selected):
    width, height = 880, 360
    left, right, top, bottom = 70, 30, 45, 70
    plot_w = width - left - right
    plot_h = height - top - bottom
    max_y = 1.0
    rows = selected.sort_values("mean_test_mcc")
    bar_w = plot_w / len(rows) * 0.52
    gap = plot_w / len(rows)

    parts = [
        f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="Average test MCC by dataset">',
        '<rect width="100%" height="100%" fill="white" rx="8"/>',
    ]
    for tick in [0, 0.25, 0.5, 0.75, 1.0]:
        y = top + plot_h - tick / max_y * plot_h
        parts.append(f'<line x1="{left}" x2="{width-right}" y1="{y:.1f}" y2="{y:.1f}" stroke="#e3ebf3"/>')
        parts.append(f'<text x="{left-12}" y="{y+4:.1f}" text-anchor="end" font-size="12" fill="#657184">{tick:.2f}</text>')
    parts.append(f'<line x1="{left}" x2="{left}" y1="{top}" y2="{top+plot_h}" stroke="#2d4058"/>')
    parts.append(f'<line x1="{left}" x2="{width-right}" y1="{top+plot_h}" y2="{top+plot_h}" stroke="#2d4058"/>')

    for i, (_, row) in enumerate(rows.iterrows()):
        cx = left + gap * i + gap / 2
        value = row["mean_test_mcc"]
        sd = row["std_test_mcc"]
        bar_h = value / max_y * plot_h
        x = cx - bar_w / 2
        y = top + plot_h - bar_h
        color = DATASET_COLORS[row["dataset"]]
        err_top = top + plot_h - min(value + sd, max_y) / max_y * plot_h
        err_bottom = top + plot_h - max(value - sd, 0) / max_y * plot_h
        parts.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{bar_h:.1f}" fill="{color}"/>')
        parts.append(f'<line x1="{cx:.1f}" x2="{cx:.1f}" y1="{err_top:.1f}" y2="{err_bottom:.1f}" stroke="#1f3348" stroke-width="2"/>')
        parts.append(f'<line x1="{cx-12:.1f}" x2="{cx+12:.1f}" y1="{err_top:.1f}" y2="{err_top:.1f}" stroke="#1f3348" stroke-width="2"/>')
        parts.append(f'<line x1="{cx-12:.1f}" x2="{cx+12:.1f}" y1="{err_bottom:.1f}" y2="{err_bottom:.1f}" stroke="#1f3348" stroke-width="2"/>')
        parts.append(f'<text x="{cx:.1f}" y="{y-10:.1f}" text-anchor="middle" font-size="18" font-weight="700" fill="#17324d">{fmt(value)}</text>')
        parts.append(f'<text x="{cx:.1f}" y="{top+plot_h+34:.1f}" text-anchor="middle" font-size="14" fill="#17324d">{DATASET_LABELS[row["dataset"]]}</text>')

    parts.append('<text x="440" y="24" text-anchor="middle" font-size="18" font-weight="700" fill="#17324d">Average test MCC across 10 stratified splits</text>')
    parts.append('</svg>')
    return "\n".join(parts)


def svg_confusion_panels(selected):
    width, height = 1000, 300
    panel_w = 300
    parts = [
        f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="Average confusion matrices">',
        '<rect width="100%" height="100%" fill="white" rx="8"/>',
    ]
    ordered = selected.sort_values("mean_test_mcc")
    max_val = selected[["mean_tn", "mean_fp", "mean_fn", "mean_tp"]].max().max()
    for idx, (_, row) in enumerate(ordered.iterrows()):
        x0 = 35 + idx * 320
        y0 = 58
        cell = 86
        parts.append(f'<text x="{x0+cell}" y="28" text-anchor="middle" font-size="16" font-weight="700" fill="#17324d">{DATASET_LABELS[row["dataset"]]}</text>')
        cells = [
            (0, 0, row["mean_tn"], "TN"),
            (1, 0, row["mean_fp"], "FP"),
            (0, 1, row["mean_fn"], "FN"),
            (1, 1, row["mean_tp"], "TP"),
        ]
        for cx, cy, value, label in cells:
            intensity = int(245 - (value / max_val) * 120)
            fill = f"rgb({intensity},{intensity+8},{255})"
            x = x0 + cx * cell
            y = y0 + cy * cell
            parts.append(f'<rect x="{x}" y="{y}" width="{cell}" height="{cell}" fill="{fill}" stroke="#d5e0eb"/>')
            parts.append(f'<text x="{x+cell/2}" y="{y+34}" text-anchor="middle" font-size="14" fill="#657184">{label}</text>')
            parts.append(f'<text x="{x+cell/2}" y="{y+62}" text-anchor="middle" font-size="24" font-weight="700" fill="#17324d">{value:.1f}</text>')
        parts.append(f'<text x="{x0+cell}" y="{y0+cell*2+32}" text-anchor="middle" font-size="13" fill="#657184">Mean counts per test split</text>')
    parts.append('</svg>')
    return "\n".join(parts)


def svg_model_mcc(summary):
    width, height = 1000, 380
    left, top, bottom = 78, 50, 80
    plot_w = width - left - 30
    plot_h = height - top - bottom
    datasets = ["pcp_aac", "pcp_aac_ctd", "pcp_aac_ctd_dpc"]
    models = ["svm_rbf", "xgboost", "random_forest", "logistic_regression", "svm_linear", "knn"]
    x_gap = plot_w / len(models)
    bar_w = x_gap / 4.2
    indexed = summary.set_index(["dataset", "model"])

    parts = [
        f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="Model MCC by dataset">',
        '<rect width="100%" height="100%" fill="white" rx="8"/>',
    ]
    for tick in [0, 0.25, 0.5, 0.75, 1.0]:
        y = top + plot_h - tick * plot_h
        parts.append(f'<line x1="{left}" x2="{width-30}" y1="{y:.1f}" y2="{y:.1f}" stroke="#e3ebf3"/>')
        parts.append(f'<text x="{left-12}" y="{y+4:.1f}" text-anchor="end" font-size="12" fill="#657184">{tick:.2f}</text>')
    parts.append(f'<line x1="{left}" x2="{left}" y1="{top}" y2="{top+plot_h}" stroke="#2d4058"/>')
    parts.append(f'<line x1="{left}" x2="{width-30}" y1="{top+plot_h}" y2="{top+plot_h}" stroke="#2d4058"/>')
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
    legend_x = 250
    for idx, dataset in enumerate(datasets):
        x = legend_x + idx * 210
        parts.append(f'<rect x="{x}" y="18" width="16" height="16" fill="{DATASET_COLORS[dataset]}"/>')
        parts.append(f'<text x="{x+22}" y="31" font-size="13" fill="#17324d">{DATASET_LABELS[dataset]}</text>')
    parts.append('</svg>')
    return "\n".join(parts)


def svg_model_selection_counts(results):
    selected = results[results["selected_by_val_mcc"]]
    counts = selected.groupby(["dataset", "model"]).size().reset_index(name="count")
    datasets = ["pcp_aac", "pcp_aac_ctd", "pcp_aac_ctd_dpc"]
    models = ["xgboost", "svm_rbf", "random_forest", "logistic_regression", "svm_linear", "knn"]
    width, height = 900, 330
    left, top, plot_h = 85, 50, 205
    bar_w = 120
    gap = 250
    parts = [f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="Model selection counts">', '<rect width="100%" height="100%" fill="white" rx="8"/>']
    for tick in [0, 2, 4, 6, 8, 10]:
        y = top + plot_h - tick / 10 * plot_h
        parts.append(f'<line x1="{left}" x2="{width-30}" y1="{y:.1f}" y2="{y:.1f}" stroke="#e3ebf3"/>')
        parts.append(f'<text x="{left-12}" y="{y+4:.1f}" text-anchor="end" font-size="12" fill="#657184">{tick}</text>')
    for d_idx, dataset in enumerate(datasets):
        x = left + d_idx * gap + 60
        y_cursor = top + plot_h
        for model in models:
            row = counts[(counts["dataset"] == dataset) & (counts["model"] == model)]
            value = int(row["count"].iloc[0]) if not row.empty else 0
            if value == 0:
                continue
            h = value / 10 * plot_h
            y_cursor -= h
            parts.append(f'<rect x="{x}" y="{y_cursor:.1f}" width="{bar_w}" height="{h:.1f}" fill="{MODEL_COLORS[model]}"/>')
            if h > 18:
                parts.append(f'<text x="{x+bar_w/2}" y="{y_cursor+h/2+4:.1f}" text-anchor="middle" font-size="12" fill="white">{MODEL_LABELS[model]} {value}</text>')
        parts.append(f'<text x="{x+bar_w/2}" y="{top+plot_h+30}" text-anchor="middle" font-size="13" fill="#17324d">{DATASET_LABELS[dataset]}</text>')
    parts.append('</svg>')
    return "\n".join(parts)


def selected_table_rows(selected):
    rows = []
    for _, row in selected.iterrows():
        rows.append(
            f"<tr><td>{DATASET_LABELS[row['dataset']]}</td><td>{fmt(row['mean_test_mcc'])} +/- {fmt(row['std_test_mcc'])}</td><td>{fmt(row['mean_test_accuracy'])}</td><td>{fmt(row['mean_test_f1'])}</td><td>{row['mean_tn']:.1f}</td><td>{row['mean_fp']:.1f}</td><td>{row['mean_fn']:.1f}</td><td>{row['mean_tp']:.1f}</td></tr>"
        )
    return "\n".join(rows)


def top_model_rows(summary):
    rows = []
    for _, row in summary.sort_values("mean_test_mcc", ascending=False).head(8).iterrows():
        rows.append(
            f"<tr><td>{DATASET_LABELS[row['dataset']]}</td><td>{MODEL_LABELS[row['model']]}</td><td>{fmt(row['mean_test_mcc'])}</td><td>{fmt(row['mean_test_precision'])}</td><td>{fmt(row['mean_test_recall'])}</td><td>{row['mean_fp']:.1f}</td><td>{row['mean_fn']:.1f}</td></tr>"
        )
    return "\n".join(rows)


def representative_plot_paths(results, selected_summary):
    selected = results[results["selected_by_val_mcc"]].copy()
    paths = []
    for _, summary_row in selected_summary.iterrows():
        dataset = summary_row["dataset"]
        subset = selected[selected["dataset"] == dataset].copy()
        subset["distance"] = (subset["test_mcc"] - summary_row["mean_test_mcc"]).abs()
        row = subset.sort_values("distance").iloc[0]
        path = f"plots/{dataset}_seed{int(row['seed'])}_{row['model']}.png"
        paths.append((DATASET_LABELS[dataset], row["model"], int(row["seed"]), row["test_mcc"], path))
    return paths


def build_html(results, summary, selected):
    selected_order = ["pcp_aac_ctd_dpc", "pcp_aac_ctd", "pcp_aac"]
    selected = selected.set_index("dataset").loc[selected_order].reset_index()
    best = selected.iloc[0]
    reps = representative_plot_paths(results, selected)

    rep_cards = "\n".join(
        f"""
        <div class="card">
          <h3>{label}</h3>
          <img class="matrix" src="{path}" alt="{label} representative confusion matrix">
          <p class="small">Seed {seed}, {MODEL_LABELS[model]}, MCC {mcc:.3f}</p>
        </div>
        """
        for label, model, seed, mcc, path in reps
    )

    html = f"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>Confusion Matrix Results</title>
<style>
body {{ margin: 0; background: #dfe7f1; font-family: Arial, Helvetica, sans-serif; color: #17324d; }}
.slide {{ width: 1280px; height: 720px; box-sizing: border-box; margin: 28px auto; padding: 54px 70px 42px; background: #f6f8fb; border: 1px solid #c9d4e1; position: relative; overflow: hidden; box-shadow: 0 8px 24px rgba(20,40,70,.14); page-break-after: always; }}
.topbar {{ position: absolute; top: 0; left: 0; right: 0; height: 16px; background: #6551d8; }}
h1 {{ font-size: 42px; margin: 0 0 8px; letter-spacing: 0; }}
h2 {{ font-size: 32px; margin: 0 0 8px; letter-spacing: 0; }}
h3 {{ margin: 0 0 10px; font-size: 20px; }}
.subtitle {{ color: #657184; font-size: 18px; margin-bottom: 22px; }}
.grid2 {{ display: grid; grid-template-columns: 1fr 1fr; gap: 30px; align-items: center; }}
.grid3 {{ display: grid; grid-template-columns: repeat(3, 1fr); gap: 20px; }}
.card {{ background: white; border: 1.5px solid #dbe4ee; border-radius: 8px; padding: 18px 20px; }}
.metric {{ font-size: 44px; font-weight: 700; margin-bottom: 6px; }}
.metric-label {{ color: #657184; font-size: 17px; line-height: 1.25; }}
.callout {{ background: #eef7f2; border-left: 6px solid #75b266; padding: 17px 22px; font-size: 22px; line-height: 1.25; }}
table {{ width: 100%; border-collapse: collapse; font-size: 17px; background: white; border: 1px solid #dbe4ee; }}
th {{ text-align: left; background: #2f7dbd; color: white; padding: 10px; font-weight: 700; }}
td {{ padding: 9px 10px; border-top: 1px solid #e6edf5; }}
tr:nth-child(even) td {{ background: #f1f5fa; }}
ul {{ font-size: 22px; line-height: 1.35; padding-left: 24px; }}
li {{ margin-bottom: 10px; }}
svg {{ width: 100%; max-height: 500px; }}
img.matrix {{ width: 100%; height: 300px; object-fit: contain; display: block; }}
.small {{ font-size: 14px; color: #657184; }}
.footer {{ position: absolute; left: 70px; right: 70px; bottom: 18px; color: #6d7989; font-size: 12px; border-top: 1px solid #cbd7e5; padding-top: 8px; }}
@media print {{ body {{ background: white; }} .slide {{ margin: 0; box-shadow: none; border: none; page-break-after: always; }} }}
</style>
</head>
<body>

<section class="slide">
<div class="topbar"></div>
<h1>Confusion Matrix Results</h1>
<div class="subtitle">Held-out test performance across 3 feature datasets, 10 stratified splits, and 6 model families</div>
<div class="grid3" style="margin-top: 48px;">
<div class="card"><div class="metric">{fmt(best['mean_test_mcc'])}</div><div class="metric-label">best average test MCC<br>{DATASET_LABELS[best['dataset']]}</div></div>
<div class="card"><div class="metric">{fmt(best['mean_test_accuracy'])}</div><div class="metric-label">average test accuracy<br>validation-selected models</div></div>
<div class="card"><div class="metric">{fmt(best['mean_test_f1'])}</div><div class="metric-label">average test F1<br>validation-selected models</div></div>
</div>
<div class="callout" style="margin-top: 58px;">The confusion matrices support the same story as the performance analysis: adding CTD improves prediction, adding DPC improves it further, and SVM RBF is the strongest model family for the richer feature sets.</div>
<div class="footer">Generated from confusion_matrices/*.csv</div>
</section>

<section class="slide">
<div class="topbar"></div>
<h2>How To Read The Matrices</h2>
<div class="subtitle">Each held-out test split has about 88 negatives and 17 positives</div>
<div class="grid2">
<div>
<table><tr><th></th><th>Predicted negative</th><th>Predicted antithrombotic</th></tr><tr><td><b>True negative</b></td><td>TN</td><td>FP</td></tr><tr><td><b>True antithrombotic</b></td><td>FN</td><td>TP</td></tr></table>
<ul style="margin-top: 28px;"><li>False positives are negatives incorrectly flagged as antithrombotic.</li><li>False negatives are true antithrombotic peptides missed by the model.</li><li>MCC is prioritized because the classes are moderately imbalanced.</li></ul>
</div>
<div class="card"><h2 style="font-size: 25px;">Evaluation design</h2><ul><li>60 percent training, 20 percent validation, 20 percent test.</li><li>Splits are stratified, preserving the positive/negative ratio.</li><li>Best model per seed is selected by validation MCC.</li><li>The held-out test split is used for the final confusion matrix.</li></ul></div>
</div>
<div class="footer">Matrix order: [[TN, FP], [FN, TP]]</div>
</section>

<section class="slide">
<div class="topbar"></div>
<h2>Validation-Selected Results</h2>
<div class="subtitle">Average held-out test MCC across 10 seeds</div>
{svg_bar_chart_selected(selected)}
<div class="footer">Error bars show standard deviation across random seeds.</div>
</section>

<section class="slide">
<div class="topbar"></div>
<h2>Average Confusion Matrices</h2>
<div class="subtitle">Mean test-set counts for validation-selected models</div>
{svg_confusion_panels(selected)}
<div class="footer">DPC has the fewest false positives and highest average true positives among validation-selected models.</div>
</section>

<section class="slide">
<div class="topbar"></div>
<h2>Selected-Model Summary Table</h2>
<div class="subtitle">These are the models chosen by validation MCC in each seed</div>
<table><tr><th>Feature dataset</th><th>Test MCC</th><th>Accuracy</th><th>F1</th><th>TN</th><th>FP</th><th>FN</th><th>TP</th></tr>{selected_table_rows(selected)}</table>
<div class="callout" style="margin-top: 28px;">PCP + AAC + CTD + DPC is the strongest dataset: higher MCC, higher F1, fewer false positives, and more true positives.</div>
<div class="footer">Counts are averages per held-out test split.</div>
</section>

<section class="slide">
<div class="topbar"></div>
<h2>Model Comparison</h2>
<div class="subtitle">Average test MCC for each model family and feature dataset</div>
{svg_model_mcc(summary)}
<div class="footer">SVM RBF has the best average MCC in all three feature spaces in this confusion-matrix audit.</div>
</section>

<section class="slide">
<div class="topbar"></div>
<h2>Model Selection Stability</h2>
<div class="subtitle">How often each model was selected by validation MCC across 10 seeds</div>
{svg_model_selection_counts(results)}
<div class="footer">SVM RBF is most often selected for the CTD and CTD + DPC feature sets.</div>
</section>

<section class="slide">
<div class="topbar"></div>
<h2>Representative Matrix Plots</h2>
<div class="subtitle">One validation-selected test matrix per feature dataset, chosen near the dataset average MCC</div>
<div class="grid3">{rep_cards}</div>
<div class="footer">Full set of matrix plots is in confusion_matrices/plots/.</div>
</section>

<section class="slide">
<div class="topbar"></div>
<h2>Top Model Rows</h2>
<div class="subtitle">Highest average test MCC combinations</div>
<table><tr><th>Feature dataset</th><th>Model</th><th>MCC</th><th>Precision</th><th>Recall</th><th>FP</th><th>FN</th></tr>{top_model_rows(summary)}</table>
<div class="footer">Precision reflects false-positive control; recall reflects how many true positives were recovered.</div>
</section>

<section class="slide">
<div class="topbar"></div>
<h2>Conclusion</h2>
<div class="subtitle">What the confusion matrices add to the story</div>
<div class="grid2">
<div class="card"><h2 style="font-size: 25px;">Main result</h2><ul><li>DPC gives the strongest overall test performance.</li><li>CTD improves over PCP + AAC alone.</li><li>SVM RBF is the most reliable model family for richer features.</li></ul></div>
<div class="card"><h2 style="font-size: 25px;">Recommended wording</h2><ul><li>Use "validation-selected held-out test results."</li><li>Report MCC plus false positives and false negatives.</li><li>Use one representative matrix on the main slide and keep the full audit in backup.</li></ul></div>
</div>
<div class="callout" style="margin-top: 32px;">Best concise statement: PCP + AAC + CTD + DPC with SVM RBF gave the strongest balance of sensitivity and specificity, with very low false-positive counts and the highest average MCC.</div>
<div class="footer">Full results: confusion_matrix_results_by_model_seed.csv</div>
</section>

</body></html>"""
    DECK_PATH.write_text(html)


def main():
    results = pd.read_csv(OUTDIR / "confusion_matrix_results_by_model_seed.csv")
    summary = pd.read_csv(OUTDIR / "confusion_matrix_summary_by_dataset_model.csv")
    selected = pd.read_csv(OUTDIR / "confusion_matrix_summary_validation_selected.csv")
    build_html(results, summary, selected)
    print(f"Saved slide deck: {DECK_PATH}")


if __name__ == "__main__":
    main()
