#!/usr/bin/env python3

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


OUTDIR = Path("confusion_matrices")
ASSET_DIR = OUTDIR / "slide_assets"
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

COLORS = {
    "pcp_aac": "#3f7fbd",
    "pcp_aac_ctd": "#4aa3a1",
    "pcp_aac_ctd_dpc": "#75b266",
}


def savefig(path):
    plt.tight_layout()
    plt.savefig(path, dpi=220, bbox_inches="tight")
    plt.close()


def style_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#e6edf5", linewidth=0.8)
    ax.set_axisbelow(True)


def make_selected_performance_chart(selected):
    labels = [DATASET_LABELS[d] for d in selected["dataset"]]
    x = np.arange(len(labels))

    fig, ax = plt.subplots(figsize=(9.5, 5.2))
    bars = ax.bar(
        x,
        selected["mean_test_mcc"],
        yerr=selected["std_test_mcc"],
        capsize=5,
        color=[COLORS[d] for d in selected["dataset"]],
        edgecolor="#22364f",
        linewidth=0.8,
    )
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("Average test MCC")
    ax.set_title("Validation-selected model performance across 10 stratified splits")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=15, ha="right")
    style_axes(ax)

    for bar, value in zip(bars, selected["mean_test_mcc"]):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value + 0.035,
            f"{value:.3f}",
            ha="center",
            va="bottom",
            fontsize=11,
            fontweight="bold",
            color="#1b334f",
        )

    out = ASSET_DIR / "selected_mcc_by_dataset.png"
    savefig(out)
    return out


def make_selected_confusion_chart(selected):
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.7))

    max_value = selected[["mean_tn", "mean_fp", "mean_fn", "mean_tp"]].to_numpy().max()
    for ax, (_, row) in zip(axes, selected.iterrows()):
        cm = np.array(
            [
                [row["mean_tn"], row["mean_fp"]],
                [row["mean_fn"], row["mean_tp"]],
            ]
        )
        im = ax.imshow(cm, cmap="Blues", vmin=0, vmax=max_value)
        ax.set_title(DATASET_LABELS[row["dataset"]], fontsize=11, fontweight="bold")
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Pred neg", "Pred pos"], fontsize=9)
        ax.set_yticks([0, 1])
        ax.set_yticklabels(["True neg", "True pos"], fontsize=9)

        for i in range(2):
            for j in range(2):
                value = cm[i, j]
                color = "white" if value > max_value * 0.5 else "#14304a"
                ax.text(j, i, f"{value:.1f}", ha="center", va="center", color=color, fontsize=13, fontweight="bold")

    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.72, label="Mean count")
    fig.suptitle("Average held-out test confusion matrix for validation-selected models", y=1.05)
    out = ASSET_DIR / "selected_average_confusion_matrices.png"
    savefig(out)
    return out


def make_model_mcc_chart(summary):
    order = ["pcp_aac", "pcp_aac_ctd", "pcp_aac_ctd_dpc"]
    models = ["svm_rbf", "xgboost", "random_forest", "logistic_regression", "svm_linear", "knn"]
    available = [m for m in models if m in set(summary["model"])]
    x = np.arange(len(available))
    width = 0.25

    fig, ax = plt.subplots(figsize=(11, 5.2))
    for idx, dataset in enumerate(order):
        subset = summary[summary["dataset"] == dataset].set_index("model")
        values = [subset.loc[m, "mean_test_mcc"] if m in subset.index else np.nan for m in available]
        ax.bar(
            x + (idx - 1) * width,
            values,
            width=width,
            label=DATASET_LABELS[dataset],
            color=COLORS[dataset],
            alpha=0.92,
        )

    ax.set_ylim(0, 1.0)
    ax.set_ylabel("Average test MCC")
    ax.set_title("Average test MCC by model and feature dataset")
    ax.set_xticks(x)
    ax.set_xticklabels([MODEL_LABELS[m] for m in available], rotation=20, ha="right")
    ax.legend(frameon=False, ncols=3, loc="upper center", bbox_to_anchor=(0.5, 1.16))
    style_axes(ax)

    out = ASSET_DIR / "model_mcc_by_dataset.png"
    savefig(out)
    return out


def make_error_tradeoff_chart(summary):
    fig, ax = plt.subplots(figsize=(8.8, 5.2))
    for dataset, subset in summary.groupby("dataset"):
        ax.scatter(
            subset["mean_fp"],
            subset["mean_fn"],
            s=110,
            color=COLORS[dataset],
            edgecolor="#203449",
            label=DATASET_LABELS[dataset],
            alpha=0.92,
        )
        for _, row in subset.iterrows():
            ax.text(row["mean_fp"] + 0.15, row["mean_fn"] + 0.05, MODEL_LABELS[row["model"]], fontsize=8)

    ax.set_xlabel("Average false positives per test split")
    ax.set_ylabel("Average false negatives per test split")
    ax.set_title("Error tradeoff: false positives vs false negatives")
    ax.legend(frameon=False, loc="upper right")
    style_axes(ax)

    out = ASSET_DIR / "false_positive_false_negative_tradeoff.png"
    savefig(out)
    return out


def make_selected_model_counts_chart(results):
    selected = results[results["selected_by_val_mcc"]].copy()
    counts = selected.groupby(["dataset", "model"]).size().reset_index(name="count")
    order = ["pcp_aac", "pcp_aac_ctd", "pcp_aac_ctd_dpc"]
    models = ["svm_rbf", "xgboost", "random_forest", "logistic_regression", "svm_linear", "knn"]
    available = [m for m in models if m in set(counts["model"])]

    fig, ax = plt.subplots(figsize=(10.4, 5.1))
    bottom = np.zeros(len(order))
    model_colors = {
        "svm_rbf": "#355c9c",
        "xgboost": "#6fbf73",
        "random_forest": "#9b6fc7",
        "logistic_regression": "#ef8a47",
        "svm_linear": "#58a6a6",
        "knn": "#d95f5f",
    }

    for model in available:
        values = []
        for dataset in order:
            row = counts[(counts["dataset"] == dataset) & (counts["model"] == model)]
            values.append(int(row["count"].iloc[0]) if not row.empty else 0)
        ax.bar(
            [DATASET_LABELS[d] for d in order],
            values,
            bottom=bottom,
            label=MODEL_LABELS[model],
            color=model_colors[model],
        )
        bottom += np.array(values)

    ax.set_ylim(0, 10)
    ax.set_ylabel("Times selected across 10 seeds")
    ax.set_title("Which model was selected by validation MCC?")
    ax.legend(frameon=False, ncols=3, loc="upper center", bbox_to_anchor=(0.5, 1.18))
    style_axes(ax)

    out = ASSET_DIR / "selected_model_counts.png"
    savefig(out)
    return out


def format_metric(value, digits=3):
    return f"{value:.{digits}f}"


def table_rows_selected(selected):
    rows = []
    for _, row in selected.iterrows():
        rows.append(
            f"""
            <tr>
              <td>{DATASET_LABELS[row['dataset']]}</td>
              <td>{format_metric(row['mean_test_mcc'])} +/- {format_metric(row['std_test_mcc'])}</td>
              <td>{format_metric(row['mean_test_accuracy'])}</td>
              <td>{format_metric(row['mean_test_f1'])}</td>
              <td>{row['mean_tn']:.1f}</td>
              <td>{row['mean_fp']:.1f}</td>
              <td>{row['mean_fn']:.1f}</td>
              <td>{row['mean_tp']:.1f}</td>
            </tr>
            """
        )
    return "\n".join(rows)


def table_rows_top_models(summary):
    rows = []
    top = summary.sort_values("mean_test_mcc", ascending=False).head(8)
    for _, row in top.iterrows():
        rows.append(
            f"""
            <tr>
              <td>{DATASET_LABELS[row['dataset']]}</td>
              <td>{MODEL_LABELS[row['model']]}</td>
              <td>{format_metric(row['mean_test_mcc'])}</td>
              <td>{format_metric(row['mean_test_precision'])}</td>
              <td>{format_metric(row['mean_test_recall'])}</td>
              <td>{row['mean_fp']:.1f}</td>
              <td>{row['mean_fn']:.1f}</td>
            </tr>
            """
        )
    return "\n".join(rows)


def build_html(selected, summary, assets):
    best = selected.iloc[0]
    top_model = summary.iloc[0]

    html = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>Confusion Matrix Results Slide Deck</title>
  <style>
    :root {{
      --navy: #17324d;
      --muted: #657184;
      --line: #dbe4ee;
      --blue: #3f7fbd;
      --teal: #4aa3a1;
      --green: #75b266;
      --bg: #f6f8fb;
    }}
    body {{
      margin: 0;
      background: #dfe7f1;
      font-family: Arial, Helvetica, sans-serif;
      color: var(--navy);
    }}
    .slide {{
      width: 1280px;
      height: 720px;
      box-sizing: border-box;
      margin: 28px auto;
      padding: 54px 70px 44px 70px;
      background: var(--bg);
      border: 1px solid #c9d4e1;
      position: relative;
      overflow: hidden;
      box-shadow: 0 8px 24px rgba(20, 40, 70, 0.14);
      page-break-after: always;
    }}
    .topbar {{
      position: absolute;
      top: 0;
      left: 0;
      right: 0;
      height: 16px;
      background: #6551d8;
    }}
    h1 {{
      font-size: 42px;
      margin: 0 0 8px 0;
      letter-spacing: 0;
    }}
    h2 {{
      font-size: 32px;
      margin: 0 0 8px 0;
      letter-spacing: 0;
    }}
    .subtitle {{
      color: var(--muted);
      font-size: 18px;
      margin-bottom: 24px;
    }}
    .grid2 {{
      display: grid;
      grid-template-columns: 1fr 1fr;
      gap: 34px;
      align-items: center;
    }}
    .grid3 {{
      display: grid;
      grid-template-columns: repeat(3, 1fr);
      gap: 20px;
    }}
    .card {{
      background: white;
      border: 1.5px solid var(--line);
      border-radius: 8px;
      padding: 20px 22px;
    }}
    .metric {{
      font-size: 42px;
      font-weight: 700;
      margin-bottom: 6px;
    }}
    .metric-label {{
      color: var(--muted);
      font-size: 17px;
    }}
    .callout {{
      background: #eef7f2;
      border-left: 6px solid var(--green);
      padding: 18px 22px;
      font-size: 22px;
      line-height: 1.26;
    }}
    img.chart {{
      width: 100%;
      max-height: 500px;
      object-fit: contain;
      display: block;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 17px;
      background: white;
      border: 1px solid var(--line);
    }}
    th {{
      text-align: left;
      background: #2f7dbd;
      color: white;
      padding: 10px;
      font-weight: 700;
    }}
    td {{
      padding: 9px 10px;
      border-top: 1px solid #e6edf5;
    }}
    tr:nth-child(even) td {{
      background: #f1f5fa;
    }}
    ul {{
      font-size: 22px;
      line-height: 1.35;
      padding-left: 24px;
    }}
    li {{
      margin-bottom: 10px;
    }}
    .small {{
      font-size: 15px;
      color: var(--muted);
    }}
    .footer {{
      position: absolute;
      left: 70px;
      right: 70px;
      bottom: 18px;
      color: #6d7989;
      font-size: 12px;
      border-top: 1px solid #cbd7e5;
      padding-top: 8px;
    }}
    @media print {{
      body {{ background: white; }}
      .slide {{
        margin: 0;
        box-shadow: none;
        border: none;
        page-break-after: always;
      }}
    }}
  </style>
</head>
<body>

<section class="slide">
  <div class="topbar"></div>
  <h1>Confusion Matrix Results</h1>
  <div class="subtitle">Held-out test performance across 3 feature datasets, 10 stratified splits, and 6 model families</div>
  <div class="grid3" style="margin-top: 48px;">
    <div class="card">
      <div class="metric">{format_metric(best['mean_test_mcc'])}</div>
      <div class="metric-label">best average test MCC<br>{DATASET_LABELS[best['dataset']]}</div>
    </div>
    <div class="card">
      <div class="metric">{format_metric(best['mean_test_accuracy'])}</div>
      <div class="metric-label">average test accuracy<br>validation-selected models</div>
    </div>
    <div class="card">
      <div class="metric">{format_metric(best['mean_test_f1'])}</div>
      <div class="metric-label">average test F1<br>validation-selected models</div>
    </div>
  </div>
  <div class="callout" style="margin-top: 58px;">
    The confusion matrices support the same story as the performance analysis: adding CTD improves prediction, adding DPC improves it further, and SVM RBF is the strongest model family for the richer feature sets.
  </div>
  <div class="footer">Generated from confusion_matrices/*.csv</div>
</section>

<section class="slide">
  <div class="topbar"></div>
  <h2>How To Read These Matrices</h2>
  <div class="subtitle">Each seed has a held-out test set with about 88 negatives and 17 positives</div>
  <div class="grid2">
    <div>
      <table>
        <tr><th></th><th>Predicted negative</th><th>Predicted antithrombotic</th></tr>
        <tr><td><b>True negative</b></td><td>TN</td><td>FP</td></tr>
        <tr><td><b>True antithrombotic</b></td><td>FN</td><td>TP</td></tr>
      </table>
      <ul style="margin-top: 28px;">
        <li><b>False positives</b>: negatives incorrectly flagged as antithrombotic.</li>
        <li><b>False negatives</b>: true antithrombotic peptides missed by the model.</li>
        <li><b>MCC</b> is prioritized because the classes are moderately imbalanced.</li>
      </ul>
    </div>
    <div class="card">
      <h2 style="font-size: 26px;">Evaluation design</h2>
      <ul>
        <li>60 percent training, 20 percent validation, 20 percent test.</li>
        <li>Splits are stratified, preserving the positive/negative ratio.</li>
        <li>Best model per seed is selected by validation MCC.</li>
        <li>The held-out test split is used for the final confusion matrix.</li>
      </ul>
    </div>
  </div>
  <div class="footer">Matrix order: [[TN, FP], [FN, TP]]</div>
</section>

<section class="slide">
  <div class="topbar"></div>
  <h2>Validation-Selected Results</h2>
  <div class="subtitle">Average held-out test performance across 10 seeds</div>
  <img class="chart" src="{assets['selected_mcc'].name}" alt="Average test MCC by dataset">
  <div class="footer">Error bars show standard deviation across random seeds.</div>
</section>

<section class="slide">
  <div class="topbar"></div>
  <h2>Average Confusion Matrices</h2>
  <div class="subtitle">Mean counts for validation-selected models across seeds</div>
  <img class="chart" src="{assets['selected_cm'].name}" alt="Average confusion matrices">
  <div class="footer">DPC has the fewest false positives and the highest true positives among validation-selected models.</div>
</section>

<section class="slide">
  <div class="topbar"></div>
  <h2>Selected-Model Summary Table</h2>
  <div class="subtitle">These are the models chosen by validation MCC in each seed</div>
  <table>
    <tr>
      <th>Feature dataset</th><th>Test MCC</th><th>Accuracy</th><th>F1</th><th>TN</th><th>FP</th><th>FN</th><th>TP</th>
    </tr>
    {table_rows_selected(selected)}
  </table>
  <div class="callout" style="margin-top: 28px;">
    PCP + AAC + CTD + DPC is the best validation-selected dataset: higher MCC, higher F1, fewer false positives, and more true positives.
  </div>
  <div class="footer">Counts are averages per held-out test split.</div>
</section>

<section class="slide">
  <div class="topbar"></div>
  <h2>Model Comparison</h2>
  <div class="subtitle">Average test MCC for each model family and feature dataset</div>
  <img class="chart" src="{assets['model_mcc'].name}" alt="Model MCC comparison">
  <div class="footer">SVM RBF has the best average MCC in all three feature spaces in this confusion-matrix audit.</div>
</section>

<section class="slide">
  <div class="topbar"></div>
  <h2>Error Tradeoff</h2>
  <div class="subtitle">Different model families make different kinds of mistakes</div>
  <img class="chart" src="{assets['error_tradeoff'].name}" alt="False positive false negative tradeoff">
  <div class="footer">Lower-left is best: fewer false positives and fewer false negatives.</div>
</section>

<section class="slide">
  <div class="topbar"></div>
  <h2>Model Selection Stability</h2>
  <div class="subtitle">How often each model was selected by validation MCC across 10 seeds</div>
  <img class="chart" src="{assets['selection_counts'].name}" alt="Model selection counts">
  <div class="footer">SVM RBF is most often selected for the CTD and CTD + DPC feature sets.</div>
</section>

<section class="slide">
  <div class="topbar"></div>
  <h2>Top Model Rows</h2>
  <div class="subtitle">Highest average test MCC combinations</div>
  <table>
    <tr>
      <th>Feature dataset</th><th>Model</th><th>MCC</th><th>Precision</th><th>Recall</th><th>FP</th><th>FN</th>
    </tr>
    {table_rows_top_models(summary)}
  </table>
  <div class="footer">Precision reflects false-positive control; recall reflects how many true positives were recovered.</div>
</section>

<section class="slide">
  <div class="topbar"></div>
  <h2>Conclusion</h2>
  <div class="subtitle">What the confusion matrices add to the story</div>
  <div class="grid2">
    <div class="card">
      <h2 style="font-size: 25px;">Main result</h2>
      <ul>
        <li>DPC gives the strongest overall test performance.</li>
        <li>CTD improves over PCP + AAC alone.</li>
        <li>SVM RBF is the most reliable model family for richer features.</li>
      </ul>
    </div>
    <div class="card">
      <h2 style="font-size: 25px;">Recommended wording</h2>
      <ul>
        <li>Use "validation-selected held-out test results."</li>
        <li>Report MCC plus false positives and false negatives.</li>
        <li>Use one representative matrix on the main slide and keep the full audit in backup.</li>
      </ul>
    </div>
  </div>
  <div class="callout" style="margin-top: 32px;">
    Best concise statement: PCP + AAC + CTD + DPC with SVM RBF gave the strongest balance of sensitivity and specificity, with very low false-positive counts and the highest average MCC.
  </div>
  <div class="footer">Full results: confusion_matrix_results_by_model_seed.csv</div>
</section>

</body>
</html>
"""
    DECK_PATH.write_text(html)


def main():
    ASSET_DIR.mkdir(parents=True, exist_ok=True)

    results = pd.read_csv(OUTDIR / "confusion_matrix_results_by_model_seed.csv")
    summary = pd.read_csv(OUTDIR / "confusion_matrix_summary_by_dataset_model.csv")
    selected = pd.read_csv(OUTDIR / "confusion_matrix_summary_validation_selected.csv")

    dataset_order = ["pcp_aac_ctd_dpc", "pcp_aac_ctd", "pcp_aac"]
    selected["dataset"] = pd.Categorical(selected["dataset"], categories=dataset_order, ordered=True)
    selected = selected.sort_values("dataset")

    summary = summary.sort_values("mean_test_mcc", ascending=False)

    assets = {
        "selected_mcc": make_selected_performance_chart(selected),
        "selected_cm": make_selected_confusion_chart(selected),
        "model_mcc": make_model_mcc_chart(summary),
        "error_tradeoff": make_error_tradeoff_chart(summary),
        "selection_counts": make_selected_model_counts_chart(results),
    }

    build_html(selected, summary, assets)
    print(f"Saved slide deck: {DECK_PATH}")
    print(f"Saved chart assets: {ASSET_DIR}")


if __name__ == "__main__":
    main()
