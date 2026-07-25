#!/usr/bin/env python3

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


OUTDIR = Path("results/figures/pca_plots")
OUTFILE = OUTDIR / "pca_3d_interactive.html"

TABLES = ["data/processed/baseline/pcp_aac.tsv", "data/processed/baseline/pcp_aac_ctd.tsv", "data/processed/baseline/pcp_aac_ctd_dpc.tsv"]
DATASET_LABELS = {
    "pcp_aac": "PCP + AAC",
    "pcp_aac_ctd": "PCP + AAC + CTD",
    "pcp_aac_ctd_dpc": "PCP + AAC + CTD + DPC",
}


def load_table(path):
    df = pd.read_csv(path, sep="\t")
    y = df["label"].astype(int)
    drop_cols = ["label"]
    if "sequence" in df.columns:
        drop_cols.append("sequence")
    X = df.drop(columns=drop_cols)
    X = X.apply(pd.to_numeric, errors="coerce")
    X = X.replace([np.inf, -np.inf], np.nan)
    X = X.fillna(X.median(numeric_only=True))
    return X, y


def build_dataset_payload():
    payload = {}
    for table in TABLES:
        dataset = Path(table).stem
        X, y = load_table(table)
        scaled = StandardScaler().fit_transform(X)
        pca = PCA(n_components=3, random_state=0)
        scores = pca.fit_transform(scaled)

        # Normalize each axis for consistent rendering in the browser.
        max_abs = np.abs(scores).max(axis=0)
        max_abs[max_abs == 0] = 1
        normalized = scores / max_abs

        payload[dataset] = {
            "label": DATASET_LABELS[dataset],
            "variance": [float(v) for v in pca.explained_variance_ratio_],
            "points": [
                {
                    "x": float(normalized[i, 0]),
                    "y": float(normalized[i, 1]),
                    "z": float(normalized[i, 2]),
                    "label": int(y.iloc[i]),
                }
                for i in range(len(y))
            ],
        }
    return payload


def write_html(payload):
    html = f"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>Interactive 3D PCA Plot</title>
<style>
body {{
  margin: 0;
  background: #dfe7f1;
  font-family: Arial, Helvetica, sans-serif;
  color: #17324d;
}}
.wrap {{
  width: 1180px;
  margin: 28px auto;
  background: #f6f8fb;
  border: 1px solid #c9d4e1;
  box-shadow: 0 8px 24px rgba(20,40,70,.14);
  padding: 44px 54px 38px;
  box-sizing: border-box;
}}
h1 {{
  margin: 0 0 6px;
  font-size: 38px;
}}
.subtitle {{
  color: #657184;
  font-size: 18px;
  margin-bottom: 22px;
}}
.toolbar {{
  display: flex;
  align-items: center;
  gap: 10px;
  margin-bottom: 18px;
}}
button {{
  border: 1px solid #b9c7d8;
  background: white;
  color: #17324d;
  padding: 9px 13px;
  border-radius: 7px;
  font-size: 15px;
  cursor: pointer;
}}
button.active {{
  background: #2f7dbd;
  color: white;
  border-color: #2f7dbd;
}}
.panel {{
  display: grid;
  grid-template-columns: 780px 1fr;
  gap: 28px;
  align-items: start;
}}
canvas {{
  width: 780px;
  height: 540px;
  background: white;
  border: 1px solid #dbe4ee;
  border-radius: 8px;
}}
.card {{
  background: white;
  border: 1px solid #dbe4ee;
  border-radius: 8px;
  padding: 18px 20px;
}}
.metric {{
  font-size: 28px;
  font-weight: 700;
  margin: 8px 0;
}}
.small {{
  color: #657184;
  font-size: 14px;
  line-height: 1.35;
}}
.legend-row {{
  display: flex;
  align-items: center;
  gap: 8px;
  margin: 10px 0;
}}
.dot {{
  width: 12px;
  height: 12px;
  border-radius: 50%;
  display: inline-block;
}}
</style>
</head>
<body>
<div class="wrap">
  <h1>Interactive 3D PCA Plot</h1>
  <div class="subtitle">PC1, PC2, and PC3 view of each peptide feature dataset. Drag the plot to rotate.</div>

  <div class="toolbar" id="buttons"></div>

  <div class="panel">
    <canvas id="plot" width="780" height="540"></canvas>
    <div class="card">
      <h2 id="datasetTitle" style="margin-top:0;"></h2>
      <div class="small">Explained variance</div>
      <div class="metric" id="varianceText"></div>
      <div class="legend-row"><span class="dot" style="background:#4a6f9f;"></span> Negative peptides</div>
      <div class="legend-row"><span class="dot" style="background:#d95f5f;"></span> Antithrombotic positives</div>
      <p class="small">
        PCA is unsupervised. It compresses all features into three combined axes so you can visually inspect whether positives and negatives occupy different regions of feature space.
      </p>
      <p class="small">
        Direction by itself is not biologically meaningful; clustering and separation are the useful parts.
      </p>
      <button id="resetBtn">Reset view</button>
    </div>
  </div>
</div>

<script>
const DATASETS = {json.dumps(payload)};
const canvas = document.getElementById("plot");
const ctx = canvas.getContext("2d");
const buttons = document.getElementById("buttons");
const title = document.getElementById("datasetTitle");
const varianceText = document.getElementById("varianceText");
const resetBtn = document.getElementById("resetBtn");

let currentKey = Object.keys(DATASETS)[0];
let angleX = -0.55;
let angleY = 0.72;
let dragging = false;
let lastX = 0;
let lastY = 0;

function resetView() {{
  angleX = -0.55;
  angleY = 0.72;
  draw();
}}

function rotatePoint(p) {{
  const cosY = Math.cos(angleY), sinY = Math.sin(angleY);
  const cosX = Math.cos(angleX), sinX = Math.sin(angleX);

  let x = p.x * cosY + p.z * sinY;
  let z = -p.x * sinY + p.z * cosY;
  let y = p.y * cosX - z * sinX;
  z = p.y * sinX + z * cosX;
  return {{x, y, z}};
}}

function project(p) {{
  const scale = 205;
  return {{
    x: canvas.width / 2 + p.x * scale,
    y: canvas.height / 2 - p.y * scale,
    z: p.z
  }};
}}

function drawAxes() {{
  const axes = [
    {{name:"PC1", x:1.18, y:0, z:0, color:"#2f7dbd"}},
    {{name:"PC2", x:0, y:1.18, z:0, color:"#6a55d8"}},
    {{name:"PC3", x:0, y:0, z:1.18, color:"#75b266"}}
  ];
  ctx.lineWidth = 2;
  ctx.font = "14px Arial";
  for (const axis of axes) {{
    const start = project(rotatePoint({{x:0,y:0,z:0}}));
    const end = project(rotatePoint(axis));
    ctx.strokeStyle = axis.color;
    ctx.beginPath();
    ctx.moveTo(start.x, start.y);
    ctx.lineTo(end.x, end.y);
    ctx.stroke();
    ctx.fillStyle = axis.color;
    ctx.fillText(axis.name, end.x + 6, end.y + 4);
  }}
}}

function draw() {{
  const data = DATASETS[currentKey];
  ctx.clearRect(0, 0, canvas.width, canvas.height);
  ctx.fillStyle = "#ffffff";
  ctx.fillRect(0, 0, canvas.width, canvas.height);

  drawAxes();

  const projected = data.points.map(p => {{
    const rotated = rotatePoint(p);
    const screen = project(rotated);
    return {{...p, sx: screen.x, sy: screen.y, depth: screen.z}};
  }}).sort((a, b) => a.depth - b.depth);

  for (const p of projected) {{
    const isPositive = p.label === 1;
    const radius = isPositive ? 5.2 : 3.5;
    const alpha = 0.58 + (p.depth + 1) * 0.16;
    ctx.globalAlpha = Math.max(0.35, Math.min(0.92, alpha));
    ctx.fillStyle = isPositive ? "#d95f5f" : "#4a6f9f";
    ctx.strokeStyle = isPositive ? "#8d3030" : "#29435f";
    ctx.lineWidth = 0.8;
    ctx.beginPath();
    ctx.arc(p.sx, p.sy, radius, 0, Math.PI * 2);
    ctx.fill();
    ctx.stroke();
  }}
  ctx.globalAlpha = 1;
}}

function setDataset(key) {{
  currentKey = key;
  for (const btn of buttons.querySelectorAll("button")) {{
    btn.classList.toggle("active", btn.dataset.key === key);
  }}
  const data = DATASETS[key];
  title.textContent = data.label;
  varianceText.textContent = `PC1 ${{(data.variance[0]*100).toFixed(1)}}% | PC2 ${{(data.variance[1]*100).toFixed(1)}}% | PC3 ${{(data.variance[2]*100).toFixed(1)}}%`;
  draw();
}}

for (const [key, data] of Object.entries(DATASETS)) {{
  const btn = document.createElement("button");
  btn.textContent = data.label;
  btn.dataset.key = key;
  btn.addEventListener("click", () => setDataset(key));
  buttons.appendChild(btn);
}}

canvas.addEventListener("mousedown", e => {{
  dragging = true;
  lastX = e.clientX;
  lastY = e.clientY;
}});
window.addEventListener("mouseup", () => dragging = false);
window.addEventListener("mousemove", e => {{
  if (!dragging) return;
  const dx = e.clientX - lastX;
  const dy = e.clientY - lastY;
  lastX = e.clientX;
  lastY = e.clientY;
  angleY += dx * 0.01;
  angleX += dy * 0.01;
  draw();
}});

resetBtn.addEventListener("click", resetView);
setDataset(currentKey);
</script>
</body>
</html>
"""
    OUTFILE.write_text(html)


def main():
    OUTDIR.mkdir(exist_ok=True)
    payload = build_dataset_payload()
    write_html(payload)
    print(f"Saved interactive 3D PCA plot: {OUTFILE}")


if __name__ == "__main__":
    main()
