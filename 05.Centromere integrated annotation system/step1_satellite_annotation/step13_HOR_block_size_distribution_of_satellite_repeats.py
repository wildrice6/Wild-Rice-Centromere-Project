#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
import matplotlib.pyplot as plt

# -------- 0) Input/Output Configuration --------
INPUT_TSV = "C:/Users/10042/Desktop/CC_EE126155.tsv"
OUT_PDF   = "C:/Users/10042/Desktop/blocksize_plot.pdf"

# -------- 1) Data Loading (Utilizing pre-calculated percentages) --------
df = pd.read_csv(INPUT_TSV, sep="\t")

# -------- 2) Style Settings: hap1 (solid line + circle), hap2 (dashed line + triangle) --------
linestyle_map = {"hap1": "-", "hap2": "--"}
marker_map = {"hap1": "o", "hap2": "^"}

# Fixed color mapping (defined by Species + CEN combination)
color_map = {
    ("O. officinalis", "CEN126"): "#1f78b4",
    ("O. officinalis", "CEN155"): "#fdbf6f",
    ("O. australiensis", "CEN155"): "#e31a1c"
}

# -------- 3) Visualization --------
fig, ax = plt.subplots(figsize=(10, 6))

for (species, hap, cen), group in df.groupby(["species", "hap", "cen"]):
    group = group.sort_values("blocksize")
    x = group["blocksize"]
    y = group["percentage(%)"]

    color = color_map.get((species, cen), "#999")
    ls = linestyle_map[hap]
    mk = marker_map[hap]

    ax.plot(
        x, y,
        color=color,
        linestyle=ls,
        marker=mk,
        linewidth=2,
        markersize=4,
        label=f"{species} {cen} {hap}"
    )

# -------- 4) Axis Specifications: X-axis starts from 3, Y-axis max at 20 --------
ax.set_xlabel("HOR block monomer length", fontsize=12)
ax.set_ylabel("% of HORs", fontsize=12)
ax.set_xlim(3, 50)          # Set X-axis lower limit to 3
ax.set_xticks(range(0, 51, 10))
ax.set_ylim(0, 20)          # Set Y-axis upper limit to 20
ax.set_yticks(range(0, 21, 4))
ax.grid(alpha=0.3)
ax.legend(fontsize=9, frameon=False)

plt.tight_layout()
plt.savefig(OUT_PDF, dpi=300)
print(f"✅ Success: Plot has been saved to {OUT_PDF}")