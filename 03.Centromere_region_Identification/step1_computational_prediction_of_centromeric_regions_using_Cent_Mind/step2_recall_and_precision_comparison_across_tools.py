import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde, false_discovery_control
from itertools import combinations

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "pdf.fonttype": 42,
})

# ===================== 1. Read Data =====================
input_file = r"C:\Users\10042\Desktop\source_data_TableA2.xlsx"
df_raw = pd.read_excel(input_file)

test_loci = df_raw[df_raw["dataset"] == "test"]["Modifi_chr"].unique()
df_all = df_raw[df_raw["Modifi_chr"].isin(test_loci)].copy()

df = df_all[
    df_all["Recall"].notna()
    & df_all["Precision"].notna()
    & np.isfinite(df_all["Recall"])
    & np.isfinite(df_all["Precision"])
].copy()

df["F1"] = (
    2 * df["Recall"] * df["Precision"] /
    (df["Recall"] + df["Precision"])
).fillna(0)

display_names = {
    "bigmodel": "Cent-Mind",
    "centier": "CentIER",
    "quartet": "quarTeT",
}

df["Method_Display"] = df["Method"].map(display_names)

ORDER = ["Cent-Mind", "CentIER", "quarTeT"]
COLORS = {
    "Cent-Mind": "#2166AC",
    "CentIER": "#D6604D",
    "quarTeT": "#4DAC26",
}
LIGHT = {
    "Cent-Mind": "#92C5DE",
    "CentIER": "#F4A582",
    "quarTeT": "#B8E186",
}

df = df[df["Method_Display"].isin(ORDER)].copy()

# ===================== 2. Permutation Test + BH =====================
np.random.seed(42)

def permutation_test(data, metric, n_perm=10000):
    results = []

    for m1, m2 in combinations(data["Method_Display"].unique(), 2):
        g1 = data[data["Method_Display"] == m1][metric].values
        g2 = data[data["Method_Display"] == m2][metric].values

        actual_diff = np.mean(g1) - np.mean(g2)
        combined = np.concatenate([g1, g2])
        n1 = len(g1)

        perm_diffs = []
        for _ in range(n_perm):
            perm = np.random.permutation(combined)
            perm_diffs.append(np.mean(perm[:n1]) - np.mean(perm[n1:]))

        perm_diffs = np.array(perm_diffs)

        p_raw = (
            np.sum(np.abs(perm_diffs) >= np.abs(actual_diff)) + 1
        ) / (n_perm + 1)

        results.append({
            "group1": m1,
            "group2": m2,
            "metric": metric,
            "p_raw": p_raw,
        })

    res = pd.DataFrame(results)
    res["p_adj"] = false_discovery_control(res["p_raw"], method="bh")

    res["sig"] = res["p_adj"].apply(
        lambda p: "***" if p < 0.001 else
        ("**" if p < 0.01 else
         ("*" if p < 0.05 else "ns"))
    )

    return res

recall_perm = permutation_test(df, "Recall", n_perm=10000)
precision_perm = permutation_test(df, "Precision", n_perm=10000)

# ===================== 3. Significance Brackets Inside 0-1 =====================
def draw_brackets(ax, sig_df, order):
    pairs_sorted = sorted(
        [(g1, g2) for g1 in order for g2 in order if order.index(g1) < order.index(g2)],
        key=lambda x: order.index(x[1]) - order.index(x[0])
    )

    y_levels = [0.84, 0.91, 0.975]
    tip = 0.018

    for level, (g1, g2) in enumerate(pairs_sorted):
        row = sig_df[
            (sig_df["group1"].isin([g1, g2])) &
            (sig_df["group2"].isin([g1, g2]))
        ]

        if row.empty:
            continue

        sig_label = row.iloc[0]["sig"]
        p_adj_val = row.iloc[0]["p_adj"]

        x1, x2 = order.index(g1), order.index(g2)
        xmid = (x1 + x2) / 2
        y = y_levels[min(level, len(y_levels) - 1)]

        ax.plot(
            [x1, x1, x2, x2],
            [y - tip, y, y, y - tip],
            lw=1.0,
            color="#222222",
            clip_on=True,
        )

        if sig_label == "ns":
            ax.text(
                xmid, y + 0.006, "ns",
                ha="center",
                va="bottom",
                fontsize=7.5,
                color="#666666",
                clip_on=True,
            )
        else:
            ax.text(
                xmid, y + 0.004, sig_label,
                ha="center",
                va="bottom",
                fontsize=11,
                color="#CC0000",
                fontweight="bold",
                clip_on=True,
            )

            p_str = "p<0.001" if p_adj_val < 0.001 else f"p={p_adj_val:.3f}"
            ax.text(
                xmid, y + 0.040, p_str,
                ha="center",
                va="bottom",
                fontsize=6.2,
                color="#555555",
                clip_on=True,
            )

# ===================== 4. Raincloud Plot =====================
def raincloud(ax, data, metric, perm, ylabel, panel_label):
    for i, method in enumerate(ORDER):
        vals = data[data["Method_Display"] == method][metric].values
        c = COLORS[method]
        lc = LIGHT[method]

        if len(vals) > 2 and vals.std() > 0:
            kde = gaussian_kde(vals, bw_method=0.35)
            y_g = np.linspace(0, 1, 300)
            dens = kde(y_g)
            dens = dens / dens.max() * 0.36

            ax.fill_betweenx(
                y_g,
                i + 0.06,
                i + 0.06 + dens,
                color=lc,
                alpha=0.88,
                linewidth=0,
            )
            ax.plot(
                i + 0.06 + dens,
                y_g,
                color=c,
                lw=0.7,
                alpha=0.65,
            )

        np.random.seed(i * 7 + 13)
        jitter = np.random.uniform(-0.20, -0.05, size=len(vals))

        ax.scatter(
            i + jitter,
            vals,
            color=c,
            alpha=0.45,
            s=9,
            linewidths=0,
            zorder=3,
        )

        ax.boxplot(
            vals,
            positions=[i],
            widths=0.11,
            patch_artist=True,
            notch=False,
            manage_ticks=False,
            boxprops=dict(facecolor=lc, color=c, linewidth=1.3),
            medianprops=dict(color=c, linewidth=2.2),
            whiskerprops=dict(color=c, linewidth=1.1),
            capprops=dict(color=c, linewidth=1.1),
            flierprops=dict(marker=""),
        )

        ax.scatter(
            i,
            vals.mean(),
            marker="D",
            color="white",
            edgecolors=c,
            s=24,
            zorder=5,
            linewidths=1.3,
        )

    ax.axhline(
        1.0,
        color="#CC0000",
        lw=0.9,
        ls="--",
        alpha=0.65,
        zorder=1,
    )

    ax.text(
        len(ORDER) - 0.45,
        0.985,
        "Perfect = 1.0",
        fontsize=6.5,
        color="#CC0000",
        ha="right",
        va="top",
        style="italic",
    )

    draw_brackets(ax, perm, ORDER)

    ax.set_xticks(np.arange(len(ORDER)))
    ax.set_xticklabels(ORDER, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9.5)
    ax.set_xlim(-0.58, len(ORDER) - 0.42)

    # Strictly restrict y-axis to 0-1.
    ax.set_ylim(0, 1.0)
    ax.yaxis.set_major_locator(plt.MultipleLocator(0.25))

    ax.tick_params(axis="x", length=0, pad=4)

    ax.text(
        -0.20,
        1.04,
        panel_label,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        va="top",
    )

# ===================== 5. Plot AB =====================
fig, (ax_A, ax_B) = plt.subplots(
    1,
    2,
    figsize=(8, 5.5),
    gridspec_kw={"wspace": 0.40},
)

raincloud(ax_A, df, "Recall", recall_perm, "Recall", "A")
raincloud(ax_B, df, "Precision", precision_perm, "Precision", "B")

legend_handles = [
    mpatches.Patch(facecolor=COLORS[m], label=m, alpha=0.85)
    for m in ORDER
] + [
    Line2D(
        [0], [0],
        color="#CC0000",
        lw=1.2,
        ls="--",
        label="Perfect score (1.0)",
    ),
    Line2D(
        [0], [0],
        marker="D",
        color="w",
        markerfacecolor="white",
        markeredgecolor="#444444",
        markersize=6,
        label="Mean",
    ),
]

fig.legend(
    handles=legend_handles,
    loc="upper center",
    bbox_to_anchor=(0.5, 1.02),
    ncol=5,
    frameon=False,
    fontsize=8.5,
    columnspacing=1.2,
    handlelength=1.5,
)

fig.suptitle(
    "Centromere Prediction Performance: Cent-Mind vs CentIER vs quarTeT\n"
    "Permutation test (10,000 permutations), BH multiple comparison correction",
    fontsize=10,
    fontweight="bold",
    y=1.10,
)

plt.savefig("AB_recall_precision_y0_1_sig_ab.pdf", bbox_inches="tight")
plt.savefig("AB_recall_precision_y0_1_sig_ab.png", dpi=300, bbox_inches="tight")
plt.close()

print("Saved: AB_recall_precision_y0_1_sig_ab.pdf / .png")