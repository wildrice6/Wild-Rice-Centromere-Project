import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
from matplotlib.gridspec import GridSpec

# ===================== Global Configuration =====================
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['figure.facecolor'] = '#ffffff'
plt.rcParams['axes.facecolor'] = '#ffffff'
plt.rcParams['font.size'] = 8
plt.rcParams['axes.titlesize'] = 10
plt.rcParams['axes.labelsize'] = 8
plt.rcParams['xtick.labelsize'] = 7
plt.rcParams['ytick.labelsize'] = 7
plt.rcParams['legend.fontsize'] = 7
plt.rcParams['lines.linewidth'] = 0.8
plt.rcParams['axes.linewidth'] = 0.4
plt.rcParams['patch.linewidth'] = 0.4

# Core Species Configuration
SPECIES_CONFIG = {
    "AA_Osat_jap": {"display_name": "O. sativa ssp. japonica"},
    "AA_Osat_ind": {"display_name": "O. sativa ssp. indica"},
    "AA_Ogla": {"display_name": "O. glaberrima"},
    "AA_Oruf": {"display_name": "O. rufipogon"},
    "AA_Oniv": {"display_name": "O. nivara"},
    "AA_Olon": {"display_name": "O. longistaminata"},
    "AA_Oglu": {"display_name": "O. glumaepatula"},
    "BB_Opun": {"display_name": "O. punctata"},
    "CC_Ooff": {"display_name": "O. officinalis"},
    "EE_Oaus": {"display_name": "O. australiensis"},
    "FF_Obra": {"display_name": "O. brachyantha"},
    "GG_Omey": {"display_name": "O. meyeriana"}
}
SPECIES_ORDER = list(SPECIES_CONFIG.keys())

# Color Mapping: Blue for upper track/baseline, Red for lower track/target
COLOR_MAP_LEFT = {
    "Non_TOP": "#1976D2",
    "Top1_HOR": "#D32F2F"
}
COLOR_MAP_RIGHT = {
    "Bottom 20%": "#1976D2",
    "Top 20%": "#D32F2F"
}

TRACK_OFFSET = 0.15
BOX_WIDTH = 0.18
BOX_OFFSET = 0.18
LEFT_XLIM = (-4, 4)
RIGHT_XLIM = (-4, 4)
N_PERM = 10000

# ===================== Path Configuration =====================
LEFT_DATA_DIR = "/share/org/YZWL/yzwl_shahd/wangdan/Oryza_Centromere/HOR_CENH3/01.process_centromere_signals/Final_Results"
RIGHT_DATA_DIR = "/share/org/YZWL/yzwl_shahd/wangdan/Oryza_Centromere/HOR_CENH3/09.HOR_score_top_bottom0.2/Intermediate_Data"
OUTPUT_DIR = "./Result_Output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ===================== Utility Functions =====================
def extract_chromosome(chrom_str):
    if pd.isna(chrom_str):
        return None
    match = re.search(r'(\d+)', str(chrom_str))
    return int(match.group(1)) if match else None

def get_sig_star(p):
    if pd.isna(p):
        return "ns"
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "ns"

def permutation_test(group1, group2, n_permutations=10000):
    group1 = group1[~np.isnan(group1)]
    group2 = group2[~np.isnan(group2)]
    if len(group1) < 2 or len(group2) < 2:
        return np.nan
    min_len = min(len(group1), len(group2))
    group1 = group1[:min_len]
    group2 = group2[:min_len]
    obs_diff = np.mean(group1) - np.mean(group2)
    combined = np.concatenate([group1, group2])
    perm_diffs = []
    for _ in range(n_permutations):
        np.random.shuffle(combined)
        perm1 = combined[:len(group1)]
        perm2 = combined[len(group1):]
        perm_diffs.append(np.mean(perm1) - np.mean(perm2))
    perm_diffs = np.array(perm_diffs)
    p_value = np.mean(np.abs(perm_diffs) >= np.abs(obs_diff))
    return p_value

def get_species_from_filename(filename):
    fname = os.path.splitext(filename)[0]
    for sp in SPECIES_ORDER:
        if sp in fname:
            return sp
    return None

# ===================== Data Loading =====================
def load_all_data():
    left_data = []
    left_files = glob.glob(os.path.join(LEFT_DATA_DIR, "*.bed"))
    print(f"🔍 Loading Left Plot data: Found {len(left_files)} BED files")
    for f in left_files:
        fname = os.path.basename(f)
        species = get_species_from_filename(fname)
        if not species:
            print(f"⚠️ Skipping {fname}: No matching species configuration")
            continue
        try:
            df = pd.read_csv(f, sep='\t', header=0, on_bad_lines='skip')
            required_cols = ['Chr', 'Region_Type', 'CENH3_Mean_log2Ratio']
            missing = [c for c in required_cols if c not in df.columns]
            if missing:
                print(f"⚠️ Skipping {fname}, missing columns: {missing}")
                continue
            df['CENH3_Mean_log2Ratio'] = pd.to_numeric(df['CENH3_Mean_log2Ratio'], errors='coerce')
            df = df.dropna(subset=['CENH3_Mean_log2Ratio'])
            df['chrom_num'] = df['Chr'].apply(extract_chromosome)
            df = df[df['chrom_num'].between(1, 12)]
            df = df[df['Region_Type'].isin(['Top1_HOR', 'Non_HOR'])]
            df['Region_Type'] = df['Region_Type'].replace({'Non_HOR': 'Non_TOP'})
            df['Species'] = species
            df = df.rename(columns={'CENH3_Mean_log2Ratio': 'Value'})
            left_data.append(df[['Species', 'Region_Type', 'Value']])
        except Exception as e:
            print(f"⚠️ Error processing {fname}: {e}")
            continue
    left_df = pd.concat(left_data, ignore_index=True) if left_data else pd.DataFrame(columns=['Species', 'Region_Type', 'Value'])

    right_data = []
    right_files = glob.glob(os.path.join(RIGHT_DATA_DIR, "*.tsv"))
    print(f"\n🔍 Loading Right Plot data: Found {len(right_files)} TSV files")
    for fp in right_files:
        fname = os.path.basename(fp)
        species = get_species_from_filename(fname)
        if not species:
            print(f"⚠️ Skipping {fname}: No matching species configuration")
            continue
        try:
            df = pd.read_csv(fp, sep='\t', header=0, on_bad_lines='skip')
            required_cols = ['Chr', 'CENH3_Mean_log2Ratio', 'Group']
            missing = [c for c in required_cols if c not in df.columns]
            if missing:
                print(f"⚠️ Skipping {fname}, missing columns: {missing}")
                continue
            df = df[df['Group'].isin(['Top 20%', 'Bottom 20%'])]
            if df.empty:
                print(f"⚠️ Skipping {fname}: No relevant grouped data")
                continue
            df['CENH3_Mean_log2Ratio'] = pd.to_numeric(df['CENH3_Mean_log2Ratio'], errors='coerce')
            df = df.dropna(subset=['CENH3_Mean_log2Ratio'])
            df['chrom_num'] = df['Chr'].apply(extract_chromosome)
            df = df[df['chrom_num'].between(1, 12)]
            df['Species'] = species
            df = df.rename(columns={'CENH3_Mean_log2Ratio': 'Value'})
            right_data.append(df[['Species', 'Group', 'Value']])
        except Exception as e:
            print(f"⚠️ Error processing {fname}: {e}")
            continue
    right_df = pd.concat(right_data, ignore_index=True) if right_data else pd.DataFrame(columns=['Species', 'Group', 'Value'])

    if not left_df.empty:
        left_agg = left_df.groupby(['Species', 'Region_Type'], dropna=False)['Value'].agg(
            mean='mean', std='std', min='min', max='max', n='count'
        ).reset_index()
        left_agg['se'] = left_agg['std'] / np.sqrt(left_agg['n'])
    else:
        left_agg = pd.DataFrame(columns=['Species', 'Region_Type', 'mean', 'std', 'min', 'max', 'n', 'se'])

    if not right_df.empty:
        right_agg = right_df.groupby(['Species', 'Group'], dropna=False)['Value'].agg(
            mean='mean', median='median', std='std', n='count'
        ).reset_index()
    else:
        right_agg = pd.DataFrame(columns=['Species', 'Group', 'mean', 'median', 'std', 'n'])

    return left_df, right_df, left_agg, right_agg

# ===================== Significance Calculation =====================
def calculate_significance(left_df, right_df):
    left_sig = []
    print("\n=== Left Plot Significance (Top1_HOR vs Non_TOP) ===")
    for sp in SPECIES_ORDER:
        d = left_df[left_df['Species'] == sp]
        if len(d['Region_Type'].unique()) < 2:
            left_sig.append({'Species': sp, 'p_value': np.nan, 'significance': 'ns'})
            print(f"👉 {sp} | Insufficient data | ns")
            continue
        g1 = d[d['Region_Type'] == 'Top1_HOR']['Value'].values
        g2 = d[d['Region_Type'] == 'Non_TOP']['Value'].values
        p = permutation_test(g1, g2, N_PERM)
        sig = get_sig_star(p)
        left_sig.append({'Species': sp, 'p_value': p, 'significance': sig})
        print(f"👉 {sp} | p={p:.4f} | {sig}")

    right_sig = []
    print("\n=== Right Plot Significance (Top 20% vs Bottom 20%) ===")
    for sp in SPECIES_ORDER:
        d = right_df[right_df['Species'] == sp]
        if len(d['Group'].unique()) < 2:
            right_sig.append({'Species': sp, 'p_value': np.nan, 'significance': 'ns', 'max_value': np.nan})
            print(f"👉 {sp} | Insufficient data | ns")
            continue
        g1 = d[d['Group'] == 'Top 20%']['Value'].values
        g2 = d[d['Group'] == 'Bottom 20%']['Value'].values
        if len(g1) < 2 or len(g2) < 2:
            right_sig.append({'Species': sp, 'p_value': np.nan, 'significance': 'ns', 'max_value': np.nan})
            continue
        p = permutation_test(g1, g2, N_PERM)
        sig = get_sig_star(p)
        maxv = d['Value'].max()
        right_sig.append({'Species': sp, 'p_value': p, 'significance': sig, 'max_value': maxv})
        print(f"👉 {sp} | p={p:.4f} | {sig}")

    return pd.DataFrame(left_sig), pd.DataFrame(right_sig)

# ===================== Plotting (Final Version: Native Boxplots + Outliers only) =====================
def plot_combined_figure(left_df, right_df, left_agg, right_agg, left_sig, right_sig):
    if left_df.empty and right_df.empty:
        print("⚠️ No data available, skipping plot generation")
        return

    print("\n🎨 Generating Visualization (Native boxplots + outliers shown, blue upper/red lower)...")
    y_labels = [SPECIES_CONFIG[sp]['display_name'] for sp in SPECIES_ORDER]
    y_pos = np.arange(len(y_labels))

    fig = plt.figure(figsize=(18, 10))
    gs = GridSpec(1, 2, width_ratios=[1, 1.5])

    # ===================== Left Plot: Mean & Error Bars =====================
    ax1 = fig.add_subplot(gs[0])
    ax1.set_xlim(LEFT_XLIM)
    for idx, sp in enumerate(SPECIES_ORDER):
        agg = left_agg[left_agg['Species'] == sp]
        
        # Non_TOP (Blue) - plotted with offset
        n = agg[agg['Region_Type'] == 'Non_TOP']
        if not n.empty:
            m = n['mean'].iloc[0]
            ax1.hlines(idx + TRACK_OFFSET, n['min'].iloc[0], n['max'].iloc[0], color=COLOR_MAP_LEFT['Non_TOP'], lw=0.8)
            ax1.scatter(m, idx + TRACK_OFFSET, color=COLOR_MAP_LEFT['Non_TOP'], marker='o', s=30, ec='k')
            ax1.errorbar(m, idx + TRACK_OFFSET, xerr=n['se'].iloc[0], fmt='none', color=COLOR_MAP_LEFT['Non_TOP'], capsize=1)
            ax1.text(m + 0.1, idx + TRACK_OFFSET, f"{m:.2f}", fontsize=6, va='center', c='darkblue')
            
        # Top1_HOR (Red) - plotted with offset
        t = agg[agg['Region_Type'] == 'Top1_HOR']
        if not t.empty:
            m = t['mean'].iloc[0]
            ax1.hlines(idx - TRACK_OFFSET, t['min'].iloc[0], t['max'].iloc[0], color=COLOR_MAP_LEFT['Top1_HOR'], lw=0.8)
            ax1.scatter(m, idx - TRACK_OFFSET, color=COLOR_MAP_LEFT['Top1_HOR'], marker='^', s=30, ec='k')
            ax1.errorbar(m, idx - TRACK_OFFSET, xerr=t['se'].iloc[0], fmt='none', color=COLOR_MAP_LEFT['Top1_HOR'], capsize=1)
            ax1.text(m + 0.1, idx - TRACK_OFFSET, f"{m:.2f}", fontsize=6, va='center', c='darkred')
            
        sig = left_sig[left_sig['Species'] == sp]['significance'].iloc[0]
        ax1.text(LEFT_XLIM[1] - 0.1, idx, sig, ha='right', va='center', fontsize=6, fontweight='bold')

    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(y_labels, fontsize=7)
    ax1.invert_yaxis()
    ax1.set_xlabel('CENH3 log2Ratio (GQ)', fontweight='bold')
    ax1.set_title('Top1_HOR vs Non_TOP', fontweight='bold')
    ax1.grid(axis='x', alpha=0.3)
    ax1.spines[['top', 'right']].set_visible(False)
    ax1.legend([
        plt.Line2D([0],[0],marker='o',c='w',mfc=COLOR_MAP_LEFT['Non_TOP'],ms=5,mec='k'),
        plt.Line2D([0],[0],marker='^',c='w',mfc=COLOR_MAP_LEFT['Top1_HOR'],ms=5,mec='k')
    ], ['Non_TOP', 'Top1_HOR'], loc='upper right', fontsize=6)

    # ===================== Right Plot: Native Boxplots (Outliers only) =====================
    ax2 = fig.add_subplot(gs[1])
    ax2.set_xlim(RIGHT_XLIM)
    for idx, sp in enumerate(SPECIES_ORDER):
        # Blue: Bottom 20% —— Upper position in track
        bot_data = right_df[(right_df['Species']==sp) & (right_df['Group']=='Bottom 20%')]['Value'].values
        if len(bot_data) > 0:
            ax2.boxplot(
                [bot_data], positions=[idx + BOX_OFFSET], widths=BOX_WIDTH,
                patch_artist=True,
                vert=False,
                showfliers=True,      # Displays outliers
                medianprops={'color':'white','lw':1},
                boxprops={'facecolor':COLOR_MAP_RIGHT['Bottom 20%'], 'alpha':0.7, 'linewidth':1},
                flierprops={'marker':'o', 'markerfacecolor':'black', 'markeredgecolor':'black', 'markersize':3}
            )

        # Red: Top 20% —— Lower position in track
        top_data = right_df[(right_df['Species']==sp) & (right_df['Group']=='Top 20%')]['Value'].values
        if len(top_data) > 0:
            ax2.boxplot(
                [top_data], positions=[idx - BOX_OFFSET], widths=BOX_WIDTH,
                patch_artist=True,
                vert=False,
                showfliers=True,      # Displays outliers
                medianprops={'color':'white','lw':1},
                boxprops={'facecolor':COLOR_MAP_RIGHT['Top 20%'], 'alpha':0.7, 'linewidth':1},
                flierprops={'marker':'o', 'markerfacecolor':'black', 'markeredgecolor':'black', 'markersize':3}
            )

        # Statistical Significance Annotation
        rs = right_sig[right_sig['Species'] == sp]
        if not rs.empty and not pd.isna(rs['max_value'].iloc[0]):
            sig = rs['significance'].iloc[0]
            ax2.text(RIGHT_XLIM[1]-0.1, idx, sig, ha='right', va='center', fontsize=6, fontweight='bold')

    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(y_labels, fontsize=7)
    ax2.invert_yaxis()
    ax2.set_xlabel('CENH3 log2Ratio (GQ)', fontweight='bold')
    ax2.set_title('Top 20% vs Bottom 20%', fontweight='bold')
    ax2.grid(axis='x', alpha=0.3)
    ax2.spines[['top', 'right']].set_visible(False)
    ax2.legend([
        plt.Rectangle((0,0),1,1,fc=COLOR_MAP_RIGHT['Bottom 20%'],alpha=0.7),
        plt.Rectangle((0,0),1,1,fc=COLOR_MAP_RIGHT['Top 20%'],alpha=0.7)
    ], ['Bottom 20%', 'Top 20%'], loc='lower right', fontsize=6)

    plt.tight_layout()
    pdf_path = os.path.join(OUTPUT_DIR, "Final_Figure_Outliers.pdf")
    plt.savefig(pdf_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✅ Final visualization (outliers only) saved: {pdf_path}")

# ===================== Main Execution =====================
def main():
    np.random.seed(42)
    left_df, right_df, left_agg, right_agg = load_all_data()
    if left_df.empty and right_df.empty:
        print("❌ Data not found")
        return

    left_sig, right_sig = calculate_significance(left_df, right_df)
    plot_combined_figure(left_df, right_df, left_agg, right_agg, left_sig, right_sig)
    print("\n🎉 Process completed successfully!")

if __name__ == "__main__":
    main()