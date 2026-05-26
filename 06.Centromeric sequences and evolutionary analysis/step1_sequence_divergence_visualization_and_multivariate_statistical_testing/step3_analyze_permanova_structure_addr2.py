--- START OF FILE 03.analyze_permanova_structure_addR2.py ---

import os
import argparse
import pandas as pd
import numpy as np
from skbio.stats.distance import mantel, DistanceMatrix, permanova
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rcParams

# 1. Environment Configuration (Adobe Illustrator Compatibility Settings)
# --------------------------------------------------------------------------
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
rcParams['font.size'] = 8

def parse_metadata(index_list):
    """
    Parse sample identifiers to extract metadata.
    Example format: AA_Osat_jap_hap1_Chr01
    """
    metadata = []
    for name in index_list:
        parts = name.split('_')
        n_p = len(parts)
        chrom = parts[-1]
        haplo = parts[-2]
        species = "_".join(parts[:-2])
        genome = parts[0]
        metadata.append({
            'ID': name,
            'Genome': genome,
            'Species': species,
            'Haplotype': haplo,
            'Chromosome': chrom
        })
    return pd.DataFrame(metadata).set_index('ID')

def calculate_r2(f_val, n, k):
    """
    PERMANOVA R2 Calculation Formula: R2 = [F * (k - 1)] / [F * (k - 1) + (n - k)]
    """
    if f_val is None or np.isnan(f_val) or (n - k) <= 0:
        return 0
    num = f_val * (k - 1)
    den = num + (n - k)
    return num / den

def get_sig_label(f_val, p_val):
    """Generate Pseudo-F value label with significance asterisks."""
    if p_val <= 0.001:
        sig = "***"
    elif p_val <= 0.01:
        sig = "**"
    elif p_val <= 0.05:
        sig = "*"
    else:
        sig = "ns"
    return f"{f_val:.1f}{sig}"

def run_permanova_logic(dm, meta, group_name, factors):
    """Execute statistical analysis logic for given factors."""
    results = []
    for factor in factors:
        if meta[factor].nunique() < 2:
            print(f"  Warning: Skipping {factor} in {group_name} (only one group present).")
            continue
            
        print(f"  Testing {group_name} - {factor}...")
        # skbio permanova requires a DistanceMatrix object as input
        res = permanova(dm, meta, column=factor, permutations=999)
        
        f_val = res['test statistic']
        p_val = res['p-value']
        n = res['sample size']
        k = meta[factor].nunique()
        
        r2 = calculate_r2(f_val, n, k)
        label = get_sig_label(f_val, p_val)
        
        results.append({
            'Group': group_name,
            'Factor': factor,
            'Sample_Size_(N)': n,
            'Num_Groups_(k)': k,
            'Pseudo-F': f_val,
            'p-value': p_val,
            'R2': r2,
            'Label': label
        })
    return results

def main():
    parser = argparse.ArgumentParser(description='PERMANOVA analysis with Haplotype factors on a global scale.')
    parser.add_argument('input_matrix', help='Path to input Mash distance matrix (TSV)')
    parser.add_argument('output_dir', help='Path to output directory')
    parser.add_argument('--width', type=float, default=6.0)
    parser.add_argument('--height', type=float, default=5.0)
    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # 1. Load Data
    print(f"Loading matrix: {args.input_matrix}")
    df_matrix = pd.read_csv(args.input_matrix, sep='\t', index_col=0)
    
    # Ensure matrix is C-contiguous to prevent skbio input errors
    data_all_raw = np.ascontiguousarray(df_matrix.values, dtype=float)
    dm_all = DistanceMatrix(data_all_raw, ids=df_matrix.index.tolist())
    
    # 2. Parse Metadata
    meta_all = parse_metadata(df_matrix.index)
    
    # 3. Statistical Testing
    print("Starting PERMANOVA calculations...")
    all_results = []
    
    # 3.1 Global Analysis (All Genomes) - Includes Haplotype
    factors_all = ['Genome', 'Species', 'Chromosome', 'Haplotype']
    all_results.extend(run_permanova_logic(dm_all, meta_all, "All Genomes", factors_all))
    
    # 3.2 Local Analysis (AA Genome Only)
    aa_samples = meta_all[meta_all['Genome'] == 'AA'].index
    if not aa_samples.empty:
        df_aa = df_matrix.loc[aa_samples, aa_samples]
        data_aa_raw = np.ascontiguousarray(df_aa.values, dtype=float)
        dm_aa = DistanceMatrix(data_aa_raw, ids=df_aa.index.tolist())
        meta_aa = meta_all.loc[aa_samples]
        factors_aa = ['Species', 'Chromosome', 'Haplotype']
        all_results.extend(run_permanova_logic(dm_aa, meta_aa, "AA Genome Only", factors_aa))
    
    # Save statistical results
    res_df = pd.DataFrame(all_results)
    stats_out = os.path.join(args.output_dir, 'permanova_stats_summary.txt')
    res_df.to_csv(stats_out, sep='\t', index=False)
    print(f"Stats result saved to {stats_out}")

    # 4. Plotting Logic (Broken Axis Barplot Style)
    # --------------------------------------------------------------------------
    print("Generating Broken Axis Barplot...")
    factor_order = ['Genome', 'Species', 'Chromosome', 'Haplotype']
    colors = {'Genome': '#66c2a5', 'Species': '#fc8d62', 'Chromosome': '#8da0cb', 'Haplotype': '#e78ac3'}

    # Define broken axis thresholds based on data distribution (gap between 12 and 120)
    low_limit = 12
    high_limit = 120
    y_max = res_df['Pseudo-F'].max() * 1.1

    # Create subplots with a height ratio of 1:2.5
    fig, (ax_top, ax_bot) = plt.subplots(2, 1, sharex=True, 
                                            gridspec_kw={'height_ratios': [1, 2.5]},
                                            figsize=(args.width, args.height))
    fig.subplots_adjust(hspace=0.1)

    # Draw barplots on both axes
    for ax in [ax_top, ax_bot]:
        sns.barplot(data=res_df, x='Group', y='Pseudo-F', hue='Factor', 
                    hue_order=factor_order, palette=colors, ax=ax, 
                    edgecolor='black', linewidth=0.8)

    # Adjust axis limits to create the "broken" effect
    ax_top.set_ylim(high_limit, y_max)
    ax_bot.set_ylim(0, low_limit)

    # Hide spines between the top and bottom plots
    ax_top.spines['bottom'].set_visible(False)
    ax_bot.spines['top'].set_visible(False)
    ax_top.tick_params(labelbottom=False, bottom=False)

    # Add break diagonal marks (//)
    d = .015 
    kwargs = dict(transform=ax_top.transAxes, color='k', clip_on=False, lw=1)
    ax_top.plot((-d, +d), (-d, +d), **kwargs)        
    ax_top.plot((1 - d, 1 + d), (-d, +d), **kwargs)  
    kwargs.update(transform=ax_bot.transAxes) 
    ax_bot.plot((-d, +d), (1 - d, 1 + d), **kwargs)  
    ax_bot.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

    # Place significance labels precisely
    groups = ['All Genomes', 'AA Genome Only']
    n_factors = len(factor_order)
    bar_width = 0.8 / n_factors # Default Seaborn bar width

    for g_idx, group in enumerate(groups):
        g_data = res_df[res_df['Group'] == group]
        for f_idx, factor in enumerate(factor_order):
            row = g_data[g_data['Factor'] == factor]
            if not row.empty:
                f_val = row['Pseudo-F'].values[0]
                label_text = row['Label'].values[0]
                # Calculate the center X-coordinate for the bar
                x_pos = g_idx - 0.4 + (f_idx + 0.5) * bar_width
                
                # Determine which subplot to annotate based on the value
                if f_val >= high_limit:
                    ax_top.text(x_pos, f_val + 2, label_text, ha='center', va='bottom', fontsize=7)
                else:
                    ax_bot.text(x_pos, f_val + 0.1, label_text, ha='center', va='bottom', fontsize=7)

    # Legend and Label Beautification
    ax_top.legend(title='Evolutionary Factor', loc='upper right', frameon=False, fontsize=8)
    if ax_bot.get_legend():
        ax_bot.get_legend().remove()
    
    ax_top.set_ylabel('')
    ax_bot.set_ylabel('')
    ax_bot.set_xlabel('')
    fig.text(0.04, 0.5, 'Pseudo-F Statistic (Effect Size)', va='center', rotation='vertical', fontsize=9)

    # Save final figure
    pdf_path = os.path.join(args.output_dir, 'permanova_effect_size_final.pdf')
    plt.savefig(pdf_path, bbox_inches='tight')
    print(f"Success! Plot saved to {pdf_path}")

if __name__ == "__main__":
    main()