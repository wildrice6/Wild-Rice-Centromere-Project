import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import pearsonr, spearmanr
import argparse
import os

# Global settings
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

def run_correlation_analysis(input_file, config_file, output_dir, fig_w, fig_h):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 1. Load data
    df_raw = pd.read_csv(input_file, sep='\t')
    
    # Aggregation: Merge mutation counts by Species and Chromosome
    # Renaming '物种' to 'Species' for internal consistency
    df_raw = df_raw.rename(columns={'物种': 'Species'})
    
    df = df_raw.groupby(['Species', 'Chromosome']).agg({
        'Cen_Count': 'sum',
        'NonCen_Count': 'sum'
    }).reset_index()

    # Save intermediate summary file (for data verification)
    df.to_csv(os.path.join(output_dir, "intermediate_summarized_counts.tsv"), sep='\t', index=False)

    # 2. Load color configuration file (Tab-separated, no header)
    config = pd.read_csv(config_file, sep='\t', header=None, names=['ID', 'Species', 'Color'])
    color_dict = dict(zip(config['Species'], config['Color']))

    # 3. Calculate correlations
    pr, pp = pearsonr(df['NonCen_Count'], df['Cen_Count'])
    sr, sp = spearmanr(df['NonCen_Count'], df['Cen_Count'])

    # Save statistical analysis result file
    with open(os.path.join(output_dir, "correlation_test_results.txt"), 'w', encoding='utf-8') as f:
        f.write("=== SV Correlation Statistical Analysis ===\n\n")
        f.write(f"Pearson Correlation:\n   r = {pr:.3f}\n   p-value = {pp:.3e}\n\n")
        f.write(f"Spearman Correlation:\n   rho = {sr:.3f}\n   p-value = {sp:.3e}\n")

    # 4. Plot preparation: add slight jitter to distinguish overlapping points
    np.random.seed(42)
    jitter_val = 0.22  # Magnitude of jitter
    df['NonCen_Jitter'] = df['NonCen_Count'] + np.random.uniform(-jitter_val, jitter_val, len(df))
    df['Cen_Jitter'] = df['Cen_Count'] + np.random.uniform(-jitter_val, jitter_val, len(df))

    # 5. Start plotting
    fig, axes = plt.subplots(1, 2, figsize=(fig_w, fig_h))
    
    # Titles showing 3 decimal places
    titles = [f"Pearson Correlation\nr = {pr:.3f}, p = {pp:.3f}", 
              f"Spearman Correlation\nρ = {sr:.3f}, p = {sp:.3f}"]
    
    for i, ax in enumerate(axes):
        # Draw regression line (calculated based on original integer data)
        sns.regplot(data=df, x='NonCen_Count', y='Cen_Count', 
                    ax=ax, scatter=False, color='red', line_kws={'linewidth': 0.8, 'zorder': 1})
        
        # Draw scatter points (using jittered data, alpha=1.0)
        # Plot by species sequentially to match colors
        for species in config['Species']:
            sub = df[df['Species'] == species]
            if not sub.empty:
                ax.scatter(sub['NonCen_Jitter'], sub['Cen_Jitter'], 
                           label=species, color=color_dict[species], 
                           alpha=1.0, s=15, edgecolors='none', zorder=2)

        # Style fine-tuning
        ax.set_title(titles[i], fontsize=7)
        ax.set_xlabel('SV Count (Arms)', fontsize=6)
        ax.set_ylabel('SV Count (Cen)', fontsize=6)
        
        # Set tick label font size
        ax.tick_params(axis='both', which='major', labelsize=5)
        ax.grid(True, linestyle=':', alpha=0.5, linewidth=0.5)

    # Extract legend and place it on the far right of the canvas
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='center left', bbox_to_anchor=(1.0, 0.5), 
               fontsize=5, frameon=False, title="Species", title_fontsize=6)

    plt.tight_layout()
    
    # Save results
    pdf_path = os.path.join(output_dir, "SV_Correlation_Final.pdf")
    fig.savefig(pdf_path, format='pdf', bbox_inches='tight', transparent=True)
    
    print(f"Analysis complete!")
    print(f"1. Vector Chart: {pdf_path}")
    print(f"2. Intermediate Summary Table: {output_dir}/intermediate_summarized_counts.tsv")
    print(f"3. Detailed Statistical Report: {output_dir}/correlation_test_results.txt")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Advanced SV Correlation Scatter Plot Analysis (Jitter + Species Coloring)")
    parser.add_argument('-i', '--input', required=True, help="Input SV count file (tsv)")
    parser.add_argument('-c', '--config', required=True, help="Species color configuration file (tsv)")
    parser.add_argument('-o', '--output_dir', required=True, help="Output directory")
    
    # Canvas size parameters
    parser.add_argument('--width', type=float, default=8.0, help="Canvas width (inches), default 8.0")
    parser.add_argument('--height', type=float, default=3.5, help="Canvas height (inches), default 3.5")
    
    args = parser.parse_args()
    run_correlation_analysis(args.input, args.config, args.output_dir, args.width, args.height)