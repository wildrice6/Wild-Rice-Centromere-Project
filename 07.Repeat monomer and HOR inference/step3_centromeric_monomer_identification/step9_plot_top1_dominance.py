import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import argparse
import re
import sys

def plot_top1_dominance(data_file, config_file, output_file, width=None, height=None):
    """
    Generates a horizontal matrix bar plot based on Top-1 repeat dominance data and a configuration file.
    The visual style accurately mimics the specified template.
    Supports custom figure dimensions.
    """
    # --- 1. Configure Matplotlib parameters ---
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 8
    plt.rcParams['figure.autolayout'] = False

    # --- 2. Load and preprocess data ---
    try:
        config_df = pd.read_csv(config_file, sep='\t')
        config_df.columns = ['Species', 'Color', 'Species_name']
        config_df = config_df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
        data_df = pd.read_csv(data_file, sep='\t')
    except Exception as e:
        print(f"Error: Failed to read file. Please check the file path and format (should be tab-separated). - {e}", file=sys.stderr)
        sys.exit(1)

    # Data processing logic
    def safe_split(material_string):
        if '_' in material_string: 
            return material_string.rsplit('_', 1)
        else: 
            return [material_string, 'N/A']
            
    split_data = data_df['Material'].apply(safe_split).tolist()
    data_df[['Species', 'Haplotype']] = pd.DataFrame(split_data, index=data_df.index)
    merged_df = pd.merge(data_df, config_df, on='Species', how='left')
    merged_df.dropna(subset=['Species_name'], inplace=True)
    
    if merged_df.empty:
        print("Error: No valid data available for plotting.", file=sys.stderr)
        sys.exit(1)

    # --- 3. Determine plotting order ---
    species_order = config_df['Species'].tolist()
    merged_df['Species'] = pd.Categorical(merged_df['Species'], categories=species_order, ordered=True)
    merged_df.sort_values(by=['Species', 'Haplotype'], inplace=True)
    material_order = merged_df['Material'].unique().tolist()
    
    y_labels_map = {}
    for _, row in merged_df.drop_duplicates(subset=['Material']).iterrows():
        hap_label = row['Haplotype'] if row['Haplotype'] != 'N/A' else ''
        # Format species name to italics
        species_name_formatted = ' '.join([f'$\\it{{{part}}}$' for part in row['Species_name'].split()])
        y_labels_map[row['Material']] = f"{species_name_formatted} {hap_label}".strip()
    
    merged_df['Chr_num'] = merged_df['Chromosome'].apply(lambda x: int(re.search(r'\d+', str(x)).group()))
    chromosome_order = sorted(merged_df['Chromosome'].unique(), key=lambda x: int(re.search(r'\d+', str(x)).group()))

    # --- 4. Create plotting grid ---
    n_rows = len(material_order)
    n_cols = len(chromosome_order)
    
    dynamic_width = max(8, n_cols * 0.8)
    dynamic_height = max(5, n_rows * 0.25)
    fig_width = width if width is not None else dynamic_width
    fig_height = height if height is not None else dynamic_height
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_width, fig_height), 
                             sharex=True, sharey=True, gridspec_kw={'hspace': 0.1, 'wspace': 0.2})
    
    if n_rows == 1 and n_cols == 1: axes = [[axes]]
    elif n_rows == 1: axes = [axes]
    elif n_cols == 1: axes = [[ax] for ax in axes]

    color_low_ows = '#C44E52'  # Red (<= 500)
    color_high_ows = '#4C72B0' # Blue (> 500)

    # --- 5. Iterate through data for plotting ---
    for i, material in enumerate(material_order):
        for j, chromosome in enumerate(chromosome_order):
            ax = axes[i][j]
            point_data = merged_df[(merged_df['Material'] == material) & (merged_df['Chromosome'] == chromosome)]
            
            if not point_data.empty:
                dominance = point_data['Dominance_at_OWS_Percent'].iloc[0]
                ows = point_data['OWS_Final'].iloc[0]
                bar_color = color_high_ows if ows > 500 else color_low_ows
                ax.barh(0, dominance, color=bar_color, height=0.75, zorder=3)

            # --- 6. Format each subplot (ax) ---
            ax.set_xlim(0, 100)
            ax.set_yticks([])
            ax.set_facecolor('none')
            for spine in ax.spines.values(): 
                spine.set_visible(False)
            
            # Modification 3: Draw a dotted vertical line at 60%
            ax.axvline(x=60, color='grey', linestyle=':', linewidth=0.8, zorder=0)

            if j == 0: 
                ax.set_ylabel(y_labels_map[material], rotation=0, ha='right', va='center', fontsize=8, labelpad=5)
            if i == 0: 
                ax.set_title(chromosome, fontsize=9, pad=10)

            # Modifications 1 & 2: Adjust X-axis tick style for the bottom row
            if i == n_rows - 1:
                # Draw a clear black baseline (from 0 to 100)
                ax.spines['bottom'].set_visible(True)
                ax.spines['bottom'].set_color('black')
                ax.spines['bottom'].set_linewidth(0.8)
                
                # Set tick positions, including 100, to ensure markers at the end of the baseline
                ax.set_xticks([0, 40, 80, 100])
                # Provide labels only for 0, 40, 80; label for 100 is an empty string
                ax.set_xticklabels(['0', '40', '80', ''])
                ax.tick_params(axis='x', which='both', bottom=True, labelbottom=True, labelsize=8, length=2, pad=2)
            else:
                ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)

    # --- 7. Add global elements ---
    fig.suptitle('Dominance of Top-1 Repeat in Centromeres', fontsize=12, y=0.99)
    legend_elements = [
        Patch(facecolor=color_low_ows, label='Optimal Window Size <= 500'),
        Patch(facecolor=color_high_ows, label='Optimal Window Size > 500')]
    fig.legend(handles=legend_elements, loc='lower right', bbox_to_anchor=(0.95, 0.01), fontsize=8, frameon=True, edgecolor='black')
    
    fig.tight_layout(rect=[0.01, 0.05, 0.99, 0.95])

    # --- 8. Save the figure ---
    try:
        plt.savefig(output_file, format='pdf', bbox_inches='tight', pad_inches=0.1, dpi=300)
        print(f"Plotting successful! Result saved to: {output_file}")
    except Exception as e:
        print(f"Error: An error occurred while saving the file - {e}", file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description="Plots a horizontal matrix bar chart of Top-1 repeat dominance, styled to match specific templates.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('data_file', type=str, help="Dominance result file for each material (tab-separated).")
    parser.add_argument('config_file', type=str, help="Configuration file for plotting (tab-separated).")
    parser.add_argument('output_file', type=str, help="Output PDF filename.")
    parser.add_argument('--width', type=float, default=None, help="Optional: Set the width of the output image (inches).")
    parser.add_argument('--height', type=float, default=None, help="Optional: Set the height of the output image (inches).")
    
    args = parser.parse_args()
    
    plot_top1_dominance(args.data_file, args.config_file, args.output_file, args.width, args.height)

if __name__ == '__main__':
    main()