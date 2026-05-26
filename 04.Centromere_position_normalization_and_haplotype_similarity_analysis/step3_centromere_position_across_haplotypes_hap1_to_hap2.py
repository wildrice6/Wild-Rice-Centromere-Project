import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

# Set Matplotlib style for better aesthetic default effects
plt.style.use('seaborn-v0_8-whitegrid')

def parse_arguments():
    """Parse command-line arguments"""
    parser = argparse.ArgumentParser(
        description="Plot normalized centromere positions faceted by chromosome, colored by haplotype."
    )
    parser.add_argument('--input', required=True, help="Path to the input TSV file with normalized centromere data.")
    parser.add_argument('--output', required=True, help="Path for the output plot (e.g., faceted_centromere_plot.png).")
    return parser.parse_args()

# Define fixed species order and corresponding labels
Y_AXIS_ORDER = [
    'AA_Ogla', 'AA_Oruf', 'AA_Oniv', 'AA_Olon', 'AA_Oglu',
    'BB_Opun', 'CC_Ooff', 'EE_Oaus', 'FF_Obra', 'GG_Omey'
]

Y_AXIS_LABELS = {
    'AA_Ogla': 'O. glaberrima (AA)',
    'AA_Oruf': 'O. rufipogon (AA)',
    'AA_Oniv': 'O. nivara (AA)',
    'AA_Olon': 'O. longistaminata (AA)',
    'AA_Oglu': 'O. glumaepatula (AA)',
    'BB_Opun': 'O. punctata (BB)',
    'CC_Ooff': 'O. officinalis (CC)',
    'EE_Oaus': 'O. australiensis (EE)',
    'FF_Obra': 'O. brachyantha (FF)',
    'GG_Omey': 'O. meyeriana (GG)'
}


def load_and_process_data(input_path):
    """Load data, extract species and haplotype names, and prepare structured data for plotting."""

    try:
        # Assume input file is in TSV format, columns similar to provided examples
        # Note: Using r'\s+' as separator to handle both TAB and multiple spaces
        df = pd.read_csv(
            input_path,
            sep=r'\s+',
            header=0,
            names=['ChrID', 'Start', 'End', 'Assembly'],
            dtype={'ChrID': str, 'Assembly': str}
        )
    except Exception as e:
        print(f"Error reading input file {input_path}: {e}", file=sys.stderr)
        sys.exit(1)

    # 2. Extract Species and Haplotype
    def extract_names(assembly_name):
        parts = assembly_name.split('_')
        # Assume Species part is everything except the last '_' delimiter
        if len(parts) < 2:
            haplotype = parts[-1] if parts else ""
            species = assembly_name
        else:
            haplotype = parts[-1]
            # Compatible with Assembly names like AA_Ogla_hap1
            species = '_'.join(parts[:-1])

            # If the species part is too long, perform secondary processing to match Y_AXIS_ORDER
            # e.g., if assembly is Chr01_AA_Ogla_hap1, species would be Chr01_AA_Ogla
            # We need to transform species to AA_Ogla
            species_parts = species.split('_')
            if len(species_parts) > 2 and species_parts[-2] in [s.split('_')[0] for s in Y_AXIS_ORDER]:
                 species = '_'.join(species_parts[-2:])

        return species, haplotype

    df[['Species', 'Haplotype']] = df['Assembly'].apply(
        lambda x: pd.Series(extract_names(x))
    )

    # Filter out species not in the specified list
    df = df[df['Species'].isin(Y_AXIS_ORDER)].copy()

    # 1. Create species index mapping
    species_map = {species: i for i, species in enumerate(Y_AXIS_ORDER)}
    df['Species_Index'] = df['Species'].map(species_map)

    # 2. Sort and number chromosomes
    def get_chr_number(chr_name):
        try:
            # Compatible with Chr01, Chr02, etc.
            return int(chr_name.replace('Chr', '').replace('chr', ''))
        except ValueError:
            return 999

    df['Chr_Index'] = df['ChrID'].apply(get_chr_number)

    # 3. Sorting
    df = df.sort_values(by=['Species_Index', 'Haplotype', 'Chr_Index']).reset_index(drop=True)

    return df

def plot_faceted_centromeres(df, output_path):
    """Plot faceted figures with Position on X-axis, Species on Y-axis, faceted by Chromosome."""

    all_species = df['Species'].unique()
    all_haplotypes = df['Haplotype'].unique()

    sorted_chrs = df[['ChrID', 'Chr_Index']].drop_duplicates().sort_values('Chr_Index')
    num_chrs = len(sorted_chrs)

    if num_chrs == 0:
        print("Error: No valid chromosome data found for plotting.", file=sys.stderr)
        return

    # ----------------------------------------------------
    # ** Modification: Calculate X-axis range for each chromosome **
    # ----------------------------------------------------
    chr_x_ranges = {}
    for chr_id in df['ChrID'].unique():
        if chr_id in ['Chr04', 'Chr09']:
            # For Chr04 and Chr09, calculate range based on original logic
            chr_data = df[df['ChrID'] == chr_id]
            min_val = min(chr_data['Start'].min(), chr_data['End'].min())
            max_val = max(chr_data['Start'].max(), chr_data['End'].max())
            x_min = max(0, min_val - 5)
            x_max = max_val + 5
        else:
            # For other chromosomes, fix X-axis range to [30, 60]
            x_min = 30
            x_max = 60
        chr_x_ranges[chr_id] = (x_min, x_max)

    # --- Setup plotting canvas ---
    # Do not share X-axis (as each range varies)
    fig, axes = plt.subplots(
        1, num_chrs,
        figsize=(2.2 * num_chrs, 7),
        sharey=True,
        sharex=False, # Set to False
        squeeze=False
    )
    axes = axes[0]

    # *** Color Assignment ***
    haplotype_colors = {}
    cmap_extra = plt.cm.get_cmap('Set2')

    for i, hap in enumerate(all_haplotypes):
        if i == 0:
            haplotype_colors[hap] = '#dfaf69'
        elif i == 1:
            haplotype_colors[hap] = '#8bb1a1'
        else:
            haplotype_colors[hap] = cmap_extra((i - 2) % cmap_extra.N)

    # Get existing species and their indices in Y_AXIS_ORDER
    species_to_plot = [s for s in Y_AXIS_ORDER if s in all_species]
    species_map = {species: i for i, species in enumerate(Y_AXIS_ORDER) if species in all_species}

    # Prepare Y-axis labels
    y_tick_indices = [species_map[s] for s in species_to_plot]
    y_tick_labels = [Y_AXIS_LABELS[s] for s in species_to_plot]


    print(f"Plotting {num_chrs} chromosomes for {len(species_to_plot)} species.")

    # --- Iterate through each chromosome to plot facets ---
    for k, (chr_id, chr_index) in enumerate(sorted_chrs.values):
        ax = axes[k]
        chr_data = df[df['ChrID'] == chr_id]

        # 1. Plot data (remains unchanged)
        for hap in all_haplotypes:
            hap_data = chr_data[chr_data['Haplotype'] == hap]
            if hap_data.empty:
                continue

            color = haplotype_colors[hap]
            hap_data_sorted = hap_data.sort_values(by='Species_Index')

            # (1) Plot Start points and lines
            ax.plot(hap_data_sorted['Start'], hap_data_sorted['Species_Index'],
                    marker='o', linestyle='-', color=color, linewidth=1.0, markersize=4, zorder=3,
                    alpha=0.8,
                    label=f'{hap} Start' if k == 0 else "")

            # (2) Plot End points and lines
            ax.plot(hap_data_sorted['End'], hap_data_sorted['Species_Index'],
                    marker='o', linestyle='--', color=color, linewidth=1.0, markersize=4, zorder=3,
                    alpha=0.8,
                    label=f'{hap} End' if k == 0 else "")


        # 2. Set facet titles and axes
        ax.set_title(chr_id, fontsize=12, fontweight='bold')

        # ** Apply dynamic X-axis range **
        x_min, x_max = chr_x_ranges[chr_id]
        ax.set_xlim(x_min, x_max)

        # Set X-axis ticks to avoid overcrowding
        if chr_id not in ['Chr04', 'Chr09']:
            ax.set_xticks([20, 40, 60]) # Fixed ticks
        else:
            ax.xaxis.set_major_locator(plt.MaxNLocator(4)) # Dynamic ticks

        if k == 0:
            ax.set_ylabel("Species (Genome Group)", fontsize=12)
            ax.set_yticks(y_tick_indices)
            ax.set_yticklabels(y_tick_labels, rotation=0, fontsize=10)
            ax.invert_yaxis()

        else:
            # Remove Y-axis ticks for non-first plots, but keep shared tick lines (due to sharey=True)
            ax.tick_params(axis='y', length=0)

        if k == num_chrs // 2:
            ax.set_xlabel("Normalized Position (0 - 100)", fontsize=12)
        else:
            # Ensure X-axis label only displays on the central plot
            ax.set_xlabel("")

        ax.tick_params(axis='x', direction='inout', length=4)
        ax.grid(axis='x', linestyle=':', alpha=0.4, color='gray')
        ax.grid(axis='y', linestyle='-', alpha=0.2, color='gray')

    # --- Unified settings and legend ---

    legend_handles = []

    for hap in all_haplotypes:
        color = haplotype_colors[hap]
        line = plt.Line2D([0], [0], marker='o', color='w', label=hap,
                          markerfacecolor=color, markersize=10, linestyle='')
        legend_handles.append(line)

    start_line = plt.Line2D([0], [0], color='black', linestyle='-', linewidth=1.5, marker='o', markersize=5, label='Start positions')
    end_line = plt.Line2D([0], [0], color='black', linestyle='--', linewidth=1.5, marker='o', markersize=5, label='End positions')

    handles = [start_line, end_line] + legend_handles

    # Ensure variable names are correct
    labels = ['Start positions', 'End positions'] + list(hapothesis for hapothesis in haplotype_colors.keys())

    fig.legend(handles, labels,
               loc='center right',
               title="Haplotypes",
               bbox_to_anchor=(1.0, 0.5),
               frameon=True)

    fig.suptitle("Normalized Centromere Position Across Species (Faceted by Chromosome)", fontsize=16, y=1.02, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 0.9, 1])

    # Save the plot
    try:
        plt.savefig(output_path, dpi=300)
        print(f"Plot successfully saved to {output_path}")
    except Exception as e:
        print(f"Error saving plot to {output_path}: {e}", file=sys.stderr)


def main():
    args = parse_arguments()

    if not os.path.exists(args.input):
        print(f"Error: Input file not found at {args.input}", file=sys.stderr)
        sys.exit(1)

    df = load_and_process_data(args.input)

    if df.empty:
        print("Warning: Dataframe is empty after filtering/processing. Check input data and species names.", file=sys.stderr)
        return

    plot_faceted_centromeres(df, args.output)

if __name__ == "__main__":
    main()