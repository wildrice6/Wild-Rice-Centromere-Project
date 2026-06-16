--- START OF FILE 02.Chromosomal landscape of CENH3 signals.py ---

#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Function: Reads bedGraph (.bdg) files and generates a linear signal heatmap for each chromosome.
          The heatmap is presented in a bar format, where the x-axis represents chromosomal 
          positions and the color scale represents signal intensity.

Usage:
      python plot_bdg_heatmap.py --input your_data.bdg --output output_directory
"""

import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
from tqdm import tqdm

def plot_chromosome_heatmap(df_chrom, chrom_name, output_dir, vmin, vmax):
    """
    Generates and saves a heatmap for a single chromosome.

    Parameters:
    - df_chrom (pd.DataFrame): DataFrame containing data for a single chromosome.
    - chrom_name (str): Name of the chromosome.
    - output_dir (str): Directory where the output images will be saved.
    - vmin (float): Minimum value for color mapping (global).
    - vmax (float): Maximum value for color mapping (global).
    """
    # Prepare plotting data
    # xranges is a list of (start, width) tuples
    xranges = list(zip(df_chrom['start'], df_chrom['end'] - df_chrom['start']))
    signals = df_chrom['signal']

    # --- Plot Configuration ---
    fig, ax = plt.subplots(figsize=(20, 2.5))

    # Configure color mapping
    # 'coolwarm' is suitable for divergent data (blue for negative, red for positive, white for zero).
    # Sequential colormaps like 'viridis' or 'plasma' are also commonly used for continuous data.
    cmap = plt.get_cmap('coolwarm')
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    
    # Use ScalarMappable to map signal values to colors
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([]) # Required for the mappable object

    # broken_barh is ideal for plotting non-continuous horizontal bars
    # It accepts a list of (start, width) and corresponding face colors
    ax.broken_barh(xranges, yrange=(0, 1), facecolors=sm.to_rgba(signals))

    # --- Refine Axes and Titles ---
    ax.set_title(f'Signal Heatmap for {chrom_name}', fontsize=16, fontweight='bold')
    ax.set_xlabel('Chromosomal Position', fontsize=12)

    # Hide the Y-axis as it carries no physical meaning in this context
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.set_ylim(0, 1)

    # Set X-axis range
    max_pos = df_chrom['end'].max()
    ax.set_xlim(0, max_pos)

    # Format X-axis tick labels for better readability (e.g., displaying in Mb)
    def format_mb(x, pos):
        return f'{x / 1e6:.1f} Mb'
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(format_mb))
    plt.xticks(rotation=45, ha='right')

    # Add Colorbar
    cbar = fig.colorbar(sm, orientation='vertical', ax=ax, pad=0.02, aspect=10)
    cbar.set_label('Signal Value', fontsize=12)

    # Adjust layout to prevent label overlapping
    plt.tight_layout()

    # --- Save Visualization ---
    output_path = os.path.join(output_dir, f'{chrom_name}_heatmap.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig) # Close figure to release memory

def main():
    """Main function for argument parsing and executing the plotting workflow."""
    parser = argparse.ArgumentParser(description="Generate signal heatmaps for each chromosome from a .bdg file.")
    parser.add_argument('--input', '-i', required=True, help="Path to the input .bdg file.")
    parser.add_argument('--output', '-o', required=True, help="Path to the output directory to save plots.")
    args = parser.parse_args()

    # 1. Validate and create output directory
    if not os.path.exists(args.output):
        print(f"Creating output directory: {args.output}")
        os.makedirs(args.output)

    # 2. Load the bedGraph file
    try:
        print(f"Reading input file: {args.input}")
        # Using '\s+' as a separator handles both spaces and tabs
        df = pd.read_csv(
            args.input,
            sep='\s+',
            header=None,
            names=['chrom', 'start', 'end', 'signal'],
            dtype={'chrom': str, 'start': int, 'end': int, 'signal': float}
        )
    except FileNotFoundError:
        print(f"Error: Input file not found at {args.input}")
        return
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return

    # 3. Calculate global signal range to ensure consistent color scaling across all plots
    global_vmin = df['signal'].min()
    global_vmax = df['signal'].max()
    print(f"Global signal range set from {global_vmin:.2f} to {global_vmax:.2f}")

    # 4. Group by chromosome and generate plots for each group
    # Utilize groupby for efficient processing
    grouped = df.groupby('chrom')
    
    # Use tqdm to display a progress bar
    print("Generating plots for each chromosome...")
    for chrom_name, group_df in tqdm(grouped, total=len(grouped)):
        # Ensure data within each group is sorted by chromosomal position
        group_df_sorted = group_df.sort_values(by='start').reset_index(drop=True)
        plot_chromosome_heatmap(group_df_sorted, chrom_name, args.output, global_vmin, global_vmax)

    print("\nAll plots have been generated successfully!")
    print(f"Output files are saved in: {args.output}")


if __name__ == '__main__':
    main()