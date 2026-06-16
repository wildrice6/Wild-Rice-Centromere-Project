#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import re
import os

def read_genome_sizes(filepath):
    """
    Reads a chromosome size file.
    Expected file format: Chr_Name<tab>Chr_Size
    """
    sizes = {}
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if line.strip():
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        chrom_name, size = parts[0], parts[1]
                        try:
                            sizes[chrom_name] = int(size)
                        except ValueError:
                            print(f"Warning: Unable to parse size in line '{line.strip()}', skipping.")
    except FileNotFoundError:
        print(f"Error: Genome size file '{filepath}' not found.")
        return None
    return sizes

def read_blast_output(filepath):
    """
    Reads a BLAST output file in outfmt 6 format.
    """
    col_names = [
        'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
    ]
    try:
        df = pd.read_csv(filepath, sep='\t', header=None, names=col_names, comment='#')
        return df
    except FileNotFoundError:
        print(f"Error: BLAST output file '{filepath}' not found.")
        return None
    except pd.errors.EmptyDataError:
        print(f"Warning: BLAST output file '{filepath}' is empty.")
        return pd.DataFrame(columns=col_names)

def read_centromere_regions(filepath):
    """
    Reads a centromere regions file.
    Expected file format: Chr_Name<tab>Start<tab>End
    
    Returns:
        dict: A dictionary where keys are chromosome names and values are (start, end) tuples.
    """
    regions = {}
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if line.strip():
                    parts = line.strip().split()
                    if len(parts) >= 3:
                        chrom, start, end = parts[0], parts[1], parts[2]
                        try:
                            regions[chrom] = (int(start), int(end))
                        except ValueError:
                            print(f"Warning: Unable to parse centromere coordinates in line '{line.strip()}', skipping.")
    except FileNotFoundError:
        print(f"Error: Centromere regions file '{filepath}' not found.")
        return None
    return regions


def natural_sort_key(s):
    """
    Key function for natural sorting (e.g., 'Chr1', 'Chr2', 'Chr10').
    """
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', s)]

def plot_hits_on_chromosomes(blast_df, genome_sizes, output_file, centromere_regions=None):
    """
    Plots the distribution of BLAST hits on chromosomes with optional centromere annotations.
    """
    if blast_df is None or genome_sizes is None:
        print("Input data error; unable to generate plot.")
        return
        
    if centromere_regions is None:
        centromere_regions = {}

    valid_chroms = [chrom for chrom in blast_df['sseqid'].unique() if chrom in genome_sizes]
    blast_df = blast_df[blast_df['sseqid'].isin(valid_chroms)].copy()

    if blast_df.empty:
        print("Warning: No alignment results found on the provided chromosomes. Proceeding with empty chromosome map.")
    
    blast_df['hit_center'] = (blast_df['sstart'] + blast_df['send']) / 2
    sorted_chroms = sorted(genome_sizes.keys(), key=natural_sort_key)
    
    fig_height = max(6, len(sorted_chroms) * 0.5)
    fig, ax = plt.subplots(figsize=(12, fig_height))

    y_positions = range(len(sorted_chroms))
    
    for i, chrom in enumerate(sorted_chroms):
        chrom_length = genome_sizes[chrom]
        
        # 1. Draw chromosome background (grey horizontal bars)
        ax.barh(y=i, width=chrom_length, height=0.5, color='lightgrey', edgecolor='black', linewidth=0.5)
        
        # 2. Draw centromere regions (red horizontal bars)
        if chrom in centromere_regions:
            cen_start, cen_end = centromere_regions[chrom]
            cen_width = cen_end - cen_start
            ax.barh(y=i, width=cen_width, left=cen_start, height=0.6, color='#D62728', edgecolor='darkred', linewidth=0.5, zorder=2)

        # 3. Filter and draw hits for the current chromosome (blue vertical short lines)
        chrom_hits = blast_df[blast_df['sseqid'] == chrom]
        if not chrom_hits.empty:
            hit_positions = chrom_hits['hit_center'].values
            ax.vlines(x=hit_positions, ymin=i-0.4, ymax=i+0.4, color='#1F77B4', alpha=0.7, linewidth=0.1, zorder=3)

    # Finalize plot aesthetics
    ax.set_yticks(y_positions)
    ax.set_yticklabels(sorted_chroms)
    ax.invert_yaxis()

    # Format X-axis to display in Megabases (Mb)
    formatter = mticker.FuncFormatter(lambda x, pos: f'{x/1e6:.0f}')
    ax.xaxis.set_major_formatter(formatter)
    
    ax.set_xlabel('Position (Mb)')
    ax.set_ylabel('Chromosome')
    ax.set_title('Distribution of BLAST Hits on Chromosomes')
    
    # Add legend
    legend_elements = [
        mpatches.Patch(facecolor='lightgrey', edgecolor='black', label='Chromosome'),
        mpatches.Patch(facecolor='#D62728', edgecolor='darkred', label='Centromere'),
        mlines.Line2D([], [], color='#1F77B4', marker='|', linestyle='None', markersize=10, markeredgewidth=1.5, label='BLAST Hit')
    ]
    ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1, 1.05))

    ax.grid(axis='x', linestyle='--', alpha=0.6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Adjust layout to accommodate the title and legend
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    try:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot successfully saved to: {output_file}")
    except Exception as e:
        print(f"Error: Unable to save plot. Reason: {e}")
    
    plt.close()


def main():
    """
    Main function to parse arguments and execute the plotting pipeline.
    """
    parser = argparse.ArgumentParser(
        description="Visualize the distribution of hits on chromosomes based on BLAST results, genome sizes, and centromere coordinates."
    )
    parser.add_argument('--blast_output', required=True, help='Path to the input BLAST outfmt 6 file.')
    parser.add_argument('--genome_size', required=True, help='Path to the input genome size file (two columns: Chr_Name Chr_Size).')
    parser.add_argument('--cen', required=False, help='(Optional) Path to the centromere regions file (three columns: Chr_Name Start End).')
    parser.add_argument('--output', required=True, help='Path to the output image file (e.g., plot.png, plot.svg).')

    args = parser.parse_args()

    # 1. Read genome sizes
    print("Reading genome size file...")
    genome_sizes = read_genome_sizes(args.genome_size)

    # 2. Read BLAST results
    print("Reading BLAST output file...")
    blast_hits_df = read_blast_output(args.blast_output)

    # 3. Read centromere regions (optional)
    centromere_regions = None
    if args.cen:
        print("Reading centromere regions file...")
        centromere_regions = read_centromere_regions(args.cen)

    # 4. Generate plot
    print("Generating distribution plot...")
    plot_hits_on_chromosomes(blast_hits_df, genome_sizes, args.output, centromere_regions)

    print("Processing complete!")


if __name__ == "__main__":
    main()