#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import re

def get_region_sort_key(region_name):
    """
    Generate a sortable key for region names.
    For example, 'p_arm_10' -> (0, 10), 'centromere_5' -> (1, 5), 'q_arm_1' -> (2, 1).
    This ensures that regions are ordered as p_arm -> centromere -> q_arm, 
    with numeric parts sorted by their numerical value.
    """
    # Define the intended order of region types
    prefix_order = {
        'L_arm': 0, 'p_arm': 0,
        'centromere': 1,
        'R_arm': 2, 'q_arm': 2
    }
    
    match = re.match(r'([a-zA-Z_]+)_(\d+)', region_name)
    if match:
        prefix = match.group(1)
        number = int(match.group(2))
        order = prefix_order.get(prefix, 99)
        return (order, number)
    
    return (99, 0)


def main():
    """
    Main function: Read data, sort, plot, and save results.
    """
    # 1. Set up command-line argument parsing
    parser = argparse.ArgumentParser(
        description="Plot LTR proportions over ordered genomic regions.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-i', '--input', required=True, help='Path to the input TSV file.')
    parser.add_argument('-o', '--output_base', required=True, help='Base name for the output plot files.')
    args = parser.parse_args()

    # 2. Read and prepare data
    try:
        print(f"--- Reading data from: {args.input} ---")
        df = pd.read_csv(args.input, sep='\t')
        df.rename(columns={df.columns[0]: 'region'}, inplace=True)
        print("Data loaded successfully. Preview:")
        print(df.head())
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return

    # 3. Perform natural sorting on data
    print("\n--- Sorting data by region ---")
    df['sort_key'] = df['region'].apply(get_region_sort_key)
    df.sort_values('sort_key', inplace=True)
    df.drop(columns='sort_key', inplace=True)
    print("Data sorted. New order preview:")
    print(df.head())

    # 4. Generate plots
    print("\n--- Generating plot ---")
    
    fig, ax = plt.subplots(figsize=(22, 8))
    df.set_index('region', inplace=True)
    
    # --- Modification 1: Add background highlighting for centromere regions ---
    # Retrieve the list of all region names to find indices
    region_list = df.index.tolist()
    # Identify numerical indices of all regions starting with 'centromere'
    centro_indices = [i for i, region in enumerate(region_list) if region.startswith('centromere')]
    
    if centro_indices:
        # Determine the start and end positions for highlighting
        # -0.5 and +0.5 are used to ensure the rectangle fully covers the data point positions
        start_idx = centro_indices[0] - 0.5
        end_idx = centro_indices[-1] + 0.5
        # Use axvspan to draw the vertical rectangular area
        ax.axvspan(start_idx, end_idx, color='lightgray', alpha=0.5, zorder=0)

    # ==================== New/Modified Section Start ====================
    # --- New: Define mapping dictionary from LTR types to colors ---
    color_map = {
        'LTR_Gypsy_Retand': '#f7ac1a',
        'LTR_Gypsy_Tekay': '#b1d238',
        'LTR_Gypsy_Ogre': '#4198b9',
        'LTR_Copia_Angela': '#afdab5',
        'LTR_Gypsy_CRM': '#ed1c24',
        'LTR_Copia_Ale': '#65082f',
        'LTR_Copia_Ivana': '#e3c7d8',
        'LTR_Gypsy_unknown': '#465867',
        'LTR_Gypsy_Reina': '#cde8f3',
        'LTR_Copia_SIRE': '#6639a6',
        'LTR_Copia_Ikeros': '#e67b79',
        'LTR_Copia_TAR': '#d1a982',
        'LTR_Copia_Bianca': '#99bcdd',
        'LTR_Copia_unknown': '#c7c7c9',
        'LTR_Copia_Tork': '#c1c758',
        'LTR_Gypsy_Athila': '#ffc808',
        'others': '#000000'
    }

    # --- Modification: Iterate through each column and plot line charts with specified colors ---
    for ltr_type in df.columns:
        # Retrieve color from dictionary; default to 'others' color if key is not found
        color = color_map.get(ltr_type, color_map['others'])
        ax.plot(df.index, df[ltr_type], 
                label=ltr_type, 
                color=color,  # Apply assigned color
                linewidth=1.5, 
                marker='o', 
                markersize=4)
    # ==================== New/Modified Section End ====================

    # 5. Beautify the plot
    ax.set_title('LTR Proportions Across Genomic Regions', fontsize=18)
    ax.set_ylabel('Proportion', fontsize=14)
    ax.grid(axis='y', linestyle='--', alpha=0.7)
    ax.legend(title='LTR Types', bbox_to_anchor=(1.02, 1), loc='upper left')

    # --- Modification 2: Remove original X-axis ticks and add custom region annotations ---
    
    # First, remove original x-axis ticks and labels
    ax.set_xticks([])
    ax.set_xlabel('')

    # Define region groups for drawing labels below the axis
    region_groups = {
        'Left Arm': ['L_arm', 'p_arm'],
        'Centromere': ['centromere'],
        'Right Arm': ['R_arm', 'q_arm']
    }

    # Set vertical positions for labels (below the x-axis)
    # Using coordinate transformation: y-coordinate represents the proportion of the axis (0=bottom, 1=top)
    # Negative values indicate positions below the axis domain
    y_line_pos = -0.05
    y_text_pos = -0.07
    transform = ax.get_xaxis_transform() # Transformation for x-data coordinates and y-axis coordinates

    for label, prefixes in region_groups.items():
        # Find indices for all regions belonging to the current group
        indices = [i for i, region in enumerate(region_list) if any(region.startswith(p) for p in prefixes)]
        
        if not indices:
            continue # Skip if this region type is not present in the data

        start_idx = indices[0]
        end_idx = indices[-1]
        mid_idx = (start_idx + end_idx) / 2

        # Draw the horizontal line
        ax.plot([start_idx, end_idx], [y_line_pos, y_line_pos],
                color='black', linewidth=1.5,
                transform=transform, clip_on=False) # clip_on=False allows drawing outside the axis area
        
        # Add the text label
        ax.text(mid_idx, y_text_pos, label,
                ha='center', va='top', fontsize=12,
                transform=transform, clip_on=False)

    # Adjust overall layout to leave space for the legend and annotations
    # Use plt.subplots_adjust to increase the bottom margin
    plt.subplots_adjust(bottom=0.15, right=0.85)


    # 6. Save plot files
    output_png = f"{args.output_base}.png"
    output_pdf = f"{args.output_base}.pdf"

    try:
        print(f"\n--- Saving plots to {output_png} and {output_pdf} ---")
        # Use bbox_inches='tight' to ensure all elements (legend, labels) are fully preserved
        plt.savefig(output_png, dpi=300, bbox_inches='tight')
        plt.savefig(output_pdf, bbox_inches='tight')
        print("Plots saved successfully!")
    except Exception as e:
        print(f"An error occurred while saving the plots: {e}")


if __name__ == "__main__":
    main()