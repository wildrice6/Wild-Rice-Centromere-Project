#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import argparse
import sys
import numpy as np

def plot_custom_grid_final_with_circles(input_file, output_file):
    """
    Script for generating a highly customized grid plot with bars and proportional circles.
    """
    # --- 1. Configure Matplotlib font to Arial ---
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.rcParams['axes.unicode_minus'] = False

    # --- Define specified color palette ---
    COLOR_RED = '#c44e52'
    COLOR_BLUE = '#4c72b0'
    COLOR_CIRCLE = '#DDBF84'

    # --- 2. Data loading and preprocessing ---
    try:
        column_names = [
            'material_name', 'haplotype', 'chromosome', 'window_size',
            'position', 'copy_number', 'proportion'
        ]
        df = pd.read_csv(input_file, sep='\t', header=None, names=column_names)
        
        for col in ['window_size', 'copy_number', 'proportion']:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        df.dropna(subset=['window_size', 'copy_number', 'proportion'], inplace=True)
        df['chr_num'] = df['chromosome'].str.extract(r'(\d+)').astype(int)
        df = df.sort_values(by=['material_name', 'chr_num'])
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred while reading or processing the file: {e}", file=sys.stderr)
        sys.exit(1)

    # --- 3. Prepare data structures for plotting ---
    custom_order = [
        'O. sativa (japonica)', 'O. sativa (indica)', 'O. glaberrima',
        'O. rufipogon', 'O. nivara', 'O. longistaminata', 'O. glumaepatula',
        'O. punctata', 'O. officinalis', 'O. australiensis', 'O. brachyantha',
        'O. meyeriana'
    ]
    all_materials_in_data = set(df['material_name'].unique())
    materials = [name for name in custom_order if name in all_materials_in_data]
    
    chromosomes = sorted(df['chromosome'].unique(), key=lambda x: int(''.join(filter(str.isdigit, x))))
    
    chromosome_spacing = 1.8 
    material_map = {name: i for i, name in enumerate(materials)}
    chromosome_map = {name: i * chromosome_spacing for i, name in enumerate(chromosomes)}

    # --- 4. Initialize plot ---
    fig_width = max(20, len(chromosomes) * 2.2 * chromosome_spacing) 
    fig_height = max(8, len(materials) * 0.7)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # --- 5. Define scaling function and iterate through data for plotting ---
    bar_container_width = 0.8
    min_cn = df['copy_number'].min()
    max_cn = df['copy_number'].max()
    
    def scale_copy_number(cn, min_val, max_val, min_size=30, max_size=1500):
        if max_val == min_val:
            return (min_size + max_size) / 2
        # Apply square root scaling to ensure circle area is proportional to the value
        normalized = np.sqrt(cn / max_val)
        return min_size + normalized * (max_size - min_size)

    for index, row in df.iterrows():
        if row['material_name'] not in material_map or row['chromosome'] not in chromosome_map:
            continue
            
        x = chromosome_map[row['chromosome']]
        y = material_map[row['material_name']]

        bar_height = 0.6
        bar_start_x = x - bar_container_width / 2
        bar_length = (row['proportion'] / 100.0) * bar_container_width
        bar_color = COLOR_RED if row['window_size'] <= 500 else COLOR_BLUE

        rect = mpatches.Rectangle(
            (bar_start_x, y - bar_height / 2), 
            bar_length, bar_height,
            facecolor=bar_color, edgecolor=None, zorder=2
        )
        ax.add_patch(rect)

        circle_x = bar_start_x + bar_length + (0.2 * chromosome_spacing) 
        circle_y = y
        circle_size = scale_copy_number(row['copy_number'], min_cn, max_cn)

        ax.scatter(
            circle_x, circle_y, s=circle_size, 
            color=COLOR_CIRCLE, alpha=0.8, edgecolors='grey',
            linewidth=0.5, zorder=3
        )
        
    # --- 6. Configure plot aesthetics ---
    ax.set_title('Dominance of Repeats in Centromeres', fontsize=16, pad=35)
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    ax.set_xticks(np.arange(len(chromosomes)) * chromosome_spacing)
    ax.set_xticklabels(chromosomes, fontsize=12, rotation=0)
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_yticks(np.arange(len(materials)))
    ax.set_yticklabels(materials, fontsize=11)
    ax.set_ylim(len(materials) - 0.5 + 1.2, -0.5) 
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.grid(False)
    ax.tick_params(axis='both', which='both', length=0)

    # Render the bottom scale bar
    scale_y_pos = len(materials) - 0.5 + 0.3
    for i in range(len(chromosomes)):
        chromosome_x_pos = i * chromosome_spacing
        base_x = chromosome_x_pos - bar_container_width / 2
        pos_0 = base_x
        pos_40 = base_x + (40 / 100.0) * bar_container_width
        pos_80 = base_x + (80 / 100.0) * bar_container_width
        pos_60 = base_x + (60 / 100.0) * bar_container_width
        
        ax.plot([pos_60, pos_60], [-0.5, len(materials) - 0.5], color='lightgray', linestyle='--', linewidth=1, zorder=1)
        ax.plot([pos_0, pos_80], [scale_y_pos, scale_y_pos], color='black', lw=1.5)
        for px, lbl in zip([pos_0, pos_40, pos_80], ['0', '40', '80']):
            ax.plot([px, px], [scale_y_pos, scale_y_pos + 0.1], color='black', lw=1.5)
            ax.text(px, scale_y_pos - 0.2, lbl, ha='center', va='top', fontsize=10)
    
    # --- 7. Create legends ---
    red_patch = mpatches.Patch(color=COLOR_RED, label='Window Size <= 500')
    blue_patch = mpatches.Patch(color=COLOR_BLUE, label='Window Size > 500')
    color_legend = ax.legend(handles=[red_patch, blue_patch], 
                             loc='upper left',
                             bbox_to_anchor=(1.02, 1),
                             frameon=True, edgecolor='black', fontsize=10,
                             title='Window Size')
    ax.add_artist(color_legend)

    # --- Logic for generating circle size legend using round values ---
    def get_nice_legend_values(max_val):
        """ Generate round values (1, 2, 5 sequences) across orders of magnitude based on max_val """
        candidates = []
        for p in range(0, 8): # Covers 1 to 10,000,000
            for multiplier in [1, 2, 5]:
                val = multiplier * (10**p)
                if val <= max_val:
                    candidates.append(val)
        # Select the last 4 largest levels if many exist, otherwise show all
        return candidates[-4:] if len(candidates) > 4 else candidates

    legend_cns = get_nice_legend_values(max_cn)
    legend_handles = []
    for cn in legend_cns:
        handle = plt.scatter([], [], s=scale_copy_number(cn, min_cn, max_cn),
                             color=COLOR_CIRCLE, alpha=0.8, edgecolors='black',
                             label=f'{cn:,}')
        legend_handles.append(handle)

    if legend_handles:
        size_legend = ax.legend(handles=legend_handles, title='Copy Number (Circle Size)',
                                loc='upper left',
                                bbox_to_anchor=(1.02, 0.85),
                                frameon=True, edgecolor='black', fontsize=10)
        ax.add_artist(size_legend)
        
    # --- 8. Adjust layout and save the figure ---
    legend_space_inches = 3.0
    required_right_margin_ratio = legend_space_inches / fig_width
    right_param = 1.0 - required_right_margin_ratio
    fig.subplots_adjust(right=min(0.95, right_param))

    try:
        plt.savefig(output_file, dpi=300, bbox_extra_artists=(color_legend, size_legend), bbox_inches='tight', pad_inches=0.1)
        print(f"Plotting successful! Figure saved to '{output_file}'")
    except Exception as e:
        print(f"An error occurred while saving the figure: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Generates a highly customized grid plot with circles based on TSV data.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('--input', required=True, help='Path to the input TSV file.')
    parser.add_argument('--output', required=True, help='Path to the output image file.')
    args = parser.parse_args()

    plot_custom_grid_final_with_circles(args.input, args.output)