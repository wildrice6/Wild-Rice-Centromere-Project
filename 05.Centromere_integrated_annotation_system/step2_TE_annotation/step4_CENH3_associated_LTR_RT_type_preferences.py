#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
from matplotlib.collections import PatchCollection
import re
import warnings

# --- Global Plotting Settings ---
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['font.size'] = 8
DPI_VALUE = 300

# Suppress warnings to maintain clean output
warnings.filterwarnings("ignore")


def parse_arguments():
    parser = argparse.ArgumentParser(description="Visualize TE delta_signal with Error Bands and Insertion Time Tracks.")

    default_input = 'merged_TE_CENH3_LTR_identity.tsv'
    default_config = 'plot_species_col.tsv'
    default_output = '06.Final_Plots_With_Tracks'

    parser.add_argument('-i', '--input_file', default=default_input)
    parser.add_argument('-o', '--output_dir', default=default_output)
    parser.add_argument('--min_points', type=int, default=50)
    parser.add_argument('--config', default=default_config)
    parser.add_argument('--width', type=float, default=15.0)
    parser.add_argument('--height', type=float, default=10.0)
    parser.add_argument('--select_tes', default='All')
    parser.add_argument('--show_context_tracks', default='yes')
    # Although lines are no longer drawn, this parameter is retained to control the normalization range of the color map
    parser.add_argument('--insertion_time_range', default='min-max', 
                       help="Insertion time range for normalization: '0-1' or 'min-max' (default)")
    return parser.parse_args()


def preprocess_data(df: pd.DataFrame) -> pd.DataFrame:
    print("--- Preprocessing data ---")
    try:
        df['TE_ID'] = df['TE_ID'].astype(str)
        id_parts = df['TE_ID'].str.split('.', expand=True)

        if id_parts.shape[1] < 3:
            df['Sample_ID'] = df['Species'] + '_' + df['Haplotype']
            df['TE_Type_Detailed'] = df['Type'].str.replace('LTR_', '')
            df['Chr_ID'] = df['Chr']
            df['Species_ID'] = df['Species']
        else:
            df['Sample_ID'] = id_parts[0]
            df['TE_Type_Detailed'] = id_parts[1]
            df['Chr_ID'] = id_parts[2]
            df['Species_ID'] = df['Sample_ID'].str.rsplit('_', n=1).str[0]

        top_n_types = 15
        type_counts = df['TE_Type_Detailed'].value_counts()
        major_types = type_counts.nlargest(top_n_types).index
        df['TE_Type_Grouped'] = df['TE_Type_Detailed'].apply(lambda x: x if x in major_types else 'Other')
        
        if 'Insertion_Time' not in df.columns:
            print("Warning: 'Insertion_Time' column not found.")
            df['Insertion_Time'] = np.nan
        
        return df
    except Exception as e:
        print(f"Error parsing TE info: {e}", file=sys.stderr)
        sys.exit(1)


def analysis_1_boxplot_with_stripplot(df, output_dir, height, selected_tes):
    print("\n--- Plotting Analysis 1: Boxplot ---")
    df_plot = df[df['TE_Type_Grouped'] != 'Other'].copy()
    if selected_tes: df_plot = df_plot[df_plot['TE_Type_Grouped'].isin(selected_tes)]

    if df_plot.empty:
        print("Warning: No data for boxplot.")
        return

    order = df_plot.groupby('TE_Type_Grouped')['delta_signal'].median().sort_values().index
    dynamic_width = max(8, len(order) * 0.5)

    plt.figure(figsize=(dynamic_width, height))

    sns.boxplot(data=df_plot, x='TE_Type_Grouped', y='delta_signal', order=order,
                     palette='coolwarm', hue='TE_Type_Grouped', legend=False,
                     showfliers=False, width=0.6)

    sns.stripplot(data=df_plot, x='TE_Type_Grouped', y='delta_signal', order=order,
                  s=2, color=".3", alpha=0.2, jitter=0.2)

    plt.axhline(0, color='grey', linestyle='--', linewidth=1)
    plt.xticks(rotation=45, ha='right')

    try:
        plt.tight_layout()
    except UserWarning:
        pass

    plt.savefig(os.path.join(output_dir, "1_TE_Type_boxplot.pdf"), dpi=DPI_VALUE)
    plt.close()


def draw_barcode_track(ax, data, color_map, track_type='categorical', norm=None):
    patches = []
    colors = []
    data_list = data.tolist()

    for i, value in enumerate(data_list):
        rect = mpatches.Rectangle((i, 0), 1, 1)
        patches.append(rect)
        if track_type == 'categorical':
            color = color_map.get(value, (0.8, 0.8, 0.8, 1))
        elif track_type == 'continuous':
            # Color mapping is applied only if the value is valid and norm exists; otherwise, use light grey
            if pd.notnull(value) and norm is not None:
                color = color_map(norm(value))
            else:
                color = (0.9, 0.9, 0.9, 1)
        else:  # Binary
            color = (0.2, 0.2, 0.2, 1.0) if value == 1 else (1, 1, 1, 0)
        colors.append(color)

    collection = PatchCollection(patches, facecolors=colors, edgecolors='none')
    ax.add_collection(collection)
    ax.set_xlim(0, len(data_list))
    ax.set_ylim(0, 1)


def analysis_2_ranked_signal_with_errorbars(df, output_dir, color_config, width, height, selected_tes, show_context, insertion_time_range='min-max'):
    print("\n--- Plotting Analysis 2: Ranked Signal with Heatmap Tracks ---")

    # Sort by delta_signal
    df_sorted = df.sort_values('delta_signal', ascending=False).reset_index(drop=True)

    x_vals = df_sorted.index
    y_vals = df_sorted['delta_signal']
    y_std = df_sorted['delta_std'].fillna(0)
    
    # Determine whether to render the Insertion Time Track
    has_insertion_time = 'Insertion_Time' in df_sorted.columns and not df_sorted['Insertion_Time'].isna().all()

    all_major = [t for t in df_sorted['TE_Type_Grouped'].unique() if t != 'Other']
    te_types_plot = [t for t in selected_tes if t in all_major] if selected_tes else all_major

    # Calculate track count: TE types + (context including LTR, Chr, Sp if shown) + (insertion time if present)
    num_tracks = len(te_types_plot)
    if show_context: num_tracks += 3
    if has_insertion_time: num_tracks += 1
    
    height_ratios = [8] + [1] * num_tracks

    fig, axes = plt.subplots(nrows=1 + num_tracks, ncols=1, sharex=True, figsize=(width, height),
                             gridspec_kw={'height_ratios': height_ratios, 'hspace': 0.1})
    if num_tracks == 0: axes = [axes]

    ax_main = axes[0]
    
    # Render primary plot: delta signal and error bands
    ax_main.fill_between(x_vals, y_vals - y_std, y_vals + y_std,
                         color='grey', alpha=0.3, label='Standard Deviation', rasterized=True)
    ax_main.scatter(x_vals, y_vals, s=1, alpha=0.6, color='steelblue', 
                   label='Mean Delta Signal', rasterized=True)
    
    ax_main.legend(loc='upper right', fontsize=7)
    ax_main.set_title('TEs Ranked by CENH3 Delta Signal (with Error Range and Attribute Tracks)')
    ax_main.set_ylabel('Delta Signal (±SD)')
    ax_main.axhline(0, color='red', linestyle='--', linewidth=0.8)
    ax_main.grid(True, alpha=0.3, linestyle='--')

    curr_ax = 1

    # 1. Render binary tracks for TE types
    for te in te_types_plot:
        ax = axes[curr_ax]
        draw_barcode_track(ax, (df_sorted['TE_Type_Grouped'] == te).astype(int), None, track_type='binary')
        ax.set_ylabel(te, rotation=0, ha='right', va='center', fontsize=6)
        ax.set_yticks([])
        curr_ax += 1

    # 2. Render Insertion Time Track (Heatmap style)
    if has_insertion_time:
        ax = axes[curr_ax]
        # Utilize gradient colormaps such as plasma or viridis
        it_cmap = plt.get_cmap('plasma') 
        it_vals = df_sorted['Insertion_Time']
        
        # Normalization configuration
        it_min, it_max = it_vals.min(), it_vals.max()
        if insertion_time_range == '0-1':
            it_norm_obj = mcolors.Normalize(vmin=0, vmax=1)
        else:
            it_norm_obj = mcolors.Normalize(vmin=it_min, vmax=it_max)
            
        draw_barcode_track(ax, it_vals, it_cmap, track_type='continuous', norm=it_norm_obj)
        ax.set_ylabel('Ins. Time', rotation=0, ha='right', va='center', fontsize=7)
        ax.set_yticks([])
        curr_ax += 1

    # 3. Render context tracks (LTR identity, Chromosome, Species)
    if show_context:
        # LTR Identity Track
        ax = axes[curr_ax]
        cmap_ltr = plt.get_cmap('magma')
        norm_ltr = mcolors.Normalize(vmin=0.8, vmax=1.0)
        draw_barcode_track(ax, df_sorted['ltr_identity'], cmap_ltr, track_type='continuous', norm=norm_ltr)
        ax.set_ylabel('LTR Id.', rotation=0, ha='right', va='center', fontsize=7)
        ax.set_yticks([])
        curr_ax += 1

        # Chromosome Track
        def natural_key(s):
            return int(re.sub(r'\D', '', str(s))) if re.sub(r'\D', '', str(s)) else 0

        uniq_chrs = sorted(df_sorted['Chr_ID'].unique(), key=natural_key)
        cmap_chr = plt.get_cmap('tab20', len(uniq_chrs))
        map_chr = {name: cmap_chr(i) for i, name in enumerate(uniq_chrs)}

        ax = axes[curr_ax]
        draw_barcode_track(ax, df_sorted['Chr_ID'], map_chr, track_type='categorical')
        ax.set_ylabel('Chr', rotation=0, ha='right', va='center', fontsize=7)
        ax.set_yticks([])
        curr_ax += 1

        # Species Track
        ax = axes[curr_ax]
        if color_config:
            map_sp = {}
            for k, v in color_config.items():
                try:
                    map_sp[k] = mcolors.to_rgba(v['color'])
                except:
                    pass
            draw_barcode_track(ax, df_sorted['Species_ID'], map_sp, track_type='categorical')
        else:
            uniq_s = sorted(df_sorted['Sample_ID'].unique())
            map_s = {n: plt.cm.viridis(i / len(uniq_s)) for i, n in enumerate(uniq_s)}
            draw_barcode_track(ax, df_sorted['Sample_ID'], map_s, track_type='categorical')
        ax.set_ylabel('Species', rotation=0, ha='right', va='center', fontsize=7)
        ax.set_yticks([])

    axes[-1].set_xlabel('TEs Ranked by Delta Signal')

    try:
        plt.tight_layout()
    except UserWarning:
        pass

    plt.savefig(os.path.join(output_dir, "2_ranked_signal_with_heatmap_tracks.pdf"), dpi=DPI_VALUE)
    plt.close()


def main():
    # Ensure execution within the script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if script_dir: os.chdir(script_dir)
    
    args = parse_arguments()
    os.makedirs(args.output_dir, exist_ok=True)

    if not os.path.exists(args.input_file):
        print(f"Error: {args.input_file} not found!")
        sys.exit(1)

    color_config = None
    if args.config and os.path.exists(args.config):
        try:
            cdf = pd.read_csv(args.config, sep='\t')
            color_config = {r['Species']: {'color': r['Color'], 'name': r['Species_name']} for i, r in cdf.iterrows()}
        except:
            pass

    df = pd.read_csv(args.input_file, sep='\t')

    # Filter out samples belonging to XX_Lhex
    initial_count = len(df)
    df = df[~df['Species'].astype(str).str.contains("Lhex", case=False, na=False)]
    df = df[~df['TE_ID'].astype(str).str.contains("Lhex", case=False, na=False)]
    print(f"--- Filtering: Removed {initial_count - len(df)} records belonging to XX_Lhex ---")

    if 'delta_std' not in df.columns:
        print("Warning: 'delta_std' column missing. Error bands will be 0.")
        df['delta_std'] = 0

    df = preprocess_data(df)

    # Filter TE types based on abundance
    cnt = df['TE_Type_Grouped'].value_counts()
    keep = cnt[(cnt >= args.min_points) & (cnt.index != 'Other')].index
    df = df[df['TE_Type_Grouped'].isin(keep) | (df['TE_Type_Grouped'] == 'Other')].copy()

    sel = [x.strip() for x in args.select_tes.split(';')] if args.select_tes.lower() != 'all' else None

    analysis_1_boxplot_with_stripplot(df, args.output_dir, args.height, sel)
    analysis_2_ranked_signal_with_errorbars(df, args.output_dir, color_config, args.width, args.height, sel,
                                            args.show_context_tracks == 'yes', args.insertion_time_range)

    print("\nDone! Check output folder.")


if __name__ == "__main__":
    main()