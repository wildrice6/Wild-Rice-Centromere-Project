import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
from collections import defaultdict

def get_multi_mapped_intervals(intervals):
    """
    Identifies nuclear genome intervals with coverage > 1 using a scan-line algorithm.
    """
    if not intervals:
        return []
    
    events = []
    for start, end in intervals:
        events.append((start, 1))
        events.append((end + 1, -1))
    
    events.sort()
    
    multi_intervals = []
    current_depth = 0
    interval_start = None
    
    for pos, change in events:
        prev_depth = current_depth
        current_depth += change
        
        # When depth changes from 1 to 2, start a multi-mapping interval
        if prev_depth == 1 and current_depth == 2:
            interval_start = pos
        # When depth changes from 2 to 1, end a multi-mapping interval
        elif prev_depth == 2 and current_depth == 1:
            if interval_start is not None:
                multi_intervals.append((interval_start, pos - 1))
                interval_start = None
                
    return multi_intervals

def plot_coverage(blast_file, size_file, output_prefix):
    # 1. Read organelle genome sizes
    print("Reading organelle genome sizes...")
    org_sizes = {}
    with open(size_file, 'r') as f:
        for line in f:
            if line.strip():
                cols = line.split('\t')
                org_sizes[cols[0]] = int(cols[1])

    # 2. First pass reading BLAST file: collect nuclear genome intervals
    print("Analyzing nuclear genome multi-mapping regions...")
    nuc_intervals = defaultdict(list)
    all_hits = [] 
    
    # Read all columns to prevent indexing errors
    try:
        df = pd.read_csv(blast_file, sep='\t', header=None)
    except Exception as e:
        print(f"Error: Unable to read input file: {e}")
        return
    
    for _, row in df.iterrows():
        org_id = str(row[0])
        chr_id = str(row[1])
        o_s, o_e = int(row[6]), int(row[7])
        n_s, n_e = int(row[8]), int(row[9])
        
        n_min, n_max = min(n_s, n_e), max(n_s, n_e)
        nuc_intervals[chr_id].append((n_min, n_max))
        all_hits.append((org_id, o_s, o_e, chr_id, n_s, n_e))

    # Calculate multi-mapping intervals for each nuclear chromosome
    nuc_multi_map = {}
    for chr_id, intervals in nuc_intervals.items():
        nuc_multi_map[chr_id] = get_multi_mapped_intervals(intervals)

    # 3. Calculate coverage arrays
    print("Calculating coverage values...")
    total_cov_dict = {sid: np.zeros(size, dtype=int) for sid, size in org_sizes.items()}
    red_cov_dict = {sid: np.zeros(size, dtype=int) for sid, size in org_sizes.items()}

    for org_id, o_s, o_e, chr_id, n_s, n_e in all_hits:
        if org_id not in total_cov_dict:
            continue
            
        os_min, os_max = min(o_s, o_e), max(o_s, o_e)
        # Update total coverage (blue section)
        total_cov_dict[org_id][os_min-1 : os_max] += 1
        
        n_min, n_max = min(n_s, n_e), max(n_s, n_e)
        # Match the multi-mapping regions on the nuclear genome for this alignment
        for m_s, m_e in nuc_multi_map[chr_id]:
            inter_s = max(n_min, m_s)
            inter_e = min(n_max, m_e)
            
            if inter_s <= inter_e:
                n_len = n_max - n_min + 1
                o_len = os_max - os_min + 1
                ratio = o_len / n_len
                
                # Map back to organelle coordinates
                if o_s < o_e: # Forward
                    rel_os = os_min + int((inter_s - n_min) * ratio)
                    rel_oe = os_min + int((inter_e - n_min) * ratio)
                else: # Reverse
                    rel_os = os_max - int((inter_e - n_min) * ratio)
                    rel_oe = os_max - int((inter_s - n_min) * ratio)
                
                rs, re = min(rel_os, rel_oe), max(rel_os, rel_oe)
                red_cov_dict[org_id][max(0, rs-1) : min(re, org_sizes[org_id])] += 1

    # 4. Plotting
    print("Rendering SVG charts...")
    # Set global font style
    plt.rcParams['font.family'] = 'sans-serif'
    # Modification: set SVG fonttype to 'none' (keeps text editable) or 'path' (converts to paths)
    # Typically set to 'none' so text remains editable when opened in AI or Inkscape
    plt.rcParams['svg.fonttype'] = 'none' 

    for org_id in org_sizes:
        if total_cov_dict[org_id].sum() == 0:
            continue
            
        fig, ax = plt.subplots(figsize=(14, 5))
        x = np.arange(1, org_sizes[org_id] + 1)
        
        # Color definitions
        color_total = 'SteelBlue'   # Represents total/unique coverage
        color_multi = 'Crimson'     # Represents multi-mapped coverage
        
        # Draw total coverage area (blue)
        ax.fill_between(x, total_cov_dict[org_id], color=color_total, 
                        label='Total/Unique Coverage', alpha=0.7, linewidth=0)
        
        # Draw multi-mapped coverage area (red, overlaid)
        ax.fill_between(x, red_cov_dict[org_id], color=color_multi, 
                        label='Multi-mapped Coverage', alpha=0.8, linewidth=0)

        # Aesthetics
        ax.set_title(f"Organelle Genome Coverage Profile: {org_id}", fontsize=16, pad=15)
        ax.set_xlabel("Position on Organelle Genome (bp)", fontsize=12)
        ax.set_ylabel("Depth of Coverage", fontsize=12)
        
        # Polishing ticks
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', linestyle='--', alpha=0.4)
        
        ax.legend(loc='upper right', frameon=True, fontsize=10)
        
        # Modification: output as SVG
        out_name = f"{output_prefix}_{org_id}.svg"
        plt.tight_layout()
        plt.savefig(out_name, format='svg')
        plt.close()
        print(f"Generated: {out_name}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plotting script: Highlight organelle coverage corresponding to nuclear multi-mapping regions")
    parser.add_argument("-i", "--input", required=True, help="Input BLASTN tsv result file")
    parser.add_argument("-s", "--size", required=True, help="Input organelle length tsv file")
    parser.add_argument("-o", "--output", default="coverage_result", help="Prefix for output files")
    
    args = parser.parse_args()
    plot_coverage(args.input, args.size, args.output)