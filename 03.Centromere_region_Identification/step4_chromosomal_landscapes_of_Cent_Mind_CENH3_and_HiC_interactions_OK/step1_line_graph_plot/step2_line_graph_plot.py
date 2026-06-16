import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os
import sys

def load_tsv(file_path):
    """Read the signal TSV file and assign column names."""
    cols = ['chr', 'type', 'pos', 'val']
    return pd.read_csv(file_path, sep='\t', names=cols)

def load_cen_file(file_path):
    """Read the centromere position file and return a dictionary: {chr: (start, end)}."""
    cen_dict = {}
    try:
        with open(file_path, 'r') as f:
            for line in f:
                if not line.strip():
                    continue
                parts = line.split()
                if len(parts) >= 3:
                    chrom = parts[0]
                    start = float(parts[1])
                    end = float(parts[2])
                    cen_dict[chrom] = (start, end)
    except Exception as e:
        print(f"Error reading centromere file: {e}")
        sys.exit(1)
    return cen_dict

def plot_chromosome(chrom, df_hic, df_cenh3, df_prob, cen_range, output_dir):
    """Generate a triple-axis plot for a single chromosome."""
    # Filter data for the current chromosome
    d_hic = df_hic[df_hic['chr'] == chrom].sort_values('pos')
    d_cenh3 = df_cenh3[df_cenh3['chr'] == chrom].sort_values('pos')
    d_prob = df_prob[df_prob['chr'] == chrom].sort_values('pos')

    if d_cenh3.empty and d_hic.empty:
        print(f"Warning: No signal data found for chromosome {chrom}. Skipping.")
        return

    # Initialize the figure and the first axis (Left: CENH3)
    fig, ax1 = plt.subplots(figsize=(12, 7), dpi=300)

    # 1. Plot CENH3 Signal (Left axis - Red)
    color_red = '#d62728' 
    ax1.plot(d_cenh3['pos'], d_cenh3['val'], color=color_red, linewidth=1.5, label='CENH3')
    ax1.set_ylabel('CENH3 Signal', color=color_red, fontweight='bold', fontsize=12)
    ax1.tick_params(axis='y', labelcolor=color_red)
    
    # 2. Plot Hi-C Contact (Right axis 1 - Blue)
    ax2 = ax1.twinx()
    color_blue = '#1f77b4'
    ax2.plot(d_hic['pos'], d_hic['val'], color=color_blue, linewidth=1.5, label='Hi-C')
    ax2.set_ylabel('Hi-C Contact', color=color_blue, fontweight='bold', fontsize=12)
    ax2.tick_params(axis='y', labelcolor=color_blue)
    
    # 3. Plot Probability (Right axis 2 - Yellow)
    ax3 = ax1.twinx()
    ax3.spines['right'].set_position(('outward', 65))
    # Use a distinct dark golden yellow for better visibility on white backgrounds
    color_yellow = '#ED8828' 
    ax3.plot(d_prob['pos'], d_prob['val'], color=color_yellow, linewidth=1.5, label='Probability')
    ax3.set_ylabel('Probability', color=color_yellow, fontweight='bold', fontsize=12)
    ax3.tick_params(axis='y', labelcolor=color_yellow)

    # Handle centromere region visualization logic
    if cen_range:
        cen_min, cen_max = cen_range
        cen_mid = (cen_min + cen_max) / 2
        
        # Plot grey shaded area for the centromere
        ax1.axvspan(cen_min, cen_max, color='gray', alpha=0.2, zorder=0)
        # Plot center dashed line
        ax1.axvline(x=cen_mid, color='gray', linestyle='--', linewidth=1.2, alpha=0.6)
        
        # Set X-axis ticks: Map 'CEN' to the centromere midpoint
        ax1.set_xticks([-1, cen_mid, 1])
    else:
        # Default ticks if no centromere information is provided
        ax1.set_xticks([-1, 0, 1])
        print(f"Warning: Chromosome {chrom} position not defined in the --cen file.")

    ax1.set_xticklabels(['TEL', 'CEN', 'TEL'], fontsize=11)
    ax1.set_xlabel('Normalized Position', fontsize=12)
    ax1.set_xlim([-1.05, 1.05])

    # Set title
    plt.title(f'Chromosome: {chrom}', fontsize=16, fontweight='bold', pad=20)
    
    # Adjust layout
    fig.tight_layout()

    # Save outputs
    out_path = os.path.join(output_dir, f"{chrom}_signal_profile")
    plt.savefig(f"{out_path}.png", bbox_inches='tight')
    plt.savefig(f"{out_path}.pdf", bbox_inches='tight')
    plt.close(fig)

def main():
    parser = argparse.ArgumentParser(description="Generate triple-axis chromosomal signal plots.")
    parser.add_argument('--hic', required=True, help="Hi-C signal TSV file")
    parser.add_argument('--cenh3', required=True, help="CENH3 signal TSV file")
    parser.add_argument('--prob', required=True, help="Probability signal TSV file")
    parser.add_argument('--cen', required=True, help="Centromere normalized position file")
    parser.add_argument('--output_dir', required=True, help="Output directory")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    print("Loading data...")
    df_hic = load_tsv(args.hic)
    df_cenh3 = load_tsv(args.cenh3)
    df_prob = load_tsv(args.prob)
    cen_info = load_cen_file(args.cen)

    # Identify all chromosomes to process (using CENH3 file as reference)
    chroms = sorted(df_cenh3['chr'].unique())
    
    print(f"Detected {len(chroms)} chromosomes. Starting visualization...")
    for chrom in chroms:
        print(f"Processing: {chrom}")
        # Retrieve the centromere range for the current chromosome
        current_cen_range = cen_info.get(chrom, None)
        plot_chromosome(chrom, df_hic, df_cenh3, df_prob, current_cen_range, args.output_dir)
    
    print(f"Task completed! Images saved to: {args.output_dir}")

if __name__ == "__main__":
    main()