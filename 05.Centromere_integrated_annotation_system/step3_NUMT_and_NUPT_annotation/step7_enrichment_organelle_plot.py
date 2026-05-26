import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import numpy as np

def get_args():
    parser = argparse.ArgumentParser(description='Plot Centromere Length and Enrichment Ratio')
    parser.add_argument('--input', '-i', required=True, help='Input TSV file')
    parser.add_argument('--output', '-o', required=True, help='Output prefix (e.g., result)')
    return parser.parse_args()

def prepare_data(df):
    # Define mapping relationships for species nomenclature
    label_map = {
        "AA_Osat_hap2": "O. sativa (japonica)",
        "AA_Oruf": "O. rufipogon",
        "AA_Osat_hap1": "O. sativa (indica)",
        "AA_Oniv": "O. nivara",
        "AA_Ogla": "O. glaberrima",
        "AA_Oglu": "O. glumaepatula",  
        "AA_Olon": "O. longistaminata",
        "BB_Opun": "O. punctata",
        "CC_Ooff": "O. officinalis",
        "EE_Oaus": "O. australiensis",
        "FF_Obra": "O. brachyantha",
        "GG_Omey": "O. meyeriana"
    }

    # Define the fundamental sorting order for species (Key list)
    species_order = [
        "AA_Osat_hap2", "AA_Osat_hap1", "AA_Ogla", "AA_Oruf", "AA_Oniv", "AA_Olon",  
        "AA_Oglu", "BB_Opun", "CC_Ooff", "EE_Oaus", "FF_Obra", "GG_Omey"
    ]

    # 1. Construct sorting weights for consistency across plots
    def get_sort_weight(assembly):
        if assembly == "AA_Osat_hap2": return 0
        if assembly == "AA_Osat_hap1": return 20 # Positioned after Oruf (10)
        
        prefix = "_".join(assembly.split("_")[:2])
        if prefix in species_order:
            base_weight = species_order.index(prefix) * 10
            hap_weight = 1 if "hap1" in assembly else 2
            return base_weight + hap_weight
        return 999

    # 2. Function to convert internal identifiers to display names
    def get_display_name(assembly):
        # Direct match for full identifiers (e.g., AA_Osat_hap2)
        if assembly in label_map:
            return label_map[assembly]
        
        # Match by prefix (e.g., AA_Oruf) and append haplotype information
        prefix = "_".join(assembly.split("_")[:2])
        if prefix in label_map:
            species_name = label_map[prefix]
            hap_info = " (Hap1)" if "hap1" in assembly else " (Hap2)"
            return f"{species_name}{hap_info}"
        
        return assembly

    df['sort_weight'] = df['assembly'].apply(get_sort_weight)
    df = df.sort_values('sort_weight').reset_index(drop=True)
    
    # Apply name mapping for publication-quality labels
    df['display_name'] = df['assembly'].apply(get_display_name)
    
    # Calculate the centromeric enrichment ratio
    df['enrichment_ratio'] = df['INcen_density'] / df['OUTcen_density']
    return df

def plot_stacked_lengths(df, output_prefix):
    # Configure plotting aesthetics
    plt.rcParams['font.family'] = 'sans-serif'
    plt.figure(figsize=(16, 8))
    sns.set_style("whitegrid", {'axes.grid': True, 'grid.linestyle': '-', 'grid.alpha': 0.2})
    
    color_out = "#4c72b0" # Morandi Blue
    color_in = "#c44e52"  # Morandi Red
    
    x = np.arange(len(df))
    plt.bar(x, df['OUTcen_length'], color=color_out, label='Outside Centromere', width=0.7, edgecolor='white', linewidth=0.5)
    plt.bar(x, df['INcen_length'], bottom=df['OUTcen_length'], color=color_in, label='Inside Centromere', width=0.7, edgecolor='white', linewidth=0.5)
    
    # Configure X-axis with italicized species names via mathtext
    plt.xticks(x, [f"$\mathit{{{n}}}$" if "O." in n else n for n in df['display_name']], 
               rotation=45, ha='right', fontsize=11)
    
    plt.ylabel('Sequence Length (bp)', fontsize=13)
    plt.title('Genomic Distribution of Sequences Relative to Centromeres', fontsize=15, pad=20)
    plt.legend(frameon=True, loc='upper right', fontsize=11)
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_lengths.png", dpi=300)
    plt.savefig(f"{output_prefix}_lengths.pdf")
    plt.close()
    print(f"Successfully saved Stacked Length Plot: {output_prefix}_lengths.png/pdf")

def plot_enrichment(df, output_prefix):
    # Sort by Enrichment Ratio in descending order to highlight variation
    df_plot = df.sort_values('enrichment_ratio', ascending=False)
    
    plt.figure(figsize=(12, 10))
    sns.set_style("whitegrid")
    
    # Assign distinct colors for enriched vs. depleted states
    colors = ['#55a868' if r >= 1 else '#d65f5f' for r in df_plot['enrichment_ratio']]
    
    # Generate horizontal bar chart
    y_pos = np.arange(len(df_plot))
    bars = plt.barh(y_pos, df_plot['enrichment_ratio'], color=colors, edgecolor='none', alpha=0.9)
    
    # Set Y-axis labels with italicization
    plt.yticks(y_pos, [f"$\mathit{{{n}}}$" if "O." in n else n for n in df_plot['display_name']], fontsize=11)
    
    # Add a reference line at Enrichment Ratio = 1
    plt.axvline(x=1, color='#444444', linestyle='--', linewidth=2, label='No enrichment (ratio=1)', alpha=0.6)
    
    # Annotate specific Ratio values at the end of each bar
    for bar in bars:
        width = bar.get_width()
        plt.text(width + 0.2, bar.get_y() + bar.get_height()/2, f'{width:.2f}', 
                 va='center', fontsize=10, color='#333333')

    plt.xlabel('Enrichment Ratio (Density Inside / Density Outside)', fontsize=13)
    plt.title('Centromeric Enrichment Ratio across Species', fontsize=15, pad=20)
    plt.legend(loc='upper right')
    
    # Adjust X-axis limit to accommodate labels
    plt.xlim(0, df_plot['enrichment_ratio'].max() * 1.1)
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_enrichment.png", dpi=300)
    plt.savefig(f"{output_prefix}_enrichment.pdf")
    plt.close()
    print(f"Successfully saved Enrichment Ratio Plot: {output_prefix}_enrichment.png/pdf")

def main():
    args = get_args()
    
    try:
        # Load input TSV data
        df = pd.read_csv(args.input, sep='\t')
    except Exception as e:
        print(f"Error: Could not read file. {e}")
        return

    # Data processing and label mapping for visualization
    df = prepare_data(df)
    
    # Execute visualization functions
    plot_stacked_lengths(df, args.output)
    plot_enrichment(df, args.output)

if __name__ == "__main__":
    main()