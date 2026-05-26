import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import sys

def main():
    # 1. Argument parsing setup
    parser = argparse.ArgumentParser(description="Plot scatter plot of differences between various assemblies and the reference (custom colors)")
    parser.add_argument("--input", required=True, help="Input TSV file containing difference results")
    parser.add_argument("--output", required=True, help="Path for the output image file (e.g., png/pdf/svg)")
    
    args = parser.parse_args()

    # 2. Define color and species name mapping
    # Establish dictionary based on the provided configuration
    mapping_data = {
        "AA_Ogla": {"color": "#76D273", "name": "O. glaberrima"},
        "AA_Oruf": {"color": "#215A20", "name": "O. rufipogon"},
        "AA_Oniv": {"color": "#3BA738", "name": "O. nivara"},
        "AA_Olon": {"color": "#51C54E", "name": "O. longistaminata"},
        "AA_Oglu": {"color": "#3D8347", "name": "O. glumaepatula"},
        "BB_Opun": {"color": "#F2AE2C", "name": "O. punctata"},
        "CC_Ooff": {"color": "#684E94", "name": "O. officinalis"},
        "EE_Oaus": {"color": "#4E84C3", "name": "O. australiensis"},
        "FF_Obra": {"color": "#D55F6F", "name": "O. brachyantha"},
        "GG_Omey": {"color": "#9D5427", "name": "O. meyeriana"},
        "AA_Osat": {"color": "#CBE54E", "name": "O. sativa ssp. indica"},
        "XX_Lhex": {"color": "#595959", "name": "Leersia hexandra"}
    }

    try:
        # 3. Load data
        df = pd.read_csv(args.input, sep='\t')
        
        if 'VS_assembly' not in df.columns:
            print("Error: Missing 'VS_assembly' column in the input file")
            sys.exit(1)

        # 4. Data processing: Extract species prefixes and map names/colors
        # Extract the first two parts, e.g., AA_Ogla_hap1 -> AA_Ogla
        df['Species_code'] = df['VS_assembly'].apply(lambda x: "_".join(x.split('_')[:2]))
        
        # Map scientific names (for the legend)
        df['Species_name'] = df['Species_code'].map(lambda x: mapping_data[x]['name'] if x in mapping_data else x)
        
        # Sort by chromosome
        df = df.sort_values(by=['ChrID', 'Species_name'])

        # 5. Prepare color palette
        # Prepare colors only for species present in the current dataset
        current_species = df['Species_name'].unique()
        palette = {}
        for code, info in mapping_data.items():
            palette[info['name']] = info['color']

        # 6. Plotting
        plt.figure(figsize=(12, 7))
        sns.set_style("ticks")

        # Plot y=0 as the reference baseline
        plt.axhline(y=0, color='black', linestyle='--', linewidth=1, alpha=0.5, label='Nipponbare (Baseline)')

        # Create scatter plot
        # Use Species_name for grouping and the palette for colors
        scatter = sns.scatterplot(
            data=df,
            x='ChrID',
            y='difference',
            hue='Species_name',
            style='Species_name', # Distinct markers for different species
            palette=palette,      # Apply custom colors
            s=130,                # Marker size
            alpha=0.9,
            edgecolor='w',
            linewidth=0.5
        )

        # 7. Plot aesthetics
        plt.title("Centromere Position Difference Relative to Nipponbare", fontsize=16, pad=20)
        plt.xlabel("Chromosome ID", fontsize=12)
        plt.ylabel("Normalized Difference (Offset)", fontsize=12)
        
        plt.xticks(rotation=45)
        plt.grid(axis='y', linestyle=':', alpha=0.6)

        # Adjust legend: display scientific names in italics
        leg = plt.legend(title='Species', bbox_to_anchor=(1.05, 1), loc='upper left', prop={'style': 'italic'})
        plt.setp(leg.get_title(), style='normal') # Legend title remains normal style

        plt.tight_layout()

        # 8. Save image
        plt.savefig(args.output, dpi=300)
        print(f"Plotting successful! Colors applied as requested. Results saved to: {args.output}")

    except Exception as e:
        print(f"An error occurred: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()