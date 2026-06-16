import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr
import argparse
import os
import matplotlib.gridspec as gridspec

# Global Vector Graphic Settings (for Illustrator compatibility)
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.size'] = 7           
plt.rcParams['axes.titlesize'] = 8      
plt.rcParams['xtick.labelsize'] = 6
plt.rcParams['ytick.labelsize'] = 6

def run_sv_master_analysis(input_file, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 1. Load and Sort Data
    df = pd.read_csv(input_file, sep='\t')
    
    # Rename Chinese column names to English for internal processing
    df = df.rename(columns={'物种': 'Species'})
    
    target_order = [
        "O. sativa ssp. japonica b", "O. sativa ssp. japonica", "O. sativa ssp. indica", 
        "O. glaberrima", "O. rufipogon", "O. nivara", "O. longistaminata", 
        "O. glumaepatula", "O. punctata", "O. officinalis", 
        "O. australiensis", "O. brachyantha", "O. meyeriana", "L. hexandra"
    ]
    
    # Get sorted list of species
    existing_species = [s for s in target_order if s in df['Species'].unique()]
    other_species = [s for s in df['Species'].unique() if s not in target_order]
    final_order = existing_species + other_species

    # 2. Multidimensional Data Summary
    # Summary A: Total counts at Species + Chromosome level
    chrom_sum = df.groupby(['Species', 'Chromosome']).agg({
        'Cen_Count': 'sum',
        'NonCen_Count': 'sum'
    }).reset_index()

    # Summary B: Species totals (for bar charts)
    species_total = chrom_sum.groupby('Species').agg({
        'Cen_Count': 'sum',
        'NonCen_Count': 'sum'
    }).reindex(final_order).fillna(0)
    species_total['Total_SV'] = species_total['Cen_Count'] + species_total['NonCen_Count']

    # Summary C: Variation type distribution
    type_sum = df.groupby('VariationType').agg({
        'Cen_Count': 'sum',
        'NonCen_Count': 'sum'
    })
    type_perc = type_sum.div(type_sum.sum(axis=0), axis=1) * 100

    
    # ------------------ Plotting Section ------------------
    # Fig 1: Paired Heatmaps
    cen_matrix = chrom_sum.pivot(index='Species', columns='Chromosome', values='Cen_Count').fillna(0).reindex(final_order)
    arm_matrix = chrom_sum.pivot(index='Species', columns='Chromosome', values='NonCen_Count').fillna(0).reindex(final_order)

    fig1 = plt.figure(figsize=(7, 4.5)) 
    gs1 = gridspec.GridSpec(1, 2, width_ratios=[1, 1], wspace=0.25)
    
    # Left Heatmap: Centromere
    ax_cen = fig1.add_subplot(gs1[0])
    sns.heatmap(cen_matrix, annot=True, fmt=".0f", cmap="YlOrRd", ax=ax_cen, cbar=False, annot_kws={"size": 5})
    ax_cen.set_title('Centromere SVs (Cen)')
    plt.setp(ax_cen.get_xticklabels(), rotation=45, ha='right')
    
    # Right Heatmap: Chromosome Arms
    ax_arm = fig1.add_subplot(gs1[1])
    sns.heatmap(arm_matrix, annot=True, fmt=".0f", cmap="YlGnBu", ax=ax_arm, cbar=False, annot_kws={"size": 5})
    ax_arm.set_title('Chromosome Arm SVs (Non-Cen)')
    ax_arm.set_yticks([])
    ax_arm.set_ylabel('')
    plt.setp(ax_arm.get_xticklabels(), rotation=45, ha='right')
    
    fig1.savefig(os.path.join(output_dir, 'Figure1_Paired_Heatmap.pdf'), bbox_inches='tight')

    # Fig 2: Combined Statistics
    fig2 = plt.figure(figsize=(11, 4))
    gs2 = gridspec.GridSpec(1, 3, wspace=0.45)

    # A. Correlation (Reserved space)
    # ax_a = fig2.add_subplot(gs2[0])
    # sns.regplot(data=oryza_df, x='NonCen_Count', y='Cen_Count', ax=ax_a, scatter_kws={'alpha':0.5, 's':15}, line_kws={'color':'red', 'linewidth':1})
    # ax_a.set_title(f'A. Correlation (Oryza only)\nPearson r={corr_p_o:.2f}, Spearman ρ={corr_s_o:.2f}', fontsize=7)

    # B. Species Variation Load (with internal value labels)
    ax_b = fig2.add_subplot(gs2[1])
    spec_plot_data = species_total[['Cen_Count', 'NonCen_Count']].iloc[::-1]
    spec_plot_data.plot(kind='barh', stacked=True, ax=ax_b, color=['#ff96a4', '#8ecfc9'], legend=False, width=0.7)
    
    for container in ax_b.containers:
        labels = [f'{int(v)}' if v > 0 else "" for v in container.datavalues]
        ax_b.bar_label(container, labels=labels, label_type='center', fontsize=5, fontweight='bold')
        
    totals = spec_plot_data.sum(axis=1)
    for i, total in enumerate(totals):
        if total > 0: 
            ax_b.text(total + 1, i, f'{int(total)}', va='center', fontsize=6, fontweight='bold')
            
    ax_b.set_title('B. SV Load by Species')
    ax_b.set_xlabel('Count of SV Events')
    ax_b.set_xlim(0, totals.max() * 1.15)

    # C. Variation Type Composition (with percentage labels)
    ax_c = fig2.add_subplot(gs2[2])
    type_perc.T.plot(kind='bar', stacked=True, ax=ax_c, color=['#7fb3d5', '#f7dc6f', '#82e0aa'], width=0.6)
    
    for container in ax_c.containers:
        labels = [f'{v:.1f}%' if v > 3 else "" for v in container.datavalues]
        ax_c.bar_label(container, labels=labels, label_type='center', fontsize=5)
        
    ax_c.set_title('C. SV Type Composition (%)')
    ax_c.legend(title="Type", fontsize=5, loc='upper left', bbox_to_anchor=(1, 1))
    
    plt.tight_layout()
    fig2.savefig(os.path.join(output_dir, 'Figure2_Statistics_Combined.pdf'), bbox_inches='tight')
    
    print(f"Analysis Complete!")
    # Note: report_path variable was used in original script print but not defined in the snippet provided
    print(f"1. Vector Charts saved in: {output_dir}/Figure1 & Figure2")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, help="Input TSV file")
    parser.add_argument('-o', '--output_dir', required=True, help="Output directory")
    args = parser.parse_args()
    run_sv_master_analysis(args.input, args.output_dir)