#!/usr/bin/env python3
"""
Generate AT content scatter plots and perform permutation tests.
Usage: python script.py --input data.tsv --output scatter_plot.png
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats


def permutation_test(group1, group2, n_permutations=10000):
    """
    Performs a permutation test to compare two groups of data.
    
    Parameters:
        group1 (array-like): First data group.
        group2 (array-like): Second data group.
        n_permutations (int): Number of permutations.
    
    Returns:
        float: Calculated p-value.
    """
    # Calculate the observed difference in means
    observed_diff = np.mean(group1) - np.mean(group2)
    
    # Combine the two groups
    combined = np.concatenate([group1, group2])
    n1 = len(group1)
    
    # Execute permutation process
    count = 0
    for _ in range(n_permutations):
        # Randomly shuffle the data
        np.random.shuffle(combined)
        
        # Split into new groups
        perm_group1 = combined[:n1]
        perm_group2 = combined[n1:]
        
        # Calculate the difference in means for the permuted groups
        perm_diff = np.mean(perm_group1) - np.mean(perm_group2)
        
        # Count occurrences of extreme values (two-tailed test)
        if abs(perm_diff) >= abs(observed_diff):
            count += 1
    
    # Calculate p-value
    p_value = count / n_permutations
    
    return p_value


def get_significance_symbol(p_value):
    """
    Returns the significance symbol corresponding to a p-value.
    
    Parameters:
        p_value (float): p-value.
    
    Returns:
        str: Significance symbol (*, **, ***, or ns).
    """
    if p_value < 0.001:
        return '***'
    elif p_value < 0.01:
        return '**'
    elif p_value < 0.05:
        return '*'
    else:
        return 'ns'


def plot_at_content_with_significance(input_file, output_file):
    """
    Plots AT content scatter plots and annotates significance.
    
    Parameters:
        input_file (str): Path to the input TSV file.
        output_file (str): Path to the output image file.
    """
    # Dictionary mapping species codes to full scientific names
    species_name_map = {
        'AA_Ojap': 'O. sativa ssp. japonica',
        'AA_Oind': 'O. sativa ssp. indica',
        'AA_Ogla': 'O. glaberrima',
        'AA_Oruf': 'O. rufipogon',
        'AA_Oniv': 'O. nivara',
        'AA_Olon': 'O. longistaminata',
        'AA_Oglu': 'O. glumaepatula',
        'BB_Opun': 'O. punctata',
        'CC_Ooff': 'O. officinalis',
        'EE_Oaus': 'O. australiensis',
        'FF_Obra': 'O. brachyantha',
        'GG_Omey': 'O. meyeriana'
    }
    
    # Define the standardized species order
    species_order = [
        'AA_Ojap', 'AA_Oind', 'AA_Ogla', 'AA_Oruf', 'AA_Oniv', 'AA_Olon',
        'AA_Oglu', 'BB_Opun', 'CC_Ooff', 'EE_Oaus', 'FF_Obra', 'GG_Omey'
    ]
    
    # Load data
    df = pd.read_csv(input_file, sep='\t')
    
    print(f"Loading data: {len(df)} rows")
    print(f"Number of species: {df['Species'].nunique()}")
    print(f"Number of types: {df['Type'].nunique()}")
    print()
    
    # Retrieve unique species and sort based on the predefined order
    all_species = df['Species'].unique()
    # Retain species that exist in the data and follow the predefined order
    species_list = [s for s in species_order if s in all_species]
    # Add species present in the data but missing from the predefined order (append to the end)
    species_list.extend([s for s in all_species if s not in species_order])
    
    # Reverse the list to display the first species at the top of the y-axis
    species_list = species_list[::-1]
    
    type_list = sorted(df['Type'].unique())
    
    print(f"Species list: {species_list}")
    print(f"Type list: {type_list}")
    print()
    
    # Set plotting style
    sns.set_style("whitegrid")
    plt.rcParams['font.size'] = 10
    
    # Initialize the figure - adjusted for a narrower/longer aspect ratio
    fig, ax = plt.subplots(figsize=(8, max(10, len(species_list) * 0.8)))
    
    # Define color palette for different types
    colors = ['#3498db', '#e74c3c']  # Blue and Red
    if len(type_list) > 2:
        colors = sns.color_palette("husl", n_colors=len(type_list))
    type_colors = {type_name: colors[i] for i, type_name in enumerate(type_list)}
    
    # List to store statistical results
    significance_results = []
    
    # Process each species
    y_positions = {}
    offset_width = 0.3  # Offset distance between different types
    
    for i, species in enumerate(species_list):
        species_data = df[df['Species'] == species]
        
        type_data = {}
        n_types = len(type_list)
        
        for j, type_name in enumerate(type_list):
            type_subset = species_data[species_data['Type'] == type_name]
            if len(type_subset) > 0:
                at_values = type_subset['Average_AT_Content'].values
                type_data[type_name] = at_values
                
                # Calculate y-axis offset for side-by-side plotting of types
                if n_types == 2:
                    y_offset = (j - 0.5) * offset_width
                else:
                    y_offset = (j - (n_types - 1) / 2) * offset_width / (n_types - 1)
                
                # Add light random jitter to prevent point overlap
                y_jitter = np.random.normal(0, 0.02, size=len(at_values))
                y_coords = i + y_offset + y_jitter
                
                # Plot scatter points
                ax.scatter(at_values, y_coords, 
                          c=[type_colors[type_name]], 
                          label=type_name if i == 0 else "",
                          alpha=0.6, 
                          s=50,
                          edgecolors='black',
                          linewidths=0.5)
                
                # Plot mean lines
                mean_val = np.mean(at_values)
                ax.vlines(mean_val, i + y_offset - 0.08, i + y_offset + 0.08, 
                         colors=type_colors[type_name], 
                         linewidths=2,
                         linestyles='solid')
        
        y_positions[species] = i
        
        # If exactly two types are present, perform a permutation test
        if len(type_data) == 2:
            types = list(type_data.keys())
            group1 = type_data[types[0]]
            group2 = type_data[types[1]]
            
            # Execute permutation test
            p_value = permutation_test(group1, group2, n_permutations=10000)
            sig_symbol = get_significance_symbol(p_value)
            
            significance_results.append({
                'Species': species,
                'Type1': types[0],
                'Type2': types[1],
                'P_value': p_value,
                'Significance': sig_symbol
            })
            
            print(f"{species}: {types[0]} vs {types[1]}")
            print(f"  Means: {np.mean(group1):.6f} vs {np.mean(group2):.6f}")
            print(f"  p-value: {p_value:.4f} ({sig_symbol})")
            print()
            
            # Annotate significance on the plot
            x_max = max(np.max(group1), np.max(group2))
            x_min = min(np.min(group1), np.min(group2))
            x_range = x_max - x_min
            
            # Set significance mark position
            x_sig = x_max + x_range * 0.08
            
            # Determine y-offsets for the comparison brackets
            if len(type_list) == 2:
                y_offset_left = -0.5 * offset_width
                y_offset_right = 0.5 * offset_width
            else:
                type_indices = {t: idx for idx, t in enumerate(type_list)}
                y_offset_left = (type_indices[types[0]] - (len(type_list) - 1) / 2) * offset_width / (len(type_list) - 1)
                y_offset_right = (type_indices[types[1]] - (len(type_list) - 1) / 2) * offset_width / (len(type_list) - 1)
            
            # Draw connection lines
            ax.plot([x_sig, x_sig], [i + y_offset_left, i + y_offset_right], 'k-', linewidth=1.5)
            ax.text(x_sig + x_range * 0.02, i, sig_symbol, 
                   ha='left', va='center', fontsize=12, fontweight='bold')
    
    # Configure y-axis
    ax.set_yticks(range(len(species_list)))
    # Use scientific names for display
    species_display_names = [species_name_map.get(s, s) for s in species_list]
    ax.set_yticklabels(species_display_names, fontsize=10, style='italic')
    ax.set_ylabel('Species', fontsize=12, fontweight='bold')
    
    # Configure x-axis
    ax.set_xlabel('Average AT Content', fontsize=12, fontweight='bold')
    
    # Set plot title
    ax.set_title('AT Content Comparison Across Species and Types\nwith Permutation Test', 
                fontsize=14, fontweight='bold', pad=20)
    
    # Add legend
    handles, labels = ax.get_legend_handles_labels()
    if len(handles) > 0:
        ax.legend(handles, labels, title='Type', loc='best', frameon=True)
    
    # Add significance key at the bottom
    sig_text = "Significance: *** p < 0.001, ** p < 0.01, * p < 0.05, ns: not significant"
    fig.text(0.5, 0.02, sig_text, ha='center', fontsize=9, style='italic')
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0.04, 1, 1])
    
    # Save the figure
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nPlot saved to: {output_file}")
    
    plt.close()
    
    # Return results as a DataFrame
    return pd.DataFrame(significance_results)


def main():
    parser = argparse.ArgumentParser(
        description='Generate AT content scatter plots and perform permutation tests.'
    )
    parser.add_argument(
        '--input',
        required=True,
        help='Path to the input TSV file'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Path to the output image file (e.g., plot.png)'
    )
    
    args = parser.parse_args()
    
    # Execute plotting and statistical analysis
    results_df = plot_at_content_with_significance(args.input, args.output)
    
    # Save statistical results to a file
    if len(results_df) > 0:
        stats_output = args.output.rsplit('.', 1)[0] + '_statistics.tsv'
        results_df.to_csv(stats_output, sep='\t', index=False)
        print(f"Statistical results saved to: {stats_output}")


if __name__ == '__main__':
    main()