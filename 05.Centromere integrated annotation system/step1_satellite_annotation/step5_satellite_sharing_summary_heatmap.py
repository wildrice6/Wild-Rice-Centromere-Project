import re
from collections import defaultdict
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

def extract_chr(filename):
    """
    Extract the 'Chr' followed by two digits pattern from the filename 
    using regular expressions.
    Example: 'AA_Ogla_hap1.sat.CEN155.Chr01' -> 'Chr01'
    """
    match = re.search(r'Chr(\d{2})', filename)
    if match:
        return match.group(0) # Returns the entire matched string, e.g., "Chr01"
    return None

def summarize_sharing_ratios(input_file, output_file):
    """
    Read the sharing ratio file, summarize by chromosome pairs, 
    and calculate the average sharing ratio.
    """
    # Use defaultdict to simplify the code; it automatically creates a default 
    # value (a list [0.0, 0]) when a key is accessed for the first time.
    # Structure: {(chrA, chrB): [ratio_sum, count]}
    summary_data = defaultdict(lambda: [0.0, 0])

    print(f"Reading input file: {input_file} ...")
    
    try:
        with open(input_file, 'r') as f_in:
            # Read and skip the header
            header = next(f_in)

            for line in f_in:
                # Use split() to handle potential multiple spaces or tabs
                parts = line.strip().split()
                if len(parts) < 7:
                    print(f"Warning: Skipping malformed line -> {line.strip()}")
                    continue

                file1 = parts[0]
                file2 = parts[1]
                
                try:
                    # Column 7 is SharingRatio (index 6)
                    ratio = float(parts[6])
                except (ValueError, IndexError):
                    print(f"Warning: Unable to parse SharingRatio, skipping line -> {line.strip()}")
                    continue

                chr1 = extract_chr(file1)
                chr2 = extract_chr(file2)

                # Ensure chromosome IDs are identified in both filenames
                if chr1 and chr2:
                    # Create a canonical key by sorting to ensure (Chr01, Chr03) 
                    # and (Chr03, Chr01) are treated as the same key.
                    canonical_key = tuple(sorted((chr1, chr2)))
                    
                    # Accumulate the sum of SharingRatio
                    summary_data[canonical_key][0] += ratio
                    # Increment the occurrence count for this chromosome pair
                    summary_data[canonical_key][1] += 1

    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.")
        return
        
    print("Data summarization complete. Calculating averages and writing to the output file...")

    # Sort results by chromosome ID for a cleaner output file
    sorted_keys = sorted(summary_data.keys())

    with open(output_file, 'w') as f_out:
        # Write the new table header
        f_out.write("Chr1\tChr2\tAverageSharingRatio\n")
        
        for key in sorted_keys:
            chr1, chr2 = key
            total_ratio, count = summary_data[key]
            
            # Calculate average, avoiding division by zero
            average_ratio = total_ratio / count if count > 0 else 0.0
            
            # Write result row, formatted to 4 decimal places
            f_out.write(f"{chr1}\t{chr2}\t{average_ratio:.4f}\n")

    print(f"Processing complete! Summary results saved to: {output_file}")

def plot_triangular_heatmap(df, figsize=(10, 8), cmap='Blues', 
                           title='Chromosomal Sharing Heatmap', 
                           annot=True, cbar_label='Average Sharing (%)'):
    """
    Generate a triangular heatmap to visualize average sharing ratios between chromosomes.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing 'chrA', 'chrB', and 'average_sharing' columns.
    figsize : tuple, default (10, 8)
        Dimensions of the figure.
    cmap : str, default 'Blues'  
        Color map for the heatmap.
    title : str
        Title of the plot.
    annot : bool, default True
        Whether to annotate cells with numerical values.
    cbar_label : str
        Label for the color bar.
        
    Returns:
    --------
    fig, ax : matplotlib figure and axes objects.
    """
    
    # Retrieve all unique chromosomes
    all_chrs = sorted(list(set(df['chrA'].tolist() + df['chrB'].tolist())))
    n_chrs = len(all_chrs)
    
    # Create a symmetric matrix
    matrix = np.zeros((n_chrs, n_chrs))
    
    # Map chromosome names to indices
    chr_to_idx = {chr_name: idx for idx, chr_name in enumerate(all_chrs)}
    
    # Populate the matrix
    for _, row in df.iterrows():
        i = chr_to_idx[row['chrA']]
        j = chr_to_idx[row['chrB']]
        matrix[i, j] = row['average_sharing']
        matrix[j, i] = row['average_sharing']  # Symmetric fill
    
    # Create a mask to display only the lower triangle (including diagonal)
    mask = np.tril(np.ones_like(matrix, dtype=bool), k=-1)
    
    # Initialize figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot the heatmap
    sns.heatmap(matrix, 
                mask=mask,
                annot=annot,
                cmap=cmap,
                square=True,
                xticklabels=all_chrs,
                yticklabels=all_chrs,
                cbar_kws={'label': cbar_label},
                ax=ax)
    
    # Set titles and axis labels
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    ax.set_xlabel('Chromosome', fontsize=12)
    ax.set_ylabel('Chromosome', fontsize=12)
    
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45)
    plt.yticks(rotation=0)
    
    # Adjust layout
    plt.tight_layout()
    # Save the figure in PDF and PNG formats
    pdf_filename = f"chr.sharing.pdf"
    png_filename = f"chr.sharing.png"
    
    fig.savefig(pdf_filename, format='pdf', dpi=300, bbox_inches='tight')
    fig.savefig(png_filename, format='png', dpi=300, bbox_inches='tight')
    
    return fig, ax


# --- Main Program Entry Point ---
if __name__ == "__main__":
    # --- Configure filenames here ---
    # Input filename (the original pairwise results)
    input_filename = "sharing_results_pairwise.tsv" 
    
    # Output filename (the generated summary file)
    output_filename = "sharing_ratio_summary_by_chr.tsv"
    
    # Execute the summarization function
    summarize_sharing_ratios(input_filename, output_filename)

    # Plotting
    temp = pd.read_csv('sharing_ratio_summary_by_chr.tsv', sep = '\t')
    temp.columns = ['chrA', 'chrB', 'average_sharing']
    temp['average_sharing'] = temp['average_sharing'] * 100
    plot_triangular_heatmap(temp)