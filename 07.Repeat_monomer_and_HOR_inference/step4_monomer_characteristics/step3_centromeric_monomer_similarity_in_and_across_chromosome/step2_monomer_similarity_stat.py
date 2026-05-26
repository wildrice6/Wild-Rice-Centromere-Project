import argparse
import pandas as pd
import numpy as np
import os

def permutation_test(group_intra, group_inter, n_permutations=10000):
    """
    Performs a permutation test on the aggregated average similarity data.
    """
    # Calculate the observed difference in means
    obs_diff = np.abs(np.mean(group_intra) - np.mean(group_inter))
    combined = np.concatenate([group_intra, group_inter])
    n_intra = len(group_intra)
    
    count = 0
    # Randomly shuffle labels to generate the null distribution
    for _ in range(n_permutations):
        np.random.shuffle(combined)
        new_intra = combined[:n_intra]
        new_inter = combined[n_intra:]
        new_diff = np.abs(np.mean(new_intra) - np.mean(new_inter))
        if new_diff >= obs_diff:
            count += 1
            
    return count / n_permutations

def calculate_ani_stats(input_path, output_path):
    try:
        # 0. Extract the assembly name from the filename
        file_name = os.path.basename(input_path)
        assembly_name = file_name.split('.')[0]

        # 1. Load the original TSV file
        df = pd.read_csv(input_path, sep='\t')

        # 2. Extract chromosome IDs and normalize chromosome pairs
        df['chr1'] = df['seq1'].str.split(':').str[0]
        df['chr2'] = df['seq2'].str.split(':').str[0]
        # Ensure chrA <= chrB to standardize pair representation
        df['chrA'] = np.where(df['chr1'] <= df['chr2'], df['chr1'], df['chr2'])
        df['chrB'] = np.where(df['chr1'] <= df['chr2'], df['chr2'], df['chr1'])
        df['type'] = np.where(df['chrA'] == df['chrB'], 'within', 'between')

        # --- A. Generate Detailed Statistics Table ---
        # Group by chromosome pairs and calculate the mean similarity for each pair
        detail_result = df.groupby(['chrA', 'chrB', 'type'])['ANI'].agg(['count', 'mean']).reset_index()
        detail_result.columns = ['chrA', 'chrB', 'type', 'count', 'avg_similarity']
        
        # Save the detailed statistics table
        detail_result.to_csv(output_path, sep='\t', index=False)

        # --- B. Perform Permutation Test Based on Aggregated Data ---
        # Extract avg_similarity for intra- and inter-chromosome groups
        intra_group_means = detail_result[detail_result['type'] == 'within']['avg_similarity'].values
        inter_group_means = detail_result[detail_result['type'] == 'between']['avg_similarity'].values

        if len(intra_group_means) > 0 and len(inter_group_means) > 0:
            # Calculate overall means (mean of the pair-wise averages)
            overall_intra_mean = np.mean(intra_group_means)
            overall_inter_mean = np.mean(inter_group_means)
            
            print(f"Detected {len(intra_group_means)} intra-chromosome records and {len(inter_group_means)} inter-chromosome records.")
            print("Performing permutation test on the means...")
            
            # Execute the permutation test
            p_val = permutation_test(intra_group_means, inter_group_means, n_permutations=10000)
            
            summary_data = {
                'inter_chromosome_mean': [overall_inter_mean],
                'intra_chromosome_mean': [overall_intra_mean],
                'p_value': [p_val],
                'assembly': [assembly_name]
            }
            summary_df = pd.DataFrame(summary_data)
            
            # Define output path for the summary statistics
            summary_output_path = output_path.replace('.tsv', '') + '.summary.tsv'
            summary_df.to_csv(summary_output_path, sep='\t', index=False)
            
            print(f"Detailed statistics saved to: {output_path}")
            print(f"Summary statistics saved to: {summary_output_path}")
        else:
            print("Error: Insufficient 'within' or 'between' records for statistical testing.")

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate mean ANI for chromosome pairs and perform permutation tests on the aggregated means.")
    parser.add_argument('--input', required=True, help="Path to the input TSV file")
    parser.add_argument('--output', required=True, help="Path to the output detailed TSV file")

    args = parser.parse_args()
    calculate_ani_stats(args.input, args.output)