import pandas as pd
import numpy as np
import argparse
import sys

def process_data(input_file, output_file):
    # 1. Load data
    print(f"Reading file: {input_file}...")
    try:
        df = pd.read_csv(input_file, sep='\t', header=None, 
                         names=['chr', 'label', 'pos', 'val'])
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)

    # 2. Process data by chromosome grouping
    result_list = []
    
    # Retrieve the list of unique chromosomes
    chromosomes = df['chr'].unique()
    
    for chrom in chromosomes:
        print(f"Processing chromosome: {chrom}")
        chr_df = df[df['chr'] == chrom].copy()
        n_rows = len(chr_df)
        
        if n_rows == 0:
            continue
            
        # 3. Create bin indices (dividing rows into 100 equal segments)
        # Map row indices to a range of 0-99
        chr_df['bin'] = np.floor(np.arange(n_rows) * 100 / n_rows).astype(int)
        
        # 4. Define aggregation logic
        # label: Assign 'CEN' if any entry in the group is 'CEN', otherwise assign 'ARM'
        def aggregate_label(x):
            return 'CEN' if 'CEN' in x.values else 'ARM'

        # Execute aggregation by grouping by bin
        binned = chr_df.groupby('bin').agg({
            'chr': 'first',
            'label': aggregate_label,
            'pos': 'mean',
            'val': 'mean'
        })
        
        result_list.append(binned)

    # 5. Concatenate results and save to file
    final_df = pd.concat(result_list)
    
    # Reorder columns to match the original TSV format
    final_df = final_df[['chr', 'label', 'pos', 'val']]
    
    print(f"Processing complete. Generated {len(final_df)} records (maximum 100 bins per chromosome).")
    print(f"Saving results to: {output_file}")
    final_df.to_csv(output_file, sep='\t', header=False, index=False)

def main():
    parser = argparse.ArgumentParser(description="Divide chromosomal signal data into 100 bins and perform aggregation.")
    parser.add_argument('--input', required=True, help="Path to the input TSV file")
    parser.add_argument('--output', required=True, help="Path to the output TSV file")
    
    args = parser.parse_args()
    process_data(args.input, args.output)

if __name__ == "__main__":
    main()