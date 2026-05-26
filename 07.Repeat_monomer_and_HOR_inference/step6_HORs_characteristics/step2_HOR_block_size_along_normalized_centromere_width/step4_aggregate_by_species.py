#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import argparse
import sys

def aggregate_haplotypes(input_file, output_file):
    """
    Reads bin data containing haplotype information, aggregates haplotypes by species,
    and calculates the mean block size for each bin.

    Args:
        input_file (str): Path to the input TSV file.
                          Format: assembly, bin_number, mean_average_max_block_size
        output_file (str): Path to the output TSV file.
                           Format: species, bin_number, average_block_size
    """
    try:
        # 1. Read input TSV file
        # Use read_csv and specify tab as the delimiter
        df = pd.read_csv(input_file, sep='\t')
        
        # Verify that the input file contains the necessary columns
        required_columns = ['assembly', 'bin_number', 'mean_average_max_block_size']
        if not all(col in df.columns for col in required_columns):
            missing = set(required_columns) - set(df.columns)
            print(f"Error: Input file '{input_file}' is missing required columns: {', '.join(missing)}", file=sys.stderr)
            sys.exit(1)

    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred while reading the file: {e}", file=sys.stderr)
        sys.exit(1)

    # 2. Extract species names from the 'assembly' column
    # .str.split('_', n=1) splits only at the first underscore
    # .str[0] retrieves the first component of the resulting split
    df['species'] = df['assembly'].str.split('_', n=1).str[0]

    # 3. Group by 'species' and 'bin_number' and calculate the mean
    # - groupby(['species', 'bin_number']) creates the groups
    # - ['mean_average_max_block_size'] selects the column for aggregation
    # - .mean() calculates the average for each group
    # - .reset_index() converts the grouping keys (species, bin_number) from indices back into columns
    result_df = df.groupby(['species', 'bin_number'])['mean_average_max_block_size'].mean().reset_index()

    # 4. Rename the result column for clarity
    result_df.rename(columns={'mean_average_max_block_size': 'average_block_size'}, inplace=True)
    
    # 5. Save results to the output file
    try:
        # index=False indicates that the DataFrame index should not be written to the file
        result_df.to_csv(output_file, sep='\t', index=False)
        print(f"Processing complete. Results saved to '{output_file}'")
    except Exception as e:
        print(f"An error occurred while saving the file: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    # --- Configure Command-Line Argument Parsing ---
    parser = argparse.ArgumentParser(
        description="Aggregates haplotypes by species and calculates the mean block size for each bin.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '--input', 
        required=True, 
        help='Path to the input TSV file.\nColumns should be: assembly, bin_number, mean_average_max_block_size'
    )
    parser.add_argument(
        '--output', 
        required=True, 
        help='Path to the output TSV file.'
    )
    
    args = parser.parse_args()
    
    # --- Invoke Main Function ---
    aggregate_haplotypes(args.input, args.output)