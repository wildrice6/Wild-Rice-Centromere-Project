#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import argparse
import sys

def aggregate_bins_per_assembly(input_file, output_file):
    """
    Groups data by 'assembly' and 'bin_number' to calculate the mean of 
    'average_max_block_size' for each bin across all chromosomes within each assembly.

    Args:
        input_file (str): Path to the input TSV file.
        output_file (str): Path to the output TSV file.
    """
    try:
        # 1. Read input file
        print(f"--- Reading input file: {input_file} ---")
        df = pd.read_csv(input_file, sep='\s+')

        # --- Modification: Add 'assembly' to the required columns check ---
        required_columns = ['assembly', 'bin_number', 'average_max_block_size']
        if not all(col in df.columns for col in required_columns):
            print(f"Error: Input file must contain the following columns: {required_columns}", file=sys.stderr)
            sys.exit(1)

        print("--- File read successfully; beginning data aggregation by assembly ---")
        
        # --- Core modification: Group by both 'assembly' and 'bin_number' simultaneously ---
        # Passing a list ['assembly', 'bin_number'] to groupby()
        # This creates a group for each unique (assembly, bin_number) combination
        aggregated_data = df.groupby(['assembly', 'bin_number'])['average_max_block_size'].mean()

        # Following the .mean() operation, the result is a Pandas Series with a MultiIndex (assembly, bin_number)
        # Using .reset_index() converts it back into a DataFrame containing 'assembly' and 'bin_number' columns
        result_df = aggregated_data.reset_index()

        # 3. (Optional) Rename the aggregated column for clarity
        result_df.rename(columns={'average_max_block_size': 'mean_average_max_block_size'}, inplace=True)
        
        # 4. Save results to a new TSV file
        print(f"--- Writing aggregated results to: {output_file} ---")
        result_df.to_csv(output_file, sep='\t', index=False, float_format='%.4f')
        print("--- Task complete! ---")

    except FileNotFoundError:
        print(f"Error: Input file not found - {input_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred during processing: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    """Main function for parsing command-line arguments"""
    # --- Modification: Update script description ---
    parser = argparse.ArgumentParser(
        description="Aggregates data by assembly, calculating the mean average_max_block_size for each bin across all chromosomes within each assembly."
    )
    parser.add_argument(
        '--input',
        required=True,
        help="Path to the input TSV file; must include 'assembly', 'bin_number', and 'average_max_block_size' columns."
    )
    parser.add_argument(
        '--output',
        required=True,
        help="Path to the output aggregated TSV file."
    )
    
    args = parser.parse_args()
    
    aggregate_bins_per_assembly(args.input, args.output)

if __name__ == '__main__':
    main()