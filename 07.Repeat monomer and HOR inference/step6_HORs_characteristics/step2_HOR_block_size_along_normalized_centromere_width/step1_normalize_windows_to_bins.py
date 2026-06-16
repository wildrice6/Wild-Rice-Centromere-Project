#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import argparse
import sys

def analyze_window_blocks(input_file, output_file):
    """
    Partitions window records into 100 bins and calculates the mean max_block_size for each bin.
    If the total number of records is less than 100, the data is centered and zero-padded at both ends.

    Args:
        input_file (str): Path to the input TSV file.
        output_file (str): Path to the output TSV file.
    """
    try:
        # 1. Read input file
        print(f"--- Reading input file: {input_file} ---")
        df = pd.read_csv(input_file, sep='\s+')
        
        if len(df.columns) < 2:
            print(f"Error: Input file '{input_file}' contains fewer than 2 columns.", file=sys.stderr)
            sys.exit(1)
        block_size_col = df.columns[1]
        
        total_windows = len(df)
        print(f"--- File read successfully; {total_windows} records identified. ---")
        
        final_averages = []
        
        # 2. Check total record count and process accordingly
        if total_windows < 100:
            # --- Modification Start ---
            print("--- Total record count is less than 100; centering data with zero-padding. ---")
            
            # Retrieve existing data
            existing_data = df[block_size_col].tolist()
            
            # Calculate total required zero-padding
            total_padding = 100 - total_windows
            
            # Calculate left padding count (using integer division)
            left_padding_count = total_padding // 2
            
            # Calculate right padding count
            # This handles cases where total_padding is odd; the extra zero is placed on the right
            right_padding_count = total_padding - left_padding_count
            
            # Create padding lists
            left_padding = [0] * left_padding_count
            right_padding = [0] * right_padding_count
            
            # Concatenate three components: left padding + actual data + right padding
            final_averages = left_padding + existing_data + right_padding
            # --- Modification End ---
            
        else:
            print("--- Total record count is 100 or greater; proceeding with binning. ---")
            # 3. Use numpy.array_split to partition the DataFrame into 100 bins
            chunks = np.array_split(df, 100)
            
            # 4. Iterate through each bin and calculate the mean
            for chunk in chunks:
                average_size = chunk[block_size_col].mean()
                final_averages.append(average_size)

        # 5. Construct the final output DataFrame
        bin_numbers = range(1, 101)
        result_df = pd.DataFrame({
            'bin_number': bin_numbers,
            'average_max_block_size': final_averages
        })

        # 6. Save results to file
        print(f"--- Writing results to: {output_file} ---")
        result_df.to_csv(output_file, sep='\t', index=False, float_format='%.4f')
        print("--- Task complete! ---")

    except FileNotFoundError:
        print(f"Error: Input file not found - {input_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred during processing: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    """Main function for parsing command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Partitions window records into 100 bins and calculates the mean max_block_size for each bin. Centers the data if there are fewer than 100 records."
    )
    parser.add_argument(
        '--input',
        required=True,
        help="Path to the input TSV file."
    )
    parser.add_argument(
        '--output',
        required=True,
        help="Path to the output TSV file."
    )
    
    args = parser.parse_args()
    
    analyze_window_blocks(args.input, args.output)

if __name__ == '__main__':
    main()