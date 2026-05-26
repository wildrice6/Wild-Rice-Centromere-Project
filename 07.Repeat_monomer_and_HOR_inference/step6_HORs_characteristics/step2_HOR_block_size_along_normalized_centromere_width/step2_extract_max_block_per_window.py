#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
find_max_block.py

This script reads a TSV file containing HOR alignment information.
It parses window names within the A_pos and B_pos columns and 
identifies the maximum block_size associated with each unique window.

Input file format (tab-separated):
HOR_ID  block_size  ...  A_pos     B_pos     ...
1       3           ...  2;3;4     3;4;5     ...
...

Output file format (tab-separated):
window_name     max_block_size
2               3
3               4
4               4
...

Usage:
python find_max_block.py --input <input_file.tsv> --output <output_file.tsv>
"""

import argparse
import sys
import csv

def find_max_block_for_windows(input_path, output_path):
    """
    Processes the input file to identify the maximum block_size for each window.

    Args:
        input_path (str): Path to the input TSV file.
        output_path (str): Path to the output TSV file.
    """
    # Dictionary to store the maximum block size encountered for each window
    # Structure: { "window_name": max_block_size }
    window_max_blocks = {}

    print(f"Reading input file: {input_path}...")
    try:
        with open(input_path, 'r', newline='') as f_in:
            # Use csv.reader for robust handling of the TSV format
            reader = csv.reader(f_in, delimiter='\t')
            
            # Read header and identify indices for required columns to enhance robustness
            header = next(reader)
            try:
                # Retrieve the index positions for required columns
                block_size_idx = header.index('block_size')
                a_pos_idx = header.index('A_pos')
                b_pos_idx = header.index('B_pos')
            except ValueError as e:
                print(f"Error: Input file header is missing required columns: {e}", file=sys.stderr)
                sys.exit(1)

            # Iterate through each row of data in the file
            for line_num, row in enumerate(reader, 2): # Counter starts from line 2
                # Skip empty rows or rows with incorrect formatting
                if not row or len(row) <= max(block_size_idx, a_pos_idx, b_pos_idx):
                    print(f"Warning: Line {line_num} is incorrectly formatted or empty; skipped.", file=sys.stderr)
                    continue
                
                try:
                    # Extract the block_size for the current row
                    current_block_size = int(row[block_size_idx])
                    
                    # Extract the window lists from A_pos and B_pos
                    a_windows = row[a_pos_idx].split(';')
                    b_windows = row[b_pos_idx].split(';')
                    
                    # Combine all windows
                    all_windows_for_row = a_windows + b_windows
                    
                    # Iterate through all windows involved in the current row
                    for window in all_windows_for_row:
                        # Skip potential empty strings resulting from terminal semicolons
                        if not window:
                            continue
                        
                        # Core logic: Update the maximum block size for the respective window
                        # .get(window, 0) retrieves the current maximum, returning 0 if the window is not present
                        if current_block_size > window_max_blocks.get(window, 0):
                            window_max_blocks[window] = current_block_size

                except (ValueError, IndexError) as e:
                    print(f"Warning: Error processing data at line {line_num}; skipped. Error: {e}", file=sys.stderr)
                    continue

    except FileNotFoundError:
        print(f"Error: Input file not found: '{input_path}'", file=sys.stderr)
        sys.exit(1)
        
    print(f"File processing complete. {len(window_max_blocks)} unique windows identified.")

    # --- Write to output file ---
    print(f"Writing results to output file: {output_path}...")
    try:
        with open(output_path, 'w', newline='') as f_out:
            writer = csv.writer(f_out, delimiter='\t')
            
            # Write header
            writer.writerow(['window_name', 'max_block_size'])
            
            # Sort window names numerically to ensure organized output
            # key=int ensures that '10' is sequenced after '9'
            try:
                sorted_windows = sorted(window_max_blocks.keys(), key=int)
            except ValueError:
                print("Warning: Window names contain non-integers; performing lexicographical sort.", file=sys.stderr)
                sorted_windows = sorted(window_max_blocks.keys())

            # Iterate through sorted windows and write to file
            for window_name in sorted_windows:
                max_size = window_max_blocks[window_name]
                writer.writerow([window_name, max_size])

    except IOError as e:
        print(f"Error: Unable to write to output file '{output_path}'. Error details: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Processing complete. Results successfully saved to {output_path}")


def main():
    """
    Main function for parsing command-line arguments and invoking the processing function.
    """
    parser = argparse.ArgumentParser(
        description="Identify the maximum block_size for each individual window.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument(
        "--input",
        required=True,
        help="Path to the input TSV file."
    )
    
    parser.add_argument(
        "--output",
        required=True,
        help="Path to the output TSV file."
    )
    
    args = parser.parse_args()
    
    find_max_block_for_windows(args.input, args.output)


if __name__ == "__main__":
    main()