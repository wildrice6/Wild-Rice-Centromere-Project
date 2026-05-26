#!/usr/bin/env python3

import argparse
import sys

# --- Classification Rule Configuration ---
# Rules are defined here for ease of future modification or expansion.
# Format: {"Classification Name": Typical Length}
CLASSIFICATION_RULES = {
    "CEN103": 103,
    "CEN126": 126,
    "CEN155": 155,
    "CEN165": 165,
    "CEN180": 180,
}
# Define length tolerance range
TOLERANCE = 5

def get_classification(length: int) -> str:
    """
    Returns the classification name based on the provided length and preset rules.

    Args:
        length (int): The length of the record.

    Returns:
        str: Classification name (e.g., 'CEN155') or 'atypical'.
    """
    # Iterate through all predefined rules
    for name, typical_len in CLASSIFICATION_RULES.items():
        # Check if the length falls within the range [typical - 5, typical + 5)
        if (typical_len - TOLERANCE) <= length < (typical_len + TOLERANCE):
            return name  # Return classification name immediately upon match
    
    # Return 'atypical' if no matches are found after the loop
    return "atypical"

def process_file(input_file: str, output_file: str):
    """
    Processes the input file, performs classification, and writes the results to an output file.

    Args:
        input_file (str): Path to the input TSV file.
        output_file (str): Path to the output TSV file.
    """
    print(f"Starting file processing: {input_file}")
    
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            # Read the input file line by line
            for line_num, line in enumerate(infile, 1):
                # Skip empty lines
                if not line.strip():
                    continue
                
                # Input is a tab-separated values (TSV) file
                parts = line.strip().split('\t')
                
                # --- Robustness Check ---
                # Verify that the line contains a sufficient number of columns
                if len(parts) < 7:
                    print(f"Warning (Line {line_num}): Column count is less than 7; writing line unchanged. Content: '{line.strip()}'")
                    outfile.write(line)
                    continue
                
                try:
                    # Extract length from the 7th column and convert to integer
                    record_length = int(parts[6])
                except ValueError:
                    print(f"Warning (Line {line_num}): 7th column is not a valid integer; writing line unchanged. Content: '{line.strip()}'")
                    outfile.write(line)
                    continue

                # --- Core Logic ---
                # 1. Determine the classification name
                new_type = get_classification(record_length)
                
                # 2. Replace content of the 4th column (index 3 in 0-based indexing)
                parts[3] = new_type
                
                # 3. Reconstruct the modified list into a tab-separated string
                output_line = "\t".join(parts) + "\n"
                
                # 4. Write to the output file
                outfile.write(output_line)
                
        print(f"Processing complete! Results successfully saved to: {output_file}")

    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred during file processing: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    """
    Main function to configure and parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Classify TSV records based on the length in the 7th column and update the 4th column identifier."
    )
    
    parser.add_argument('--input', '-i',
                        required=True,
                        help="Path to the input TSV file.")
                        
    parser.add_argument('--output', '-o',
                        required=True,
                        help="Path to the output classified TSV file.")
                        
    args = parser.parse_args()
    
    process_file(args.input, args.output)

if __name__ == "__main__":
    main()