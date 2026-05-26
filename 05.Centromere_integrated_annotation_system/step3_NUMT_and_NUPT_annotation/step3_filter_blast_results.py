import pandas as pd
import argparse
import sys

def filter_blast_results(blast_file, output_file, identity_threshold=85.0, evalue_threshold=1e-5):
    """
    Filter BLAST outfmt 6 results using pandas.

    Args:
        blast_file (str): Input BLAST filename (outfmt 6).
        output_file (str): Output filtered filename.
        identity_threshold (float): Minimum threshold for percent identity.
        evalue_threshold (float): Maximum threshold for E-value.
    """
    # Define standard column names for BLAST outfmt 6 for reference
    column_names = [
        'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
    ]

    try:
        print(f"Reading BLAST file: {blast_file} ...")
        # Read the tab-separated file using pandas and specify column names
        df = pd.read_csv(
            blast_file, 
            sep='\t',         # Tab separator
            header=None,      # File contains no header row
            names=column_names # Assign defined column names
        )
        
        num_original_records = len(df)
        print(f"-> Read complete; found {num_original_records} alignment records.")

    except FileNotFoundError:
        print(f"Error: Input file '{blast_file}' not found.", file=sys.stderr)
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print(f"Error: Input file '{blast_file}' is empty and cannot be processed.", file=sys.stderr)
        sys.exit(1)

    # --- Core Filtering Logic ---
    print(f"\nStarting filtration...")
    print(f"Filter criterion 1: identity (pident) > {identity_threshold}")
    print(f"Filter criterion 2: E-value (evalue) < {evalue_threshold}")

    # Use boolean indexing for filtration; '&' represents the 'AND' operator
    filtered_df = df[
        (df['pident'] > identity_threshold) & 
        (df['evalue'] < evalue_threshold)
    ]
    
    num_filtered_records = len(filtered_df)
    print(f"-> Filtration complete; {num_filtered_records} records remaining.")

    # --- Save Results ---
    try:
        print(f"\nSaving filtered results to: {output_file} ...")
        # Save the filtered DataFrame as a new tab-separated file
        # index=False: Do not save pandas row indices
        # header=False: Do not write column names to maintain consistency with the original outfmt 6 format
        filtered_df.to_csv(output_file, sep='\t', index=False, header=False)
        
        print("-" * 30)
        print("Processing successful!")
        print(f"Original record count: {num_original_records}")
        print(f"Filtered record count: {num_filtered_records}")
        print(f"Results saved to: {output_file}")
        print("-" * 30)
        
        # Print the first 5 rows as a preview
        print("Filtered results preview (first 5 rows):")
        # Using to_string() for better output formatting
        print(filtered_df.head().to_string(index=False, header=False))

    except Exception as e:
        print(f"Error: An error occurred while saving the file: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    """
    Main function to parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Filter BLAST outfmt 6 result files using pandas.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '-i', "--blast_file", 
        help="Input BLAST result file (outfmt 6).\nExample: blast.txt"
    )
    parser.add_argument(
        '-o', "--output_file", 
        help="Output filtered filename.\nExample: filtered_blast.txt"
    )
    parser.add_argument(
        "--identity", 
        type=float, 
        default=85.0, 
        help="Minimum threshold for Identity (default: 85.0)"
    )
    parser.add_argument(
        "-e", "--evalue", 
        type=float, 
        default=1e-5, 
        help="Maximum threshold for E-value (default: 1e-5)"
    )
    
    args = parser.parse_args()
    
    filter_blast_results(args.blast_file, args.output_file, args.identity, args.evalue)

if __name__ == "__main__":
    main()