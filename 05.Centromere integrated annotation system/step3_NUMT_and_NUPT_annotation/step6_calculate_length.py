import pandas as pd
import argparse
import sys
import subprocess
import os

def check_bedtools_installed():
    """Check if bedtools is installed and accessible in the system path."""
    try:
        # Execute 'bedtools --version' and redirect output to DEVNULL to maintain a clean terminal
        subprocess.run(['bedtools', '--version'], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False

def run_bedtools_summary(input_file, output_summary_file):
    """
    Merge intervals using bedtools merge, then calculate total length per chromosome using pandas.

    Args:
        input_file (str): Path to the input BED file.
        output_summary_file (str): Path to the output summary file for chromosomal lengths.
    """
    # --- Step 0: Dependency Check ---
    if not check_bedtools_installed():
        print("Error: 'bedtools' not found or not executable.", file=sys.stderr)
        print("Ensure bedtools is properly installed and its path is added to the system's PATH environment variable.", file=sys.stderr)
        sys.exit(1)
        
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found.", file=sys.stderr)
        sys.exit(1)
        
    # Define a temporary intermediate filename
    merged_bed_file = f"{os.path.splitext(input_file)[0]}.merged.tmp.bed"

    try:
        # --- Step 1: Execute bedtools merge ---
        print(f"[1/4] Merging intervals using bedtools merge...")
        bedtools_cmd = ['bedtools', 'merge', '-i', input_file]
        
        # Open a file handle to write the bedtools output
        with open(merged_bed_file, 'w') as f_out:
            # Execute command and redirect standard output to the specified file
            subprocess.run(bedtools_cmd, stdout=f_out, check=True)
        print(f"-> Merged intervals temporarily saved to: {merged_bed_file}")

        # --- Step 2: Load and process the merged file using pandas ---
        print(f"\n[2/4] Loading merged file and calculating lengths...")
        # Define column names for the merged file
        merged_columns = ['chrom', 'start', 'end']
        df_merged = pd.read_csv(
            merged_bed_file, 
            sep='\t', 
            header=None, 
            names=merged_columns
        )
        
        # Calculate the length of each merged interval
        df_merged['length'] = df_merged['end'] - df_merged['start']
        
        # Group by chromosome and aggregate the total length
        length_stats = df_merged.groupby('chrom')['length'].sum().reset_index()
        
        # Rename columns to meet output specifications
        length_stats.columns = ['chromosome', 'total_length']
        print("-> Length calculation complete.")

        # --- Step 3: Save final summary results ---
        print(f"\n[3/4] Saving summary results to: {output_summary_file} ...")
        length_stats.to_csv(
            output_summary_file, 
            sep='\t', 
            header=True, # Include header for improved readability
            index=False
        )
        print("\nProcessing successful!")
        print(f"Summary file saved to: {output_summary_file}")

    except subprocess.CalledProcessError as e:
        print(f"Error: 'bedtools merge' execution failed. Return code: {e.returncode}", file=sys.stderr)
        print(f"Verify that the input file '{input_file}' is in a valid BED format.", file=sys.stderr)
    except Exception as e:
        print(f"An error occurred during processing: {e}", file=sys.stderr)
    finally:
        # --- Step 4: Cleanup temporary files ---
        if os.path.exists(merged_bed_file):
            print(f"\n[4/4] Cleaning up temporary file: {merged_bed_file} ...")
            os.remove(merged_bed_file)


def main():
    """
    Main function to parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Merge BED intervals using bedtools merge and calculate the total coverage length per chromosome.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '-i', "--input_file", 
        help="Input BED filename.\nExample: my_regions.bed"
    )
    parser.add_argument(
        '-o', "--output_file", 
        help="Output length summary filename (two columns, tab-separated).\nExample: coverage_summary.txt"
    )
    
    args = parser.parse_args()
    run_bedtools_summary(args.input_file, args.output_file)

if __name__ == "__main__":
    main()