import argparse
import subprocess
import os
import sys

def run_blast_pipeline(query_file, target_file):
    """
    Executes the complete pipeline for makeblastdb and blastn, 
    and organizes database files into a separate folder.

    Args:
        query_file (str): Path to the query FASTA file.
        target_file (str): Path to the target/reference FASTA file.
    """
    print("--- Starting BLAST Pipeline Execution ---")

    # --- Step 0: Verify existence of input files and dependencies ---
    if not os.path.exists(query_file):
        print(f"Error: Query file '{query_file}' not found!", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(target_file):
        print(f"Error: Target file '{target_file}' not found!", file=sys.stderr)
        sys.exit(1)

    # --- Step 1: Prepare filenames and directory paths ---
    # Automatically generate folder name from target filename (e.g., NIP-T2T.fa -> NIP-T2T)
    target_basename = os.path.basename(target_file)
    db_folder = os.path.splitext(target_basename)[0]

    # Construct the full path prefix for database files. 
    # Database files will be located within db_folder, sharing the folder's name.
    # Example: NIP-T2T/NIP-T2T
    db_path_prefix = os.path.join(db_folder, db_folder)
    
    # Logic for output filename remains unchanged; results are saved in the current directory.
    query_base = os.path.splitext(os.path.basename(query_file))[0]
    output_file = f"{query_base}.{db_folder}.blast.txt"

    print(f"Query File: {query_file}")
    print(f"Target File: {target_file}")
    print(f"Database Directory: {db_folder}/")
    print(f"Database Path Prefix: {db_path_prefix}")
    print(f"Output File: {output_file}")
    print("-" * 20)

    # --- Step 2: Create the database directory ---
    print(f"\n[1/3] Creating database directory...")
    try:
        os.makedirs(db_folder, exist_ok=True)
        print(f"Directory '{db_folder}' created or already exists.")
    except OSError as e:
        print(f"Error: Failed to create directory '{db_folder}': {e}", file=sys.stderr)
        sys.exit(1)

    # --- Step 3: Execute makeblastdb command ---
    print("\n[2/3] Constructing BLAST database...")
    makeblastdb_cmd = [
        'makeblastdb',
        '-in', target_file,
        '-dbtype', 'nucl',
        '-parse_seqids',
        '-out', db_path_prefix  # <--- Use path including the folder
    ]
    print(f"Executing command: {' '.join(makeblastdb_cmd)}")
    
    try:
        subprocess.run(makeblastdb_cmd, check=True, capture_output=True, text=True)
        print("Database construction successful!")
    except FileNotFoundError:
        print("\nError: 'makeblastdb' command not found.", file=sys.stderr)
        print("Please ensure NCBI BLAST+ is installed and its path is added to your system's PATH environment variable.", file=sys.stderr)
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"\nError: 'makeblastdb' execution failed. Return code: {e.returncode}", file=sys.stderr)
        print(f"Error message:\n{e.stderr}", file=sys.stderr)
        sys.exit(1)

    # --- Step 4: Execute blastn command ---
    print("\n[3/3] Running blastn alignment...")
    blastn_cmd = [
        'blastn',
        '-query', query_file,
        '-db', db_path_prefix,  
        '-out', output_file,
        '-outfmt', '6'
    ]
    print(f"Executing command: {' '.join(blastn_cmd)}")

    try:
        subprocess.run(blastn_cmd, check=True)
        print("blastn alignment successful!")
    except FileNotFoundError:
        print("\nError: 'blastn' command not found.", file=sys.stderr)
        print("Please ensure NCBI BLAST+ is installed and its path is added to your system's PATH environment variable.", file=sys.stderr)
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"\nError: 'blastn' execution failed. Return code: {e.returncode}", file=sys.stderr)
        print(f"Error message:\n{e.stderr}", file=sys.stderr)
        sys.exit(1)
    
    print("\n--- Pipeline Execution Completed ---")
    print(f"Alignment results saved to: {output_file}")


def main():
    """
    Main function to parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="A Python script to automate makeblastdb and blastn while organizing database files into a directory.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "-q", "--query", 
        required=True, 
        help="[Required] Query sequence file (FASTA format)."
    )
    parser.add_argument(
        "-t", "--target", 
        required=True, 
        help="[Required] Target/Reference sequence file (FASTA format) for database construction."
    )
    
    args = parser.parse_args()
    
    run_blast_pipeline(args.query, args.target)

if __name__ == "__main__":
    main()