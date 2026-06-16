#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import argparse
from collections import Counter

def parse_fasta(file_handle):
    """
    A simple FASTA file parser.
    Reads sequences one by one from an open file handle.
    Yields sequences as (header, sequence) tuples.
    Correctly handles multi-line sequences.
    """
    header = None
    sequence_parts = []
    for line in file_handle:
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if header is not None:
                yield header, ''.join(sequence_parts)
            header = line
            sequence_parts = []
        else:
            sequence_parts.append(line.upper()) # Normalize sequence to uppercase to ignore case sensitivity
    
    if header is not None:
        yield header, ''.join(sequence_parts)

def analyze_sequences(directory, output_file):
    """
    Analyzes sequences in all .fa files within the directory and generates a statistical report.
    """
    fasta_files = sorted(glob.glob(os.path.join(directory, '*.fa')) + glob.glob(os.path.join(directory, '*.fasta')))
    
    if not fasta_files:
        print(f"Error: No .fa or .fasta files found in directory '{directory}'.")
        return

    print(f"[*] Found {len(fasta_files)} FASTA files for analysis.")

    # --- Step 1: Data Collection ---
    # all_sequences_by_file: Stores the sequence counts for each file.
    # Structure: {'file1.fa': Counter({'ACGT': 2, 'GCTA': 1}), ...}
    all_sequences_by_file = {}

    for filepath in fasta_files:
        filename = os.path.basename(filepath)
        print(f"    - Collecting data from: {filename}")
        
        try:
            with open(filepath, 'r') as f:
                # Use Counter to efficiently count occurrences of each sequence in the file
                file_sequence_counter = Counter(seq for _, seq in parse_fasta(f))
                all_sequences_by_file[filename] = file_sequence_counter
        except Exception as e:
            print(f"Warning: Error processing file {filename}, skipped. Error: {e}")
            continue

    if not all_sequences_by_file:
        print("Warning: No valid sequences found in any of the FASTA files.")
        open(output_file, 'w').close()
        return

    print("\n[*] Data collection complete. Calculating and generating report...")

    # --- Step 2: Calculation and Output ---
    try:
        with open(output_file, 'w') as out_f:
            out_f.write("sequence_name\tcopies_in_file\tshared_across_files\n")
            
            # Re-iterate through files to generate output in the specified order
            for filepath in fasta_files:
                filename = os.path.basename(filepath)
                if filename not in all_sequences_by_file:
                    continue # Skip files that failed during Step 1

                # **Key mechanism**: Use a set to track unique sequence contents already processed for this file
                processed_sequences_in_this_file = set()

                print(f"    - Generating report for: {filename}")
                with open(filepath, 'r') as f:
                    for header, seq in parse_fasta(f):
                        # Skip if this specific sequence content has already been recorded for this file
                        if seq in processed_sequences_in_this_file:
                            continue
                        
                        # First encounter of this sequence content in the current file; process and output
                        # 1. Sequence name (extracted from the first occurring entry)
                        sequence_name = header.lstrip('>')
                        
                        # 2. Total copy number of this sequence in the current file (retrieved from pre-computed counts)
                        copies_in_this_file = all_sequences_by_file[filename][seq]
                        
                        # 3. Number of files in which this sequence is present
                        shared_file_count = 0
                        for file_counter in all_sequences_by_file.values():
                            if seq in file_counter:
                                shared_file_count += 1
                        
                        # Write result row
                        out_f.write(f"{sequence_name}\t{copies_in_this_file}\t{shared_file_count}\n")
                        
                        # **Mark**: Add sequence content to the set to prevent duplicate output for subsequent identical sequences
                        processed_sequences_in_this_file.add(seq)

    except IOError as e:
        print(f"Error: Unable to write to output file '{output_file}'. Error: {e}")
        return

    print(f"\n[+] Processing complete! Results saved to: {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generates statistical reports for the first occurrence of unique sequences in each FASTA file.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument(
        '-o', '--output',
        default='sequence_first_occurrence_stats.tsv',
        help='Name of the output TSV file (default: sequence_first_occurrence_stats.tsv)'
    )
    
    args = parser.parse_args()
    
    current_directory = '.'
    analyze_sequences(current_directory, args.output)