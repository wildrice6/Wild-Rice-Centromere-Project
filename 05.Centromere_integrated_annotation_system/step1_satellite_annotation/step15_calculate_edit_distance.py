#!/usr/bin/env python3

import argparse
import sys
from itertools import combinations

# Attempt to import the Levenshtein library and provide a prompt if it fails
try:
    import Levenshtein
except ImportError:
    print("Error: The 'Levenshtein' library is missing.", file=sys.stderr)
    print("Please install it using the command: 'pip install python-Levenshtein'.", file=sys.stderr)
    sys.exit(1)

def parse_fasta(file_path):
    """
    Parses a FASTA file and returns a dictionary mapping sequence names to sequences.
    
    Args:
        file_path (str): Path to the FASTA file.
        
    Returns:
        dict: A dictionary where keys are sequence headers and values are uppercase DNA sequences.
    """
    sequences = {}
    current_header = None
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    current_header = line[1:]
                    sequences[current_header] = []
                elif current_header:
                    sequences[current_header].append(line)

        for header, seq_list in sequences.items():
            sequences[header] = "".join(seq_list).upper()
            
        return sequences
    except FileNotFoundError:
        print(f"Error: Input file '{file_path}' not found.", file=sys.stderr)
        sys.exit(1)

def calculate_circular_min_distance(seq1, seq2):
    """
    Calculates the minimum circular edit distance between two sequences.
    
    Performs all possible cyclic shifts on seq2, calculating the edit distance 
    against a fixed seq1 for each shift, and returns the minimum value found. 
    This operation is performed even if the sequences have different lengths.
    
    Args:
        seq1 (str): The first sequence (fixed reference sequence).
        seq2 (str): The second sequence (the sequence to be cyclically shifted).
        
    Returns:
        int: The minimum circular edit distance.
    """
    min_distance = float('inf')
    temp_seq2 = seq2
    
    if not temp_seq2:
        return Levenshtein.distance(seq1, seq2)

    for _ in range(len(seq2)):
        dist = Levenshtein.distance(seq1, temp_seq2)
        if dist < min_distance:
            min_distance = dist
        if min_distance == 0:
            break
        temp_seq2 = temp_seq2[1:] + temp_seq2[0]
        
    return min_distance

def main():
    """Main execution function"""
    parser = argparse.ArgumentParser(
        description="Calculates the minimum circular edit distance between sequence pairs in a FASTA file (supports sequences of different lengths and utilizes streaming output to optimize memory)."
    )
    parser.add_argument(
        "--input", 
        required=True, 
        help="Path to the input FASTA file."
    )
    parser.add_argument(
        "--output", 
        required=True, 
        help="Path to the output TSV results file."
    )
    
    args = parser.parse_args()
    
    # 1. Parse the FASTA file
    print(f" Reading and parsing file: {args.input}")
    sequences = parse_fasta(args.input)
    if not sequences or len(sequences) < 2:
        print("Error: The input file must contain at least two valid sequences for comparison.", file=sys.stderr)
        sys.exit(1)
    print(f" Successfully loaded {len(sequences)} sequences.")
    
    print(" Calculating and streaming results to: {}...".format(args.output))
    
    try:
        # 2. Open output file before starting the calculation loop
        with open(args.output, 'w') as f_out:
            # 3. Write the header first
            f_out.write("Sequence_A\tSequence_B\tMin_Edit_Distance\n")
            
            sequence_names = list(sequences.keys())
            pair_iterator = combinations(sequence_names, 2)
            
            # Use a counter to estimate the total number of pairs and display progress
            # For very large N, calculating the exact length of combinations might be slow or memory-intensive
            # Thus, a simple arithmetic estimation is used
            num_seqs = len(sequence_names)
            total_pairs = num_seqs * (num_seqs - 1) // 2
            
            # 4. Enter the calculation loop
            for i, (name1, name2) in enumerate(pair_iterator):
                seq1 = sequences[name1]
                seq2 = sequences[name2]
                
                # Core calculation logic
                min_dist = calculate_circular_min_distance(seq1, seq2)
                
                # 5. Stream each result to the file immediately after calculation
                f_out.write(f"{name1}\t{name2}\t{min_dist}\n")
                
                # Print progress updates
                if (i + 1) % 100 == 0 or (i + 1) == total_pairs:
                    progress = (i + 1) / total_pairs * 100
                    print(f"  ...Processed {i+1}/{total_pairs} pairs ({progress:.1f}%)", end='\r')

    except IOError as e:
        print(f"\nError: Unable to write to output file '{args.output}': {e}", file=sys.stderr)
        sys.exit(1)
        
    print("\n\n All tasks completed successfully!")

if __name__ == "__main__":
    main()