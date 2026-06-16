#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from collections import Counter

def read_fasta_sequence(filepath):
    """
    Reads a sequence from a FASTA file containing a single record.
    The function ignores lines starting with '>' and concatenates all other 
    lines into a single string.
    
    Args:
        filepath (str): Path to the input FASTA file.
        
    Returns:
        str: The concatenated DNA sequence string; returns None if the file 
             is empty or the sequence cannot be found.
    """
    sequence_parts = []
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    # Remove trailing newline characters and whitespace
                    sequence_parts.append(line.strip())
    except FileNotFoundError:
        print(f"Error: Input file '{filepath}' not found.")
        return None
        
    if not sequence_parts:
        return None
        
    return "".join(sequence_parts)

def generate_and_count_kmers(sequence, k_size):
    """
    Generates k-mers from a given sequence and calculates their frequencies.
    
    Args:
        sequence (str): DNA sequence.
        k_size (int): The length of the k-mer.
        
    Returns:
        collections.Counter: A Counter object containing k-mers and their counts.
    """
    if not sequence or len(sequence) < k_size:
        return Counter()
        
    # Efficiently generate all k-mers using list comprehension
    kmers = [sequence[i:i+k_size] for i in range(len(sequence) - k_size + 1)]
    
    # Perform high-performance counting using collections.Counter
    return Counter(kmers)

def main():
    """
    Main function to parse command-line arguments and execute k-mer analysis.
    """
    parser = argparse.ArgumentParser(
        description="Extract sequences from a FASTA file, calculate 25-mer frequencies, and output the top 100 k-mers."
    )
    parser.add_argument(
        '--input', 
        required=True, 
        help='Path to the input FASTA file.'
    )
    parser.add_argument(
        '--output', 
        required=True, 
        help='Path to the output results file.'
    )
    
    args = parser.parse_args()
    
    KMER_SIZE = 25
    TOP_N = 100
    
    # 1. Read FASTA sequence
    print(f"Reading sequence from '{args.input}'...")
    sequence = read_fasta_sequence(args.input)
    
    if sequence is None:
        print("Failed to retrieve sequence from the input file. Terminating script.")
        return

    print(f"Sequence successfully retrieved. Total length: {len(sequence)} bp.")
    
    # 2. Generate and calculate k-mer frequencies
    print(f"Generating and counting {KMER_SIZE}-mers...")
    kmer_counts = generate_and_count_kmers(sequence, KMER_SIZE)
    
    if not kmer_counts:
        print(f"Sequence length is insufficient to generate {KMER_SIZE}-mers. Terminating script.")
        return

    print(f"Identified {len(kmer_counts)} unique {KMER_SIZE}-mers.")
    
    # 3. Retrieve top N k-mers by frequency
    # Counter.most_common(n) returns a list of the n most frequent elements sorted in descending order
    top_kmers = kmer_counts.most_common(TOP_N)
    
    # 4. Write results to the output file
    print(f"Writing the top {len(top_kmers)} k-mers to '{args.output}'...")
    try:
        with open(args.output, 'w') as f_out:
            for rank, (kmer, count) in enumerate(top_kmers, 1):
                # Format output header: >probe-rank[Rank]-number[Count]
                header = f">probe-rank{rank}-number{count}"
                f_out.write(header + '\n')
                f_out.write(kmer + '\n')
    except IOError:
        print(f"Error: Unable to write to the output file '{args.output}'.")
        return

    print("Processing complete!")

if __name__ == "__main__":
    main()