import argparse
import os
import sys
from itertools import combinations
import csv

# Define the length of k-mers
KMER_SIZE = 15

def generate_k_mers(sequence, k):
    """
    Generate the set of all k-mers for a given sequence.
    """
    k_mers = set()
    n = len(sequence)
    if n < k:
        return k_mers
    
    for i in range(n - k + 1):
        k_mers.add(sequence[i:i + k])
    return k_mers

def calculate_jaccard(set_a, set_b):
    """
    Calculate the Jaccard similarity between two sets.
    J(A, B) = |A ∩ B| / |A ∪ B|
    """
    if not set_a and not set_b:
        return 0.0
    
    # Calculate sizes of intersection and union
    intersection = len(set_a.intersection(set_b))
    union = len(set_a.union(set_b))
    
    if union == 0:
        return 0.0
        
    return intersection / union

def parse_fasta(input_path):
    """
    Read a FASTA file and generate a dictionary mapping sequence names to k-mer sets.
    """
    seq_data = {}
    current_name = None
    current_sequence = []

    print(f"Reading file and generating {KMER_SIZE}-mer sets...")
    
    try:
        with open(input_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                if line.startswith('>'):
                    # Handle a new sequence
                    if current_name is not None:
                        full_seq = "".join(current_sequence)
                        seq_data[current_name] = generate_k_mers(full_seq, KMER_SIZE)
                    
                    # Extract the sequence name, typically the first word after '>'
                    current_name = line[1:].split()[0]
                    current_sequence = []
                else:
                    # Collect sequence data
                    current_sequence.append(line)

            # Handle the last sequence in the file
            if current_name is not None:
                full_seq = "".join(current_sequence)
                seq_data[current_name] = generate_k_mers(full_seq, KMER_SIZE)

    except FileNotFoundError:
        print(f"Error: Input file {input_path} not found.")
        sys.exit(1)
    
    print(f"Successfully parsed {len(seq_data)} sequences.")
    return seq_data

def main():
    parser = argparse.ArgumentParser(
        description="Calculate k-mer Jaccard similarity between all sequence pairs in a FASTA file."
    )
    parser.add_argument('--input', required=True, 
                        help='Path to the input FASTA file')
    parser.add_argument('--output', required=True, 
                        help='Path to the output TSV file')
    
    args = parser.parse_args()

    # 1. Parse FASTA and generate k-mer sets
    kmer_sets = parse_fasta(args.input)
    
    if len(kmer_sets) < 2:
        print("Insufficient data: at least two sequences are required for comparison.")
        sys.exit(0)

    # 2. Perform iterative comparisons and output results
    
    # Prepare the list of sequence names
    names = list(kmer_sets.keys())

    print(f"Starting calculation of Jaccard similarity for {len(names) * (len(names) - 1) // 2} sequence pairs...")

    try:
        # Use csv.writer to ensure correct output format (TSV: delimiter='\t')
        with open(args.output, 'w', newline='') as outfile:
            writer = csv.writer(outfile, delimiter='\t')
            # Write header row
            writer.writerow(['Sequence_A', 'Sequence_B', 'Jaccard_Similarity'])
            
            # Use combinations to iterate over all unique sequence pairs (A, B), 
            # avoiding redundant (B, A) pairs and self-comparison (A, A)
            total_pairs = len(names) * (len(names) - 1) // 2
            counter = 0

            for name_a, name_b in combinations(names, 2):
                set_a = kmer_sets[name_a]
                set_b = kmer_sets[name_b]
                
                similarity = calculate_jaccard(set_a, set_b)
                
                # Write results to the output file during calculation
                writer.writerow([name_a, name_b, f"{similarity:.6f}"])
                
                counter += 1
                if counter % 1000 == 0:
                    print(f"  Processed {counter}/{total_pairs} pairs...", end='\r', flush=True)

        print(f"\nCalculation complete. Results saved to {args.output}")

    except Exception as e:
        print(f"An error occurred while writing the output file: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()