import os
import glob
from collections import Counter

def read_fasta_to_list(filename):
    """
    Reads a FASTA format file and returns a list containing all sequences.
    """
    sequences = []
    try:
        with open(filename, 'r') as f:
            current_sequence = []
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if current_sequence:
                        sequences.append("".join(current_sequence))
                    current_sequence = []
                else:
                    current_sequence.append(line.upper()) 
            
            if current_sequence:
                sequences.append("".join(current_sequence))
                
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return []
        
    return sequences

def calculate_sharing_value(list1, list2):
    """
    Calculates the sharing value (sum of occurrences) between two sequence lists.
    """
    counts1 = Counter(list1)
    counts2 = Counter(list2)
    
    total_sharing_value = 0
    common_sequences_keys = counts1.keys() & counts2.keys()
    
    for seq in common_sequences_keys:
        total_sharing_value += counts1[seq] + counts2[seq]
        
    return total_sharing_value

def main():
    """
    Main function: automatically finds .fa/.fasta files, performs pairwise 
    comparisons, and outputs the results in a tabular format.
    """
    print("Searching for .fa and .fasta files in the current directory...")
    fasta_files = glob.glob('*.fa') + glob.glob('*.fasta')
    
    if not fasta_files:
        print("Error: No .fa or .fasta files found in the current directory.")
        return

    # Sort the file list to ensure consistent ordering of results across runs
    fasta_files.sort()
    
    print(f"Found {len(fasta_files)} FASTA files: {fasta_files}")
    
    # --- Optimization: Load all files into memory at once ---
    print("\nPreloading all file sequences into memory, please wait...")
    all_sequences_data = {}
    for f in fasta_files:
        print(f"  - Reading {f}...")
        all_sequences_data[f] = read_fasta_to_list(f)
    print("All files loaded successfully!")

    # Define output filename
    output_filename = "sharing_results_pairwise.tsv"
    
    with open(output_filename, 'w') as out_file:
        # Write header
        header = (
            "File1\tFile2\tTotalSeqs_File1\tTotalSeqs_File2\t"
            "SharedValue\tTotalSeqs_In_Pair\tSharingRatio\n"
        )
        out_file.write(header)
        
        print(f"\nStarting pairwise comparisons; results will be written to {output_filename} ...")
        
        num_files = len(fasta_files)
        # --- Core Logic Update ---
        # Use nested index loops to avoid self-comparison (i != j) and redundant comparisons (j > i).
        # For example, compare (A, B) but not (B, A).
        for i in range(num_files):
            for j in range(i + 1, num_files):
                file1 = fasta_files[i]
                file2 = fasta_files[j]
                
                # Retrieve sequences from preloaded data
                sequences1 = all_sequences_data[file1]
                sequences2 = all_sequences_data[file2]
                
                total_seqs1 = len(sequences1)
                total_seqs2 = len(sequences2)
                
                # Calculate sharing value
                shared_value = calculate_sharing_value(sequences1, sequences2)
                
                # Calculate total sequences and ratio
                total_sequences_in_pair = total_seqs1 + total_seqs2
                ratio = (shared_value / total_sequences_in_pair) if total_sequences_in_pair > 0 else 0
                
                # Prepare data row for writing
                result_line = (
                    f"{file1}\t{file2}\t{total_seqs1}\t{total_seqs2}\t"
                    f"{shared_value}\t{total_sequences_in_pair}\t{ratio:.4f}\n"
                )
                
                # Write to file
                out_file.write(result_line)

    print(f"\nProcessing complete! All results have been saved to the tabular file {output_filename}.")


if __name__ == "__main__":
    main()