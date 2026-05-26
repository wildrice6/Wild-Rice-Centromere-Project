#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys

def get_chrom_sort_key(chrom_name):
    """
    Generates a sort key for chromosome names to achieve natural sorting 
    (e.g., Chr1, Chr2, ..., Chr10, ChrX). 
    Numerical and non-numerical chromosomes are sorted separately.
    """
    name = chrom_name.lower().replace('chr', '')
    if name.isdigit():
        return (0, int(name))  # (Group, numerical value) -> Numerical chromosomes sorted by magnitude
    else:
        return (1, name)      # (Group, string value) -> Alphanumeric chromosomes sorted alphabetically

def process_bed_file(filepath):
    """
    Reads a BED file and calculates the total length for each chromosome.

    Args:
        filepath (str): Path to the BED file.

    Returns:
        dict: A dictionary where keys are chromosome names and values are the total lengths.
    """
    lengths = {}
    print(f"[*] Processing file: {filepath}")
    try:
        with open(filepath, 'r') as f:
            for i, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split()
                if len(parts) < 3:
                    print(f"Warning: {filepath} line {i} is malformed and will be skipped: '{line}'", file=sys.stderr)
                    continue

                chrom = parts[0]
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                    length = end - start
                    # Accumulate length
                    lengths[chrom] = lengths.get(chrom, 0) + length
                except ValueError:
                    print(f"Warning: {filepath} line {i} coordinates are not valid integers and will be skipped: '{line}'", file=sys.stderr)
                    continue
    except FileNotFoundError:
        print(f"Error: File not found: {filepath}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred while processing file {filepath}: {e}", file=sys.stderr)
        sys.exit(1)
        
    print(f"[+] File {filepath} processing complete.")
    return lengths

def main():
    # 1. Setup command-line argument parser
    parser = argparse.ArgumentParser(
        description="Calculates the coverage of elements within CEN regions, providing per-chromosome and aggregate statistics.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    # 2. Add required arguments
    parser.add_argument('--element', required=True, help='Path to the element BED file (chr, start, end, ...)')
    parser.add_argument('--cen', required=True, help='Path to the CEN/target region BED file (chr, start, end, ...)')
    parser.add_argument('--output_tsv', required=True, help='Path for the output TSV file containing per-chromosome statistics')
    parser.add_argument('--output_summary', required=True, help='Path for the output TXT file containing aggregate summary statistics')
    
    args = parser.parse_args()

    # 3. Process input files
    cen_lengths = process_bed_file(args.cen)
    element_lengths = process_bed_file(args.element)

    if not cen_lengths:
        print("Error: No valid data found in the CEN BED file; terminating process.", file=sys.stderr)
        sys.exit(1)

    # 4. Generate per-chromosome TSV file
    print(f"[*] Generating TSV output file: {args.output_tsv}")
    try:
        with open(args.output_tsv, 'w') as tsv_file:
            # Write header
            tsv_file.write("chr\tcen_len\telement_len\tpercent\n")
            
            # Sort by natural chromosome order
            sorted_chroms = sorted(cen_lengths.keys(), key=get_chrom_sort_key)
            
            for chrom in sorted_chroms:
                cen_len = cen_lengths[chrom]
                # If the chromosome is absent in the element file, length is set to 0
                element_len = element_lengths.get(chrom, 0)
                
                # Calculate percentage; handle cases where cen_len is zero to avoid division by zero errors
                if cen_len > 0:
                    percentage = (element_len / cen_len) * 100
                else:
                    percentage = 0.0
                
                # Write row; percentage formatted to four decimal places
                tsv_file.write(f"{chrom}\t{cen_len}\t{element_len}\t{percentage:.4f}\n")
    except IOError as e:
        print(f"Error: Unable to write TSV file: {e}", file=sys.stderr)
        sys.exit(1)
    
    print(f"[+] TSV file {args.output_tsv} generated successfully.")

    # 5. Calculate aggregate statistics and generate summary file
    print(f"[*] Generating summary output file: {args.output_summary}")
    
    total_cen_len = sum(cen_lengths.values())
    # Note: Total element length is calculated across all records, not restricted to chromosomes present in the CEN file
    total_element_len = sum(element_lengths.values())
    
    if total_cen_len > 0:
        overall_percentage = (total_element_len / total_cen_len) * 100
    else:
        overall_percentage = 0.0
        
    try:
        with open(args.output_summary, 'w') as summary_file:
            summary_file.write("====================\n")
            summary_file.write(f"Source: {args.element}\n")
            summary_file.write("Overall Statistics\n")
            summary_file.write("====================\n")
            summary_file.write(f"Total CEN Length: {total_cen_len}\n")
            summary_file.write(f"Total Element Length: {total_element_len}\n")
            summary_file.write(f"Overall Percentage (Element / CEN): {overall_percentage:.4f}%\n")
            summary_file.write('-' * 100)
            summary_file.write(f"\n")
    except IOError as e:
        print(f"Error: Unable to write summary file: {e}", file=sys.stderr)
        sys.exit(1)
        
    print(f"[+] Summary file {args.output_summary} generated successfully.")
    print("\n[+] All tasks completed successfully!")


if __name__ == "__main__":
    main()