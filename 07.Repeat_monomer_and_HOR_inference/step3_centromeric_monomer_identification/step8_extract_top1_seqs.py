#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import subprocess
import sys
from shutil import which

def check_tool_exists(name):
    """Check if the required command-line tool exists in the system PATH."""
    if which(name) is None:
        print(f"Error: Required tool '{name}' not found or not in your system PATH.")
        print("Please install samtools and ensure it is available in your PATH.")
        sys.exit(1)

def parse_pos_filename(filepath):
    """
    Parse the accession name, chromosome, and window information from the coordinate filename.
    Example: AA_Ogla_hap1.Chr01.w100_top20_regions_list.txt -> ('AA_Ogla_hap1', 'Chr01', 'w100')
    """
    filename = os.path.basename(filepath)
    parts = filename.split('.')
    if len(parts) < 3:
        print(f"Warning: Filename format does not meet expectations, skipping file: {filename}")
        return None, None, None
        
    material_name = parts[0]
    chromosome = parts[1]
    # Extract 'w100' from 'w100_top20_regions_list'
    window_info = parts[2].split('_')[0]
    
    return material_name, chromosome, window_info

def extract_top1_regions(pos_file):
    """
    Read the coordinate file and return all regions where the second column is 'Top1'.
    """
    top1_regions = []
    try:
        with open(pos_file, 'r') as f:
            # Skip header
            next(f)
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                columns = line.split() # Split by whitespace or tab
                if len(columns) >= 2 and columns[1] == 'Top1':
                    region = columns[0]
                    top1_regions.append(region)
    except FileNotFoundError:
        print(f"Error: Coordinate file not found: {pos_file}")
        return []
    except Exception as e:
        print(f"Error reading coordinate file {pos_file}: {e}")
        return []
        
    return top1_regions

def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description="Extract 'Top1' region sequences from reference genomes based on a specified list of coordinate files.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('--pos_list', required=True, 
                        help="Required: List file containing paths to coordinate files.\n"
                             "Each line should be an absolute or relative path to a coordinate file.")
    parser.add_argument('--ref_dir', required=True, 
                        help="Required: Directory containing reference genomes (.fasta and .fasta.fai) for all materials.")
    parser.add_argument('--out_dir', required=True, 
                        help="Required: Directory to save the output FASTA files.")
    
    args = parser.parse_args()

    # 1. Verify dependencies
    check_tool_exists('samtools')

    # 2. Check and create output directory
    os.makedirs(args.out_dir, exist_ok=True)

    # 3. Read the list of coordinate files
    try:
        with open(args.pos_list, 'r') as f_list:
            pos_files = [line.strip() for line in f_list if line.strip()]
    except FileNotFoundError:
        print(f"Error: Coordinate list file not found: {args.pos_list}")
        sys.exit(1)

    print(f"Total coordinate files found for processing: {len(pos_files)}")

    # 4. Iterate and process each coordinate file
    for pos_file_path in pos_files:
        print(f"\n--- Processing: {os.path.basename(pos_file_path)} ---")
        
        # Parse metadata from filename
        material, chrom, window = parse_pos_filename(pos_file_path)
        if not all([material, chrom, window]):
            continue

        # Identify Top1 regions for extraction
        regions_to_extract = extract_top1_regions(pos_file_path)
        if not regions_to_extract:
            print(f"No 'Top1' regions found in {os.path.basename(pos_file_path)}, skipping.")
            continue
        
        print(f"Found {len(regions_to_extract)} 'Top1' regions.")

        # Construct paths for reference genome and output file
        ref_fasta_path = os.path.join(args.ref_dir, f"{material}.fasta")
        output_fasta_path = os.path.join(args.out_dir, f"{material}.{chrom}.{window}.Top1.fa")

        # Verify reference genome existence
        if not os.path.exists(ref_fasta_path):
            print(f"Warning: Corresponding reference genome not found: {ref_fasta_path}, skipping.")
            continue
        
        # Verify FASTA index file existence
        if not os.path.exists(ref_fasta_path + '.fai'):
            print(f"Warning: Reference genome index file not found: {ref_fasta_path}.fai")
            print("Please create an index using 'samtools faidx' before proceeding. Skipping.")
            continue

        # Construct samtools command
        # Passing all regions in a single call for higher efficiency
        command = ['samtools', 'faidx', ref_fasta_path] + regions_to_extract
        
        print(f"Executing sequence extraction...")
        try:
            # Execute command and capture output
            result = subprocess.run(
                command, 
                capture_output=True, 
                text=True, 
                check=True  # Raise exception if samtools returns an error code
            )
            
            # Write extracted sequences to the output file
            with open(output_fasta_path, 'w') as f_out:
                f_out.write(result.stdout)
            
            print(f"Successfully extracted sequences to: {output_fasta_path}")

        except subprocess.CalledProcessError as e:
            print(f"Error: 'samtools faidx' command failed.")
            print(f"Return code: {e.returncode}")
            print(f"Standard Error output:\n{e.stderr}")
        except Exception as e:
            print(f"An unexpected error occurred while executing samtools: {e}")

    print("\nAll files processed successfully.")

if __name__ == "__main__":
    main()