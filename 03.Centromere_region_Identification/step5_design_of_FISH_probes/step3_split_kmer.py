#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import re

def sanitize_filename(header):
    """
    Converts FASTA header information into a safe filename.
    Removes or replaces characters that are illegal in most file systems.
    
    Args:
        header (str): The FASTA sequence header information (excluding '>').
        
    Returns:
        str: A sanitized string suitable for use as a filename.
    """
    # Remove illegal characters: / \ ? * " < > |
    # Replace spaces with underscores to enhance readability and compatibility
    filename = re.sub(r'[\\/*?:"<>|]', "_", header)
    filename = filename.replace(' ', '_')
    return filename

def split_fasta(input_file, output_dir):
    """
    Splits a multi-sequence FASTA file into multiple single-sequence FASTA files.
    
    Args:
        input_file (str): Path to the input multi-sequence FASTA file.
        output_dir (str): Path to the directory where output files will be stored.
    """
    # 1. Ensure the output directory exists; create it if it does not
    if not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
            print(f"Created output directory: {output_dir}")
        except OSError as e:
            print(f"Error: Failed to create output directory '{output_dir}'. Reason: {e}")
            return

    try:
        with open(input_file, 'r') as f_in:
            output_file_handle = None
            file_count = 0
            
            for line in f_in:
                # Remove leading and trailing whitespace
                line = line.strip()
                
                # Skip empty lines
                if not line:
                    continue
                
                # Identify header lines
                if line.startswith('>'):
                    # Close the current file handle if one is already open
                    if output_file_handle:
                        output_file_handle.close()
                    
                    # Generate a filename based on the header information
                    header = line[1:]  # Remove the leading '>'
                    filename = sanitize_filename(header) + '.fasta'
                    output_filepath = os.path.join(output_dir, filename)
                    
                    # Open a new file for writing the current sequence record
                    output_file_handle = open(output_filepath, 'w')
                    output_file_handle.write(line + '\n')
                    file_count += 1
                    print(f"Creating file: {filename}")
                
                # Write sequence lines to the active file handle
                elif output_file_handle:
                    # Ensure a file handle is open (i.e., a header has been processed)
                    output_file_handle.write(line + '\n')

            # Close the final file handle
            if output_file_handle:
                output_file_handle.close()
                
            print(f"\nProcessing complete. Total files generated: {file_count}")

    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


def main():
    """
    Main function to parse command-line arguments and invoke the FASTA splitting function.
    """
    parser = argparse.ArgumentParser(
        description="Split a multi-sequence FASTA file into multiple single-sequence FASTA files and save them to a specified directory."
    )
    parser.add_argument(
        '--input', 
        required=True, 
        help='Path to the input FASTA file.'
    )
    parser.add_argument(
        '--output', 
        required=True, 
        help='Path to the output directory for storing split files.'
    )
    
    args = parser.parse_args()
    
    split_fasta(args.input, args.output)


if __name__ == "__main__":
    main()