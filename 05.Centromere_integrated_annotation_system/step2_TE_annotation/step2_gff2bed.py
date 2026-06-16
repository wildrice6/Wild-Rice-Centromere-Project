#!/usr/bin/env python3

import argparse
import sys

def parse_attributes(attr_string):
    """
    Parses the GFF attribute string (9th column) to find the 'id' value.

    Args:
        attr_string (str): The semicolon-separated attribute string.

    Returns:
        str: The value associated with the 'id' key, or '.' if not found.
    """
    attributes = attr_string.split(';')
    for attr in attributes:
        # Ensure there is an '=' to split on
        if '=' in attr:
            key, value = attr.split('=', 1)
            if key.lower() == 'id':
                return value
    return '.' # Return a placeholder if no ID is found

def convert_gff_to_bed(input_file, output_file):
    """
    Converts a GFF3 file to a BED file with custom column mapping.

    BED format will be:
    1. chrom
    2. start (0-based)
    3. end
    4. name (from GFF id=)
    5. score (from GFF type/column 3)
    6. strand
    7. extra (the full GFF attribute string)
    """
    print(f"Starting conversion from {input_file} to {output_file}...")
    
    try:
        with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
            for line in f_in:
                # Skip header/comment lines
                if line.startswith('#'):
                    continue
                
                # Skip empty lines
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split('\t')
                
                # A valid GFF line must have 9 columns
                if len(parts) != 9:
                    print(f"Warning: Skipping malformed line with {len(parts)} columns: {line}", file=sys.stderr)
                    continue
                
                # Assign GFF columns to variables for clarity
                chrom = parts[0]
                gff_type = parts[2]
                gff_start = int(parts[3])
                gff_end = int(parts[4])
                strand = parts[6]
                attributes = parts[8]
                
                # --- Coordinate Conversion ---
                # GFF is 1-based, inclusive. BED is 0-based, half-open.
                # bed_start = gff_start - 1
                # bed_end = gff_end
                bed_start = gff_start - 1
                bed_end = gff_end
                
                # --- Column Mapping as per user request ---
                # BED col 4: name (from GFF 'id=')
                bed_name = parse_attributes(attributes)
                
                # BED col 5: score (from GFF col 3 'type')
                bed_score = gff_type
                
                # BED col 7: extra info (the full GFF attribute string)
                extra_info = attributes
                
                # Assemble the output BED line
                bed_line_parts = [
                    str(chrom),
                    str(bed_start),
                    str(bed_end),
                    str(bed_name),
                    str(bed_score),
                    str(strand),
                    str(extra_info)
                ]
                
                # Write to output file
                f_out.write("\t".join(bed_line_parts) + "\n")
                
    except FileNotFoundError:
        print(f"Error: Input file not found at {input_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)

    print("Conversion completed successfully.")


def main():
    """Main function to parse arguments and run the script."""
    parser = argparse.ArgumentParser(
        description="Convert a GFF3 file to a 7-column BED file with specific column mapping.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "--input", 
        required=True, 
        help="Path to the input GFF3 file."
    )
    parser.add_argument(
        "--output", 
        required=True, 
        help="Path for the output BED file."
    )
    
    args = parser.parse_args()
    
    convert_gff_to_bed(args.input, args.output)

if __name__ == "__main__":
    main()
