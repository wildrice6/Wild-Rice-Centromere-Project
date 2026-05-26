import argparse
import sys

def blast_to_bed(blast_file, bed_file):
    """
    Converts BLAST outfmt 6 result files into BED format files.

    Args:
        blast_file (str): Path to the input BLAST (outfmt 6) file.
        bed_file (str): Path to the output BED file.
    """
    print(f"Reading BLAST file: {blast_file}")
    
    lines_processed = 0
    lines_skipped = 0
    
    try:
        with open(blast_file, 'r') as infile, open(bed_file, 'w') as outfile:
            for line in infile:
                # Skip empty lines or comment lines
                if not line.strip() or line.startswith('#'):
                    continue
                
                try:
                    # Split the line; outfmt 6 is tab or space delimited
                    fields = line.strip().split()
                    
                    # Ensure the line contains at least 12 columns
                    if len(fields) < 12:
                        lines_skipped += 1
                        continue

                    # --- 1. Extract required information from BLAST results ---
                    # qseqid (Query sequence ID) -> BED name
                    qseqid = fields[0]
                    # sseqid (Subject sequence ID) -> BED chrom
                    sseqid = fields[1]
                    # sstart, send (Alignment start/stop on subject sequence, 1-based)
                    sstart = int(fields[8])
                    send = int(fields[9])
                    # bitscore -> BED score
                    bitscore = fields[11]

                    # --- 2. Convert to BED format fields ---
                    # BED chrom: Subject sequence ID
                    bed_chrom = sseqid
                    
                    # BED name: Query sequence ID
                    bed_name = qseqid
                    
                    # BED score: bitscore
                    bed_score = bed_score_val = bitscore
                    
                    # BED strand:
                    # If sstart < send, it indicates the forward strand '+'
                    # If sstart > send, it indicates the reverse strand '-'
                    if sstart < send:
                        bed_strand = '+'
                        # BED chromStart is 0-based, so subtract 1 from sstart
                        bed_start = sstart - 1
                        bed_end = send
                    else: # sstart > send
                        bed_strand = '-'
                        # BED start must always be less than end, even on the reverse strand
                        bed_start = send - 1
                        bed_end = sstart
                    
                    # --- 3. Format and write to BED file ---
                    # BED format is tab-delimited
                    bed_line = f"{bed_chrom}\t{bed_start}\t{bed_end}\t{bed_name}\t{bed_score}\t{bed_strand}\n"
                    outfile.write(bed_line)
                    
                    lines_processed += 1
                    
                except (ValueError, IndexError) as e:
                    # Skip malformed lines (e.g., non-numeric coordinates) and print a warning
                    print(f"Warning: Skipping malformed line: '{line.strip()}' -> Error: {e}", file=sys.stderr)
                    lines_skipped += 1
                    continue
                    
    except FileNotFoundError:
        print(f"Error: Input file '{blast_file}' not found.", file=sys.stderr)
        sys.exit(1)
        
    print("\nConversion complete!")
    print(f"Successfully processed {lines_processed} lines.")
    if lines_skipped > 0:
        print(f"Skipped {lines_skipped} malformed lines.")
    print(f"BED file saved to: {bed_file}")

def main():
    """
    Main function to parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Converts BLASTn outfmt 6 result files into standard 6-column BED files.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '-i', "--blast_file", 
        help="Input BLAST result file (outfmt 6 format).\nExample: blast.txt"
    )
    parser.add_argument(
        '-o', "--bed_file", 
        help="Output BED filename.\nExample: output.bed"
    )
    
    args = parser.parse_args()
    
    blast_to_bed(args.blast_file, args.bed_file)

if __name__ == "__main__":
    main()