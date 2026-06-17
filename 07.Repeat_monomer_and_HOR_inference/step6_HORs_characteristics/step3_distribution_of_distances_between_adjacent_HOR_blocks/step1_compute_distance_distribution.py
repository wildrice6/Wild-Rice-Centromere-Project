import argparse
from collections import Counter
import sys
import os

def get_material_name(filepath):
    """
    Extracts the material name from the file path.
    Rule: The first substring of the filename separated by a period.
    Example: Extracts 'SampleA' from 'SampleA.data.txt'.
    """
    # 1. Retrieve the filename excluding the directory path
    filename = os.path.basename(filepath)
    # 2. Split by the period '.' and select the first component
    material_name = filename.split('.')[0]
    return material_name

def main():
    """
    Main function for parameter parsing, frequency statistics, and exporting a multi-column TSV file.
    """
    # 1. Configure command-line argument parsing
    parser = argparse.ArgumentParser(
        description="Calculates numerical frequencies from a single-column file, performs conversions, and generates formatted output."
    )
    parser.add_argument("--input", 
                        required=True, 
                        help="Path to the input single-column numerical file.")
    parser.add_argument("--output", 
                        required=True, 
                        help="Path to the output TSV file.")
    
    args = parser.parse_args()

    # 2. Extract material name
    material_name = get_material_name(args.input)
    print(f"Material name extracted from filename '{os.path.basename(args.input)}': {material_name}")

    # 3. Read file and calculate frequency statistics
    try:
        numbers = []
        print(f"Reading input file: {args.input}...")
        
        with open(args.input, 'r') as f_in:
            for line in f_in:
                clean_line = line.strip()
                if not clean_line:
                    continue
                try:
                    numbers.append(int(clean_line))
                except ValueError:
                    print(f"Warning: Non-integer value '{clean_line}' detected and skipped.", file=sys.stderr)

        if not numbers:
            print("Warning: Input file is empty or contains no valid numbers. Generating an empty output file.")
            # Create an empty output file and exit
            open(args.output, 'w').close()
            sys.exit(0)

        frequency_counts = Counter(numbers)
        total_count = len(numbers) # Calculate total population size
        
        print(f"Statistics complete: {len(frequency_counts)} unique numbers identified, total count: {total_count}.")

    except FileNotFoundError:
        print(f"Error: Input file does not exist -> {args.input}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unknown error occurred during file reading: {e}", file=sys.stderr)
        sys.exit(1)

    # 4. Sorting, computation, and writing to output file
    try:
        print(f"Writing results to output file: {args.output}...")
        
        with open(args.output, 'w') as f_out:
            # Write header
            header = "material_name\tnumber_div_1M\tcount\tpercentage\n"
            f_out.write(header)
            
            # Iterate through results sorted by numerical magnitude
            for number, count in sorted(frequency_counts.items()):
                # Requirement 1: Divide numerical value by one million
                number_div_1m = number / 1000000.0
                
                # Requirement 2: Calculate percentage
                percentage = (count / total_count) * 100
                
                # Format output row
                output_line = (
                    f"{material_name}\t"
                    f"{number_div_1m:.6f}\t"
                    f"{count}\t"
                    f"{percentage:.6f}\n"
                )
                f_out.write(output_line)
                
        print("Processing complete.")
        
    except Exception as e:
        print(f"An unknown error occurred during file writing: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()