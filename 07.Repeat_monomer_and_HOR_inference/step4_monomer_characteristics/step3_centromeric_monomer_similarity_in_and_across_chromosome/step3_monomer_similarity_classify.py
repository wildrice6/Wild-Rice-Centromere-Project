import argparse
import sys
import re

def get_category(window_str):
    """
    Determines whether a window size belongs to the 'short' or 'long' category 
    based on the provided string.
    100, 150 -> short
    >= 500   -> long
    Others   -> None
    """
    try:
        # Extract the first sequence of digits from the string
        size = int(re.findall(r'\d+', window_str)[0])
    except (IndexError, ValueError):
        return None
    
    if size == 100 or size == 150:
        return "short"
    elif size >= 500:
        return "long"
    else:
        # Does not belong to the two defined categories
        return None

def process_files(sim_file, win_file, output_file):
    # 1. Read window size mapping file
    win_dict = {}
    with open(win_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            assembly, chrom, win_size_str = parts[0], parts[1], parts[2]
            category = get_category(win_size_str)
            win_dict[(assembly, chrom)] = category

    # 2. Process similarity statistics file
    with open(sim_file, 'r') as f_in, open(output_file, 'w') as f_out:
        header = f_in.readline().strip()
        if not header:
            return
        # Append the new column header
        f_out.write(header + "\tcategory_tag\n")
        
        for line in f_in:
            line = line.strip()
            if not line:
                continue
            
            parts = line.split()
            chrA, chrB = parts[0], parts[1]
            comp_type = parts[2]  # 'within' or 'between'
            assembly = parts[5]
            
            catA = win_dict.get((assembly, chrA))
            catB = win_dict.get((assembly, chrB))
            
            # Assign tags only when both chromosomes have defined categories (short/long)
            if catA is None or catB is None:
                tag = "unclassified" 
            elif catA == "short" and catB == "short":
                tag = "short"
            elif catA == "long" and catB == "long":
                tag = "long"
            elif (catA == "short" and catB == "long") or (catA == "long" and catB == "short"):
                # Theoretically, only 'between-chromosome' comparisons should reach this branch
                tag = "mix"
            else:
                # Fallback logic
                tag = "other"

            # Additional safety check: skip if a 'within' comparison is classified as 'mix'
            if comp_type == "within" and tag == "mix":
                # This scenario implies contradictory window size definitions for the same chromosome
                continue 

            f_out.write(f"{line}\t{tag}\n")

def main():
    parser = argparse.ArgumentParser(description="Classify similarity records based on chromosome window sizes.")
    parser.add_argument("--similarity", required=True, help="Path to the similarity statistics file")
    parser.add_argument("--windows", required=True, help="Path to the window size mapping file")
    parser.add_argument("--output", required=True, help="Path to the output file")
    
    args = parser.parse_args()
    
    try:
        process_files(args.similarity, args.windows, args.output)
        print(f"Processing complete! Results saved to: {args.output}")
    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()