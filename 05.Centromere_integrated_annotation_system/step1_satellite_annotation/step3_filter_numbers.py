import os
import glob

def count_records_in_bed(filepath):
    """
    Calculates the number of records (i.e., line count) in a BED file.
    """
    count = 0
    try:
        # Utilizing the 'with' statement ensures proper file closure.
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            # Iterating through the file line by line is more memory-efficient 
            # than loading the entire file into memory.
            for _ in f:
                count += 1
    except Exception as e:
        print(f"  -> Error: Unable to read file {filepath}: {e}")
        return -1  # Return an error code indicating a read failure.
    return count

def process_and_clean_bed_files():
    """
    Main function: Locate, inspect, and prune BED files in the current directory 
    based on record count.
    """
    # Retrieve the current working directory.
    current_directory = os.getcwd()
    print(f"Commencing directory scan: {current_directory}\n")

    # Utilize glob to identify all files ending in .bed (case-insensitive).
    bed_files = glob.glob('*.bed') + glob.glob('*.BED')
    
    # Eliminate potential duplicates and sort the file list.
    bed_files = sorted(list(set(bed_files)))

    if not bed_files:
        print("No *.bed files were found in the current directory.")
        return

    print(f"Found {len(bed_files)} BED files; initiating inspection...\n" + "-"*30)

    # Iterate through all identified BED files.
    for filename in bed_files:
        filepath = os.path.join(current_directory, filename)

        # 1. Calculate the record count within the file.
        record_count = count_records_in_bed(filepath)

        # Skip the file if the read operation fails.
        if record_count == -1:
            print("-" * 30)
            continue

        print(f"Inspecting file: '{filename}' ... contains {record_count} records.")

        # 2. Determine if the record count is below 20 and execute the corresponding action.
        if record_count < 20:
            try:
                os.remove(filepath)
                print(f"  -> Deleted. (Reason: Record count {record_count} < 20)\n")
            except OSError as e:
                print(f"  -> Deletion failed! Error: {e}\n")
        else:
            print(f"  -> Retained. (Reason: Record count {record_count} >= 20)\n")
        
        print("-" * 30)

    print("All file inspections completed.")

# Execute the main function when the script is run directly.
if __name__ == "__main__":
    process_and_clean_bed_files()