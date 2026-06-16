import os
import argparse
import pandas as pd
import sys
from collections import defaultdict
import shutil
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import numpy as np

# --- Argument Parser Setup ---
parser = argparse.ArgumentParser(description="Process HiC-Pro output matrix files")
parser.add_argument('--matrix', type=str, required=True, help='Matrix file path')
parser.add_argument('--abs', type=str, required=True, help='Bin information file (abs file)')
parser.add_argument('--cen', type=str, required=True, help='Centromere regions file')

args = parser.parse_args()

def filter_intra_chromosomal_contacts(matrix_file, abs_file, output_dir):
    """
    Extracts intra-chromosomal contact records for each chromosome from the Hi-C contact matrix.

    Args:
        matrix_file (str): Path to the matrix file.
                           Format: bin1_id  bin2_id  contact_value (tab or space separated)
        abs_file (str): Path to the abs file.
                        Format: chr_name  start  end  bin_id (tab or space separated)
        output_dir (str): Output directory path to store results.
    """
    # --- Step 1: Read abs file to build bin -> chromosome mapping ---
    print(f"Reading abs file: {abs_file}...")
    bin_to_chr = {}
    try:
        with open(abs_file, 'r') as f:
            for line in f:
                # Use split() to handle tabs or multiple spaces
                parts = line.strip().split()
                if len(parts) < 4:
                    continue  # Skip malformed lines
                chrom = parts[0]
                bin_id = int(parts[3])
                bin_to_chr[bin_id] = chrom
    except FileNotFoundError:
        print(f"Error: abs file not found -> {abs_file}", file=sys.stderr)
        return
    except ValueError:
        print(f"Error: bin_id '{parts[3]}' in abs file is not a valid integer.", file=sys.stderr)
        return
        
    print(f"Successfully built mapping for {len(bin_to_chr)} bins.")

    # --- Step 2: Create output directory ---
    print(f"Preparing output directory: {output_dir}...")
    os.makedirs(output_dir, exist_ok=True)
    
    # --- Step 3: Process matrix file and write output ---
    print(f"Processing matrix file: {matrix_file}...")
    
    # Store file handles for each chromosome to avoid frequent re-opening
    output_file_handlers = {}
    
    line_count = 0
    intra_contact_count = 0

    try:
        with open(matrix_file, 'r') as f_matrix:
            for line in f_matrix:
                line_count += 1
                if line_count % 5000000 == 0:  # Print progress every 5 million lines
                    print(f"  ...Processed {line_count} lines...")
                    
                parts = line.strip().split()
                if len(parts) < 3:
                    continue # Skip malformed lines

                try:
                    bin1 = int(parts[0])
                    bin2 = int(parts[1])
                except ValueError:
                    print(f"Warning: Skipping invalid bin ID at line {line_count}: {line.strip()}", file=sys.stderr)
                    continue

                # Look up chromosome for both bins
                chr1 = bin_to_chr.get(bin1)
                chr2 = bin_to_chr.get(bin2)
                
                # Core logic: if both bins exist in mapping and belong to the same chromosome
                if chr1 and chr1 == chr2:
                    intra_contact_count += 1
                    output_chr = chr1
                    
                    # Open file handle if it doesn't exist
                    if output_chr not in output_file_handlers:
                        output_path = os.path.join(output_dir, f"{output_chr}.tsv")
                        # Open in write mode 'w'
                        output_file_handlers[output_chr] = open(output_path, 'w')
                        print(f"  -> Created file: {output_path}")

                    # Write the original line to the corresponding file
                    output_file_handlers[output_chr].write(line)

    except FileNotFoundError:
        print(f"Error: matrix file not found -> {matrix_file}", file=sys.stderr)
        return
    finally:
        # --- Step 4: Close all open output files ---
        print("Closing all output files...")
        for handler in output_file_handlers.values():
            handler.close()
            
    print("\nProcessing complete!")
    print(f"Total matrix records processed: {line_count}")
    print(f"Intra-chromosomal contact records found: {intra_contact_count}")
    print(f"Results saved in directory: {output_dir}")


def calculate_marginal_sums(input_dir):
    """
    Iterates through all TSV files in a directory and calculates the Marginal Sum for each bin.
    Marginal Sum = sum of contact values between a bin and all other bins (including itself) on the same chromosome.

    Args:
        input_dir (str): Path to folder containing contact TSV files.
                         Input format: bin1  bin2  contact_value (tab separated)
    
    Returns:
        None. Results are written directly to the input folder as [filename].eachbin_contact.tsv.
    """
    # --- Step 1: Check if directory exists ---
    if not os.path.isdir(input_dir):
        print(f"Error: Directory '{input_dir}' does not exist.", file=sys.stderr)
        return

    print(f"Starting processing in directory: {input_dir}")

    # --- Step 2: Iterate through all files in the directory ---
    # Filter out result files to prevent recursive processing
    files_to_process = [
        f for f in os.listdir(input_dir)
        if f.endswith('.tsv') and not f.endswith('.eachbin_sum.tsv') and not f.endswith('.eachbin_contact.tsv')
    ]
    
    if not files_to_process:
        print("No .tsv files found to process in the directory.")
        return

    for filename in files_to_process:
        input_file_path = os.path.join(input_dir, filename)
        
        if not os.path.isfile(input_file_path):
            continue
            
        print(f"\n--- Processing file: {filename} ---")

        # --- Step 3: Calculate Marginal Sums ---
        bin_contact_sums = defaultdict(float)

        try:
            with open(input_file_path, 'r') as f:
                for i, line in enumerate(f):
                    parts = line.strip().split()
                    
                    if len(parts) < 3:
                        continue
                    
                    try:
                        bin1 = int(parts[0])
                        bin2 = int(parts[1])
                        contact_value = float(parts[2])
                    except ValueError:
                        print(f"  Warning: Data type error at line {i+1}, skipping.")
                        continue

                    # Core logic: Accumulate interaction values (considering matrix symmetry)
                    # 1. Accumulate for bin1
                    bin_contact_sums[bin1] += contact_value
                    
                    # 2. Accumulate for bin2 (exclude diagonal self-interaction to avoid double counting)
                    if bin1 != bin2:
                        bin_contact_sums[bin2] += contact_value
        
        except FileNotFoundError:
            print(f"  Error: Cannot read file {input_file_path}", file=sys.stderr)
            continue
            
        print(f"  ...Read complete, generating results.")

        # --- Step 4: Prepare output data ---
        if not bin_contact_sums:
            print("  No valid data in file, skipping output generation.")
            continue

        results = []
        for bin_id, total_sum in bin_contact_sums.items():
            results.append((bin_id, total_sum))
            
        # Sort by bin ID
        results.sort(key=lambda x: x[0])
        
        # --- Step 5: Write to new TSV file ---
        base_name = os.path.splitext(filename)[0]
        output_filename = f"{base_name}.eachbin_contact.tsv" 
        output_file_path = os.path.join(input_dir, output_filename)
        
        try:
            with open(output_file_path, 'w') as f_out:
                f_out.write("bin_id\tmarginal_sum\n")
                # Write data
                for bin_id, sum_val in results:
                    f_out.write(f"{bin_id}\t{sum_val:.6f}\n")
            print(f"  -> Results successfully saved to: {output_filename}")
        except IOError as e:
            print(f"  Error: Could not write output file {output_file_path}. Reason: {e}", file=sys.stderr)

    print("\nAll files processed!")

def map_bin_locations(input_dir, matrix_file):
    """
    Merges bin average contact files with bin absolute position information.

    Args:
        input_dir (str): Folder containing '.eachbin_contact.tsv' files.
        matrix_file (str): File containing bin absolute positions.
                           Format: chr start end bin_id (tab or space separated).
    """
    # --- Step 1: Check input paths ---
    if not os.path.isdir(input_dir):
        print(f"Error: Input directory '{input_dir}' does not exist.", file=sys.stderr)
        return
    if not os.path.isfile(matrix_file):
        print(f"Error: Position file '{matrix_file}' does not exist.", file=sys.stderr)
        return

    # --- Step 2: Read position file, build bin_id -> [chr, start, end] mapping ---
    print(f"Reading position file: {matrix_file} ...")
    bin_to_location = {}
    try:
        with open(matrix_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) < 4:
                    continue
                try:
                    # bin_id is the 4th column, locations are first 3 columns
                    bin_id = int(parts[3])
                    location_info = parts[0:3]  # [chr, start, end]
                    bin_to_location[bin_id] = location_info
                except ValueError:
                    print(f"Warning: Invalid bin_id found in position file, skipping: '{line.strip()}'", file=sys.stderr)
                    continue
    except Exception as e:
        print(f"Error: Failed to read position file: {e}", file=sys.stderr)
        return
    
    if not bin_to_location:
        print("Warning: No data loaded from position file. Check format and content.", file=sys.stderr)
        return
        
    print(f"Position mapping complete. Information for {len(bin_to_location)} bins loaded.")

    # --- Step 3: Iterate directory and process each .eachbin_contact.tsv file ---
    print(f"\nScanning directory: {input_dir}")
    
    files_processed = 0
    for filename in os.listdir(input_dir):
        if filename.endswith(".eachbin_contact.tsv"):
            files_processed += 1
            input_path = os.path.join(input_dir, filename)
            print(f"--- Processing: {filename} ---")
            
            # Generate output filename
            base_name = os.path.splitext(filename)[0] 
            output_filename = f"{base_name}.abs"
            output_path = os.path.join(input_dir, output_filename)

            try:
                # Open both files using 'with' statement
                with open(input_path, 'r') as f_in, open(output_path, 'w') as f_out:
                    # Write header
                    f_out.write("chr\tstart\tend\tavg_contact\tbin_id\n")
                    
                    # Skip header of input file
                    try:
                        next(f_in)
                    except StopIteration:
                        print("  File is empty, skipping.")
                        continue

                    # Process lines
                    for line in f_in:
                        parts = line.strip().split()
                        if len(parts) < 2:
                            continue
                        
                        try:
                            bin_id = int(parts[0])
                            contact_val = parts[1]
                        except ValueError:
                            print(f"  Warning: Invalid data line skipped: '{line.strip()}'", file=sys.stderr)
                            continue
                        
                        # Find position in mapping
                        location_info = bin_to_location.get(bin_id)
                        
                        if location_info:
                            chrom, start, end = location_info
                            f_out.write(f"{chrom}\t{start}\t{end}\t{contact_val}\t{bin_id}\n")
                        else:
                            print(f"  Warning: bin_id '{bin_id}' not found in position file, skipping.", file=sys.stderr)

                print(f"  -> Results successfully saved to: {output_filename}")

            except Exception as e:
                print(f"  Error processing file {filename}: {e}", file=sys.stderr)

    if files_processed == 0:
        print("No files ending with '.eachbin_contact.tsv' found in the directory.")
    
    print("\nAll files processed!")

def split_chromosome_arms(cen_file, genome_size_file, output_tsv):
    """
    Divides chromosome arms into 50 equal bins based on centromere positions.

    Args:
        cen_file (str): BED file path containing centromere positions.
                        Format: chr start end
        genome_size_file (str): .fai file generated by samtools faidx.
                                Format: chr length ...
        output_tsv (str): Output TSV file path.
    """
    # --- Step 1: Read and parse input files ---
    print("Step 1: Reading input files...")

    # Read genome sizes (.fai)
    genome_lengths = {}
    try:
        with open(genome_size_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    chrom = parts[0]
                    length = int(parts[1])
                    genome_lengths[chrom] = length
    except FileNotFoundError:
        print(f"Error: Genome size file '{genome_size_file}' not found.", file=sys.stderr)
        return
    except ValueError:
        print(f"Error: Incorrect format in genome size file '{genome_size_file}'.", file=sys.stderr)
        return
    print(f"  -> Successfully loaded length info for {len(genome_lengths)} chromosomes.")

    # Read centromere positions (BED)
    centromeres = {}
    try:
        with open(cen_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 3:
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    centromeres[chrom] = (start, end)
    except FileNotFoundError:
        print(f"Error: Centromere file '{cen_file}' not found.", file=sys.stderr)
        return
    except ValueError:
        print(f"Error: Incorrect format in centromere file '{cen_file}'.", file=sys.stderr)
        return
    print(f"  -> Successfully loaded centromere positions for {len(centromeres)} chromosomes.")

    # --- Step 2: Calculate and generate all regions ---
    print("\nStep 2: Calculating chromosome arm partitions...")
    
    all_regions = []
    sorted_chroms = sorted(genome_lengths.keys())

    for chrom in sorted_chroms:
        if chrom not in centromeres:
            print(f"  Warning: Chromosome '{chrom}' missing in centromere file, skipping.", file=sys.stderr)
            continue

        chr_len = genome_lengths[chrom]
        cen_start, cen_end = centromeres[chrom]
        
        print(f"  Processing {chrom}...")

        # --- Left arm (p-arm) ---
        left_arm_len = cen_start
        if left_arm_len > 0:
            left_bin_size = left_arm_len / 50.0
            for i in range(50):
                start = int(i * left_bin_size)
                end = int((i + 1) * left_bin_size)
                if i == 49:
                    end = cen_start
                region_id = f"left_arm_{i + 1}"
                all_regions.append((chrom, start, end, region_id))

        # --- Centromere ---
        all_regions.append((chrom, cen_start, cen_end, "cen"))

        # --- Right arm (q-arm) ---
        right_arm_len = chr_len - cen_end
        if right_arm_len > 0:
            right_bin_size = right_arm_len / 50.0
            for i in range(50):
                start = int(cen_end + (i * right_bin_size))
                end = int(cen_end + ((i + 1) * right_bin_size))
                if i == 49:
                    end = chr_len
                region_id = f"right_arm_{i + 1}"
                all_regions.append((chrom, start, end, region_id))

    # --- Step 3: Write output file ---
    print(f"\nStep 3: Writing results to '{output_tsv}'...")
    try:
        output_dir = os.path.dirname(output_tsv)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            
        with open(output_tsv, 'w') as f_out:
            f_out.write("chr\tstart\tend\tregion_id\n")
            for region in all_regions:
                f_out.write(f"{region[0]}\t{region[1]}\t{region[2]}\t{region[3]}\n")
    except IOError as e:
        print(f"Error: Could not write output file. Reason: {e}", file=sys.stderr)
        return

    print("\nProcessing complete!")

def merge_abs_files(input_dir, output_file):
    """
    Merges all files ending in '.eachbin_contact.abs' within a folder.

    Args:
        input_dir (str): Folder containing '.eachbin_contact.abs' files.
        output_file (str): Full path and name for the merged output file.
    """
    if not os.path.isdir(input_dir):
        print(f"Error: Directory '{input_dir}' does not exist.", file=sys.stderr)
        return

    print(f"Scanning directory: {input_dir}")
    suffix = ".eachbin_contact.abs"

    # Find matching files
    files_to_merge = [
        os.path.join(input_dir, f)
        for f in os.listdir(input_dir)
        if f.endswith(suffix) and os.path.isfile(os.path.join(input_dir, f))
    ]

    if not files_to_merge:
        print(f"No files found ending in '{suffix}'.")
        return
    
    files_to_merge.sort()

    print(f"Found {len(files_to_merge)} files to merge:")
    for f_path in files_to_merge:
        print(f"  - {os.path.basename(f_path)}")

    # Security check: prevent merging output into itself
    output_filepath = output_file
    if os.path.abspath(output_filepath) in [os.path.abspath(f) for f in files_to_merge]:
        print(f"Warning: Output file '{os.path.basename(output_filepath)}' is in merge list, removing it.")
        files_to_merge = [f for f in files_to_merge if os.path.abspath(f) != os.path.abspath(output_filepath)]
        if not files_to_merge:
            print("No other files left to merge.")
            return

    print(f"\nMerged file will be saved to: {output_filepath}")

    # Perform merge
    try:
        output_dir = os.path.dirname(output_filepath)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)

        with open(output_filepath, 'w') as f_out:
            is_header_written = False
            for file_path in files_to_merge:
                with open(file_path, 'r') as f_in:
                    if not is_header_written:
                        shutil.copyfileobj(f_in, f_out)
                        is_header_written = True
                    else:
                        next(f_in) # Skip header for subsequent files
                        shutil.copyfileobj(f_in, f_out)
                f_out.write('\n') # Ensure line break between files

    except IOError as e:
        print(f"Error: File operation failed. Reason: {e}", file=sys.stderr)
        return

    print("\nMerge complete!")


def bin_chromosome_coordinates(contact_abs, cen, num_bins_per_arm=50):
    """
    Dividing the left and right arms into a specified number of windows using the 
    centromere midpoint as the origin, and calculating average contact for each.

    Args:
        contact_abs (str): Path to input TSV.
        cen (str): Centromere BED file path.
        num_bins_per_arm (int): Number of windows per arm. Total windows = 2 * num_bins_per_arm.
    """
    total_bins = num_bins_per_arm * 2
    
    # --- Step 1: Read Centromere info ---
    print("Step 1: Reading centromere info...")
    centromeres = {}
    if not os.path.isfile(cen):
        print(f"Error: File '{cen}' does not exist.", file=sys.stderr)
        return

    try:
        with open(cen, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 3:
                    chrom = parts[0]
                    start, end = int(parts[1]), int(parts[2])
                    cen_mid = (start + end) / 2.0
                    centromeres[chrom] = cen_mid
    except Exception as e:
        print(f"Error reading centromere file: {e}", file=sys.stderr)
        return

    # --- Step 2: Scan for chromosome lengths ---
    print("Step 2: Scanning chromosome lengths...")
    chr_lengths = {}
    if not os.path.isfile(contact_abs):
        print(f"Error: File '{contact_abs}' does not exist.", file=sys.stderr)
        return

    try:
        with open(contact_abs, 'r') as f:
            next(f) # Skip header
            for line in f:
                parts = line.strip().split()
                if len(parts) < 3: continue
                chrom = parts[0]
                try:
                    end_pos = int(parts[2])
                    if chrom not in chr_lengths or end_pos > chr_lengths[chrom]:
                        chr_lengths[chrom] = end_pos
                except ValueError: continue
    except Exception as e:
        print(f"Error reading contact file: {e}", file=sys.stderr)
        return

    # --- Step 3: Aggregate data into windows ---
    print(f"Step 3: Aggregating data into {total_bins} windows ({num_bins_per_arm} left + {num_bins_per_arm} right)...")
    
    bin_data = defaultdict(lambda: defaultdict(lambda: {'sum': 0.0, 'count': 0}))
    
    try:
        with open(contact_abs, 'r') as f_in:
            next(f_in)
            for line in f_in:
                parts = line.strip().split()
                if len(parts) < 5: continue
                
                chrom = parts[0]
                if chrom not in centromeres or chrom not in chr_lengths:
                    continue

                try:
                    start, end = int(parts[1]), int(parts[2])
                    contact_val = float(parts[3]) 
                except ValueError:
                    continue

                cen_mid = centromeres[chrom]
                chr_len = chr_lengths[chrom]
                bin_midpoint = (start + end) / 2.0
                
                window_idx = -1

                # --- Window mapping logic ---
                if bin_midpoint < cen_mid:
                    # Left arm: ratio maps to index 0 ~ 49
                    ratio = bin_midpoint / cen_mid if cen_mid > 0 else 0
                    idx = int(ratio * num_bins_per_arm)
                    if idx >= num_bins_per_arm:
                        idx = num_bins_per_arm - 1
                    window_idx = idx
                else:
                    # Right arm: ratio maps to index 50 ~ 99
                    right_len = chr_len - cen_mid
                    ratio = (bin_midpoint - cen_mid) / right_len if right_len > 0 else 0
                    idx = num_bins_per_arm + int(ratio * num_bins_per_arm)
                    if idx >= total_bins:
                        idx = total_bins - 1
                    window_idx = idx

                if 0 <= window_idx < total_bins:
                    bin_data[chrom][window_idx]['sum'] += contact_val
                    bin_data[chrom][window_idx]['count'] += 1

    except Exception as e:
        print(f"Error processing data: {e}", file=sys.stderr)
        return

    # --- Step 4: Calculate averages and write output ---
    base_name = os.path.splitext(os.path.basename(contact_abs))[0]
    output_file = os.path.join(os.path.dirname(contact_abs), f"{base_name}.normalized.tsv")

    print(f"Step 4: Writing results to '{output_file}'...")

    try:
        with open(output_file, 'w') as f_out:
            f_out.write("Chrom\tWindow_Index\tNormalized_Coordinate\tAvg_Contact_Value\n")
            for chrom in sorted(bin_data.keys()):
                for i in range(total_bins):
                    stats = bin_data[chrom].get(i, {'sum': 0.0, 'count': 0})
                    avg_val = stats['sum'] / stats['count'] if stats['count'] > 0 else 0.0
                    
                    # Compute normalized coordinate for X-axis plotting (-1 to 1)
                    norm_coord = (i - num_bins_per_arm) / float(num_bins_per_arm)
                    f_out.write(f"{chrom}\t{i}\t{norm_coord:.4f}\t{avg_val:.6f}\n")
                    
    except Exception as e:
        print(f"Error writing file: {e}", file=sys.stderr)
        return

    print("Processing complete!")


def plot_chromosome_contacts(input_abs_contact, output_dir):
    """
    Plots contact profile line charts for each chromosome in the TSV file.
    """
    if not os.path.isfile(input_abs_contact):
        print(f"Error: Input file '{input_abs_contact}' not found.", file=sys.stderr)
        return
        
    try:
        os.makedirs(output_dir, exist_ok=True)
        print(f"Output directory '{output_dir}' ready.")
    except OSError as e:
        print(f"Error: Failed to create output directory '{output_dir}'. Reason: {e}", file=sys.stderr)
        return

    print(f"Reading data file: {input_abs_contact}")
    try:
        df = pd.read_csv(input_abs_contact, sep='\t')
        df['Normalized_Coordinate'] = pd.to_numeric(df['Normalized_Coordinate'])
        df['Avg_Contact_Value'] = pd.to_numeric(df['Avg_Contact_Value'])
    except Exception as e:
        print(f"Error parsing file. Reason: {e}", file=sys.stderr)
        return

    chromosomes = df['Chrom'].unique()
    print(f"Found {len(chromosomes)} chromosomes. Starting plotting...")

    for chrom in chromosomes:
        chrom_df = df[df['Chrom'] == chrom].copy()
        chrom_df.sort_values('Normalized_Coordinate', inplace=True)
        
        print(f"  -> Plotting: {chrom}")

        try:
            plt.style.use('seaborn-v0_8-whitegrid')
        except OSError:
            plt.style.use('default') 

        fig, ax = plt.subplots(figsize=(10, 5)) 

        sns.lineplot(
            x='Normalized_Coordinate', 
            y='Avg_Contact_Value', 
            data=chrom_df, 
            ax=ax,
            linewidth=2.5,
            color='black' 
        )
        
        ax.set_title(f"{chrom}", fontsize=14, weight='bold')
        ax.set_xlabel(None)
        ax.set_ylabel("Average Contact Value", fontsize=12)
        
        if ax.get_legend() is not None:
            ax.get_legend().remove()
        
        ax.set_xlim(-1.1, 1.1)
        y_max = chrom_df['Avg_Contact_Value'].max()
        ax.set_ylim(bottom=0, top=y_max * 1.1 if y_max > 0 else 1)
        
        # Set specific ticks for Telomere and Centromere labels
        ax.set_xticks([-1.0, 0.0, 1.0])
        ax.set_xticklabels(['TEL', 'CEN', 'TEL'], fontsize=12, fontweight='bold')

        base_filename = os.path.join(output_dir, chrom)
        plt.savefig(f"{base_filename}.png", dpi=300, bbox_inches='tight')
        plt.savefig(f"{base_filename}.pdf", bbox_inches='tight')
        plt.close(fig)

    print("\nAll plotting tasks finished!")



def batch_create_heatmaps(input_dir: str, output_dir: str):
    """
    Batch processes TSV files to generate heatmaps with log10(x+1) transformation.
    """
    if not os.path.isdir(input_dir):
        print(f"ERROR: Input folder '{input_dir}' missing.")
        return

    try:
        os.makedirs(output_dir, exist_ok=True)
        print(f"INFO: Output saved to '{output_dir}'")
    except OSError as e:
        print(f"ERROR: Cannot create output folder. Error: {e}")
        return

    search_pattern = os.path.join(input_dir, '*.tsv')
    tsv_files = glob.glob(search_pattern)
    
    if not tsv_files:
        search_pattern = os.path.join(input_dir, '*.txt')
        tsv_files = glob.glob(search_pattern)

    if not tsv_files:
        print(f"WARNING: No .tsv or .txt files found in '{input_dir}'.")
        return
        
    print(f"INFO: Found {len(tsv_files)} files to process.")

    for input_file_path in tsv_files:
        try:
            base_filename = os.path.splitext(os.path.basename(input_file_path))[0]
            print(f"\n--- Processing: {base_filename}.tsv ---")

            df = pd.read_csv(input_file_path, sep='\s+', header=None, names=['bin1', 'bin2', 'value'])
            
            if df.empty:
                print(f"WARNING: File '{input_file_path}' is empty, skipping.")
                continue

            # Reshape to matrix
            heatmap_df = df.pivot(index='bin2', columns='bin1', values='value')
            heatmap_df = heatmap_df.sort_index(axis=0).sort_index(axis=1)
            heatmap_df = heatmap_df.fillna(0)

            # Apply log transformation for Hi-C data dynamic range compression
            heatmap_log_values = np.log10(heatmap_df.values + 1)

            fig, ax = plt.subplots(figsize=(10, 10))
            ax.imshow(heatmap_log_values, cmap='Reds', origin='lower', aspect='auto')
            ax.axis('off')
            
            output_png = os.path.join(output_dir, f"{base_filename}.png")
            output_pdf = os.path.join(output_dir, f"{base_filename}.pdf")

            fig.savefig(output_png, bbox_inches='tight', pad_inches=0, dpi=300)
            print(f"SUCCESS: PNG saved to '{output_png}'")
            fig.savefig(output_pdf, bbox_inches='tight', pad_inches=0)
            print(f"SUCCESS: PDF saved to '{output_pdf}'")

        except Exception as e:
            print(f"ERROR: Error processing file '{input_file_path}': {e}")
        finally:
            if 'fig' in locals():
                plt.close(fig)
    
    print("\nINFO: All files processed.")

def main(matrix, abs, cen):
    # Get base name for output file management
    matrix_prefix = matrix.split('.')[0]
    
    # 1. Extract intra-chromosomal contacts
    filter_intra_chromosomal_contacts(matrix, abs, matrix_prefix)
    
    # 2. Create heatmaps for each chromosome
    batch_create_heatmaps(matrix_prefix, f"{matrix_prefix}_heatmap")
    
    # 3. Calculate marginal sums (total contact per bin)
    calculate_marginal_sums(matrix_prefix)
    
    # 4. Map genomic positions to contact values
    map_bin_locations(matrix_prefix, abs)
    
    # 5. Merge processed abs files
    merge_abs_files(matrix_prefix, f"{matrix_prefix}.abs_contact")
    
    # 6. Normalize coordinates relative to centromeres
    bin_chromosome_coordinates(f"{matrix_prefix}.abs_contact", cen)
    
    # 7. Generate final line profile plots
    plot_chromosome_contacts(f"{matrix_prefix}.normalized.tsv", f"{matrix_prefix}_plot")


if __name__ == "__main__":
    main(args.matrix, args.abs, args.cen)