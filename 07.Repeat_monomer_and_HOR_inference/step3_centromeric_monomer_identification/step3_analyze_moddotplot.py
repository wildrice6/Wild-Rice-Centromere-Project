#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
from collections import Counter
import sys

def parse_region_string(region_str):
    """
    Parses a string in 'Chr:start-end' format.
    Returns (chromosome, start, end).
    """
    match = re.match(r'([^:]+):(\d+)-(\d+)', region_str)
    if match:
        chrom, start, end = match.groups()
        return chrom, int(start), int(end)
    return None, None, None

def format_region_tuple(region_tuple):
    """
    Formats a (chromosome, start, end) tuple back into a string.
    """
    return f"{region_tuple[0]}:{region_tuple[1]}-{region_tuple[2]}"

def analyze_moddotplot(bed_file, top_n, output_prefix):
    """
    Main analysis function to identify dominant repeat families from moddotplot output.
    """
    print(f"--- Starting analysis for file: {bed_file} ---")

    # 1. Reading and filtering data (Identity >= 86%)
    print("Step 1/5: Reading and filtering input file (Identity >= 86%)...")
    high_identity_records = []
    cen_info = None # To store centromere region information (chr, start, end)

    try:
        with open(bed_file, 'r') as f:
            for i, line in enumerate(f, 1):
                if line.startswith('#') or not line.strip():
                    continue
                
                parts = line.strip().split()
                if len(parts) < 7:
                    print(f"Warning: Line {i} format is incorrect, skipping: {line.strip()}", file=sys.stderr)
                    continue

                try:
                    identity = float(parts[6])
                    if identity >= 86.0:
                        if cen_info is None:
                            # Extract centromere global coordinates from the first valid record
                            cen_info = parse_region_string(parts[0])
                            if cen_info[0] is None:
                                print(f"Error: Unable to parse centromere region format in the first column: {parts[0]}", file=sys.stderr)
                                sys.exit(1)
                        
                        record = {
                            "q_name": parts[0], "q_start": int(parts[1]), "q_end": int(parts[2]),
                            "r_name": parts[3], "r_start": int(parts[4]), "r_end": int(parts[5]),
                            "identity": identity
                        }
                        high_identity_records.append(record)
                except (ValueError, IndexError):
                    print(f"Warning: Data at line {i} could not be parsed, skipping: {line.strip()}", file=sys.stderr)
                    continue
    except FileNotFoundError:
        print(f"Error: Input file not found: {bed_file}", file=sys.stderr)
        sys.exit(1)

    if not high_identity_records:
        print("Error: No valid data lines found with Identity >= 86%.", file=sys.stderr)
        sys.exit(1)
        
    cen_chr, cen_start, cen_end = cen_info
    cen_length = cen_end - cen_start + 1
    print(f"Found {len(high_identity_records)} high-identity records.")
    print(f"Centromere region: {format_region_tuple(cen_info)}, Total length: {cen_length} bp")

    # 2. Coordinate transformation and Result File 1 generation
    abs_coords_file = f"{output_prefix}_abs_coords.txt"
    print(f"\nStep 2/5: Transforming coordinates and generating file -> {abs_coords_file}")
    
    abs_coord_pairs = []
    with open(abs_coords_file, 'w') as f_out:
        f_out.write("region1\tregion2\tidentity\n")
        for record in high_identity_records:
            # Calculate absolute genome coordinates based on the centromere offset
            abs_q_start = cen_start + record["q_start"] - 1
            abs_q_end = cen_start + record["q_end"] - 1
            abs_r_start = cen_start + record["r_start"] - 1
            abs_r_end = cen_start + record["r_end"] - 1
            q_region_tuple = (cen_chr, abs_q_start, abs_q_end)
            r_region_tuple = (cen_chr, abs_r_start, abs_r_end)
            abs_coord_pairs.append((q_region_tuple, r_region_tuple))
            f_out.write(f"{format_region_tuple(q_region_tuple)}\t{format_region_tuple(r_region_tuple)}\t{record['identity']:.4f}\n")
    print("Coordinate transformation complete.")

    # 3. Iterative calculation of Top N repeats
    print(f"\nStep 3/5: Iteratively calculating Top {top_n} repeat families...")
    top_repeats_results = []
    all_top_family_regions = []
    
    remaining_pairs = list(abs_coord_pairs)
    cumulative_percentage_total = 0.0

    for i in range(1, top_n + 1):
        if not remaining_pairs:
            print(f"Warning: No remaining data to analyze after calculating Top {i-1}.")
            break

        print(f"  -> Calculating Top {i}...")
        
        all_regions_in_scope = []
        for q_region, r_region in remaining_pairs:
            all_regions_in_scope.append(q_region)
            all_regions_in_scope.append(r_region)
        
        if not all_regions_in_scope:
            break 

        # Identify the most frequent region as the 'seed' to define the current family
        region_counts = Counter(all_regions_in_scope)
        top_seed_tuple, _ = region_counts.most_common(1)[0]
        
        # Build the family: include the seed and all associated regions
        current_family_regions = set()
        for q_region, r_region in remaining_pairs:
            if q_region == top_seed_tuple or r_region == top_seed_tuple:
                current_family_regions.add(q_region)
                current_family_regions.add(r_region)

        # Logic Update: Copy number is defined as the total number of unique segments in the family
        copy_number = len(current_family_regions)

        # Collect current family region information for Result File 3
        top_id_str = f"Top{i}"
        for region_tuple in current_family_regions:
            all_top_family_regions.append({'region': region_tuple, 'top_id': top_id_str})

        # Calculate family cumulative length and centromere coverage
        cumulative_length = sum(end - start + 1 for _, start, end in current_family_regions)
        percentage_of_cen = (cumulative_length / cen_length) * 100
        cumulative_percentage_total += percentage_of_cen
        
        result_line = {
            "region": format_region_tuple(top_seed_tuple),
            "copy_number": copy_number,
            "cumulative_length": cumulative_length,
            "percentage_of_cen": percentage_of_cen,
            "cumulative_percentage": cumulative_percentage_total
        }
        top_repeats_results.append(result_line)

        # Remove pairs already assigned to this family for the next iteration
        remaining_pairs = [
            (q, r) for q, r in remaining_pairs 
            if q not in current_family_regions and r not in current_family_regions
        ]
    
    print("Top N calculation complete.")

    # 4. Result File 2 generation (Statistics)
    top_repeats_file = f"{output_prefix}_top{top_n}_repeats.txt"
    print(f"\nStep 4/5: Generating Top N statistics file -> {top_repeats_file}")

    with open(top_repeats_file, 'w') as f_out:
        header = "Top\tSeed_Region\tCopy_Number\tFamily_Cumulative_Length\tPercentage_of_Centromere\tCumulative_Percentage\n"
        f_out.write(header)
        for i, result in enumerate(top_repeats_results, 1):
            f_out.write(
                f"Top{i}\t"
                f"{result['region']}\t"
                f"{result['copy_number']}\t"
                f"{result['cumulative_length']}\t"
                f"{result['percentage_of_cen']:.4f}%\t"
                f"{result['cumulative_percentage']:.4f}%\n"
            )

    # 5. Result File 3 generation (Region List)
    regions_list_file = f"{output_prefix}_top{top_n}_regions_list.txt"
    print(f"\nStep 5/5: Generating Top N regions list file -> {regions_list_file}")

    # Sort the region list for better readability
    all_top_family_regions.sort(key=lambda x: x['region'])

    with open(regions_list_file, 'w') as f_out:
        f_out.write("Region\tTop_Family\n")
        for item in all_top_family_regions:
            region_str = format_region_tuple(item['region'])
            f_out.write(f"{region_str}\t{item['top_id']}\n")

    print(f"\n--- Analysis complete ---")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Quantify moddotplot results to identify major repeat families within centromeres.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("bed_file", help="Input BED file containing moddotplot results.")
    parser.add_argument("top_n", type=int, help="Number of Top N repeat families to identify (e.g., 10).")
    parser.add_argument("output_prefix", help="Prefix for the generated output files (e.g., 'Chr1_cen_analysis').")
    args = parser.parse_args()
    
    analyze_moddotplot(args.bed_file, args.top_n, args.output_prefix)