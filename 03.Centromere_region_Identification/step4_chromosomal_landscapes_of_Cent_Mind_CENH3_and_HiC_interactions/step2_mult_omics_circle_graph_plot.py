import os
import time
import argparse
import matplotlib.pyplot as plt
from pycirclize import Circos
import numpy as np
import pandas as pd
import tempfile
import shutil

parser = argparse.ArgumentParser(description="Plot Multi-omics Circular Diagram")
parser.add_argument('--hap1', type=str, required=True, help='Haplotype 1, e.g., AA_Ogla_hap1')
parser.add_argument('--hap2', type=str, required=True, help='Haplotype 2, e.g., AA_Ogla_hap2')

args = parser.parse_args()


### ------------------------------------- Public Functions -------------------------------------
def expand_bed(input_bed, output_bed, expand_size=500000):
    """
    Expand each interval in the BED file by a specified size (default 500 kb) on both sides.

    Args:
        input_bed: Path to the input BED file.
        output_bed: Path to the output BED file.
        expand_size: Number of bases to expand (default 500,000 = 500 kb).
    """
    with open(input_bed, 'r') as infile, open(output_bed, 'w') as outfile:
        for line in infile:
            if line.strip() == "" or line.startswith("#"):
                continue
            fields = line.strip().split()
            if len(fields) < 3:
                continue

            chr_name = fields[0]
            start = int(fields[1])
            end = int(fields[2])

            # Expand by expand_size on both sides, ensuring the start position is not less than 0
            new_start = max(0, start - expand_size)
            new_end = end + expand_size

            outfile.write(f"{chr_name}\t{new_start}\t{new_end}\n")
    print(f"✅ Expansion completed. Results saved to: {output_bed}")

def centromere_size_in_mb(bed_file):
    """
    Read the BED file, calculate the size of each centromeric region (in 0.1 Mb units),
    and return a dictionary, e.g., {"Chr01": 16.68, "Chr02": 16.24, ...}.

    Args:
        bed_file: str, Path to the input BED file.
    Returns:
        dict: { Chromosome ID: Region length (0.1 Mb units) }
    """
    result = {}

    with open(bed_file, "r") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.strip().split()
            if len(fields) < 3:
                continue  # Skip malformed lines

            chr_name = fields[0]
            start = int(fields[1])
            end = int(fields[2])

            size_mb = (end - start) / 100_000  # Convert to 0.1 Mb units
            result[chr_name] = round(size_mb, 2)  # Round to two decimal places

    return result

def bedtools_intersect(element, cen, output_file=None):
    if output_file is None:
        output_file = f"{element}.INcen"
    command_intersect = f"bedtools intersect -a {element} -b {cen} > {output_file}"
    os.system(command_intersect)
    print(f"***** Obtained records for {element} within the centromeric region *****")

def set_negative_to_zero(file_path):
    """
    Modify the TSV file by setting values in the fourth column that are less than 0 to 0, 
    overwriting the original file.
    
    Args:
        file_path (str): Path to the TSV file.
    """
    # Read TSV file without header
    df = pd.read_csv(file_path, sep='\t', header=None)
    
    # Set values in the fourth column that are less than 0 to 0
    df[3] = df[3].apply(lambda x: max(x, 0))
    
    # Save back to the original file in TSV format without index
    df.to_csv(file_path, sep='\t', header=False, index=False)

def convert_element_relative(cen_file, element_file, output_file):
    """
    Convert absolute coordinates of the element file to relative coordinates based on the cen file,
    divide by 100,000, and output to a new file.
    
    Args:
        cen_file (str): Path to the cen file (three columns: chr, start, end).
        element_file (str): Path to the element file (at least three columns: chr, start, end, ...).
        output_file (str): Path to the output file.
    """
    # Read cen file
    cen_df = pd.read_csv(cen_file, sep='\t', header=None, names=['chr', 'cen_start', 'cen_end'])
    
    # Create a dictionary to store cen_start for each chromosome
    cen_start_dict = dict(zip(cen_df['chr'], cen_df['cen_start']))
    
    # Read element file
    element_df = pd.read_csv(element_file, sep='\t', header=None)
    
    # Verify if chromosomes in the first 3 columns exist in the cen file
    for chr_name in element_df[0].unique():
        if chr_name not in cen_start_dict:
            raise ValueError(f"Chromosome {chr_name} not found in the cen file!")
    
    # Calculate relative coordinates and divide by 100,000
    element_df[1] = (element_df[1] - element_df[0].map(cen_start_dict)) / 100000
    element_df[2] = (element_df[2] - element_df[0].map(cen_start_dict)) / 100000
    
    # Save as new file, preserving original format
    element_df.to_csv(output_file, sep='\t', header=False, index=False)

def reshape_methylation(file_path, columns_to_remove=[4, 6]):
    """
    Remove specified columns from the TSV file (5th and 7th columns, indices 4 and 6),
    modifying the original file directly.
    
    Args:
        file_path: Path to the TSV file.
        columns_to_remove: List of column indices to remove (0-indexed), default [4, 6].
    
    Returns:
        int: Number of processed lines.
    """
    # Read all content from the file
    with open(file_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    # Check if the file is empty
    if not lines:
        print("File is empty")
        return 0
    
    # Process each line
    processed_lines = []
    for line in lines:
        # Split line into columns using tab delimiter
        columns = line.strip().split('\t')
        
        # Remove specified columns (in reverse order to avoid index shifts)
        for col_idx in sorted(columns_to_remove, reverse=True):
            if col_idx < len(columns):
                columns.pop(col_idx)
        
        # Reassemble the line
        processed_lines.append('\t'.join(columns) + '\n')
    
    # Write modified content back to the original file
    with open(file_path, 'w', encoding='utf-8') as f:
        f.writelines(processed_lines)
    
    print(f"Removed 5th and 7th columns; processed {len(processed_lines)} lines in total.")
    return len(processed_lines)

def remove_header(file_path):
    """
    Check if the first line of the TSV file contains 'num_sites_in_window';
    if so, remove the line, modifying the original file directly.
    
    Args:
        file_path: Path to the TSV file.
    
    Returns:
        bool: True if header was removed, False otherwise.
    """
    # Read all content from the file
    with open(file_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    # Check if the file is empty
    if not lines:
        print("File is empty")
        return False
    
    # Check if the first line contains 'num_sites_in_window'
    if 'num_sites_in_window' in lines[0]:
        # Remove the first line and retain the rest
        lines = lines[1:]
        
        # Write modified content back to the original file
        with open(file_path, 'w', encoding='utf-8') as f:
            f.writelines(lines)
        
        print("Removed header line containing 'num_sites_in_window'.")
        return True
    else:
        print("First line does not contain 'num_sites_in_window'; no removal necessary.")
        return False

def merge_bed_in_place(input_file: str):
    """
    Merge closely adjacent intervals within the BED file and write results back to the original file.
    The output file will retain only three columns (chrom, start, end).

    Merge conditions:
    1. Both records reside on the same chromosome.
    2. End position of the first record + 1 == Start position of the second record.

    Args:
        input_file (str): Path to the BED file to be processed.

    Raises:
        FileNotFoundError: If the input file is not found.
    """
    # 1. Check if the file exists
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Error: Input file not found -> '{input_file}'")

    # 2. Read all records into memory
    records = []
    print(f"Reading data from '{input_file}'...")
    try:
        with open(input_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split('\t')
                if len(parts) < 3:
                    continue  # Skip malformed lines
                
                # --- *** Critical Modification Point 1 *** ---
                # Extract only the first three columns and convert start/end to integers
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                records.append([chrom, start, end])
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    if not records:
        print("File is empty or contains no valid data; no action required.")
        return

    # 3. Sort records: primary by chromosome, secondary by start position
    records.sort(key=lambda x: (x[0], x[1]))

    # 4. Execute merge logic
    print("Merging adjacent intervals...")
    if not records:
        return
        
    merged_records = []
    # Initialize the current merge interval with the first record
    # current_record is now a list containing only three columns [chrom, start, end]
    current_record = records[0]

    for i in range(1, len(records)):
        next_record = records[i]
        
        # Unpack for clarity
        current_chrom, _, current_end = current_record
        next_chrom, next_start, _ = next_record

        # Condition: same chromosome and closely adjacent
        if current_chrom == next_chrom and current_end + 1 == next_start:
            # Merge interval: update the end position of the current record
            current_record[2] = next_record[2]
        else:
            # Merge condition not met: append the current merged interval to results
            merged_records.append(current_record)
            # Start a new merge interval
            current_record = next_record
            
    # Append the last interval processed after the loop
    merged_records.append(current_record)

    # 5. Write merged results to a temporary file, then replace the original file
    try:
        temp_fd, temp_path = tempfile.mkstemp(dir=os.path.dirname(input_file))
        
        with os.fdopen(temp_fd, 'w') as temp_f:
            for record in merged_records:
                # --- *** Critical Modification Point 2 *** ---
                # record now contains only 3 elements; write three columns correctly
                str_record = [str(item) for item in record]
                temp_f.write("\t".join(str_record) + "\n")
        
        os.replace(temp_path, input_file)
        
        print(f"Operation complete! '{input_file}' successfully merged and updated.")
        print(f"Original record count: {len(records)}, Merged record count: {len(merged_records)}")

    except Exception as e:
        print(f"Error writing to file: {e}")
        if 'temp_path' in locals() and os.path.exists(temp_path):
            os.remove(temp_path)

def process_element_file(assembly):
    # Expand the centromeric region by 500 kb on both sides
    cenRegion = f"{assembly}.cenRegion.bed"
    cenRegion_expand500k = f"{assembly}.cenRegion.expand500k.bed"
    expand_bed(cenRegion, cenRegion_expand500k)

    # CENH3 signal
    cenh3_samples = ['sample1', 'sample2', 'sample3']
    for cenh3_sample in cenh3_samples:
        cenh3 = f"{assembly}.{cenh3_sample}.CENH3.bdg"
        cenh3_INcen = f"{assembly}.{cenh3_sample}.CENH3.bdg.INcen"
        cenh3_INcen_Coortrans = f"{assembly}.{cenh3_sample}.CENH3.bdg.INcen.Coortrans"
        bedtools_intersect(cenh3, cenRegion_expand500k)
        # Set negative records to zero
        set_negative_to_zero(cenh3_INcen)
        # Coordinate transformation
        convert_element_relative(cenRegion_expand500k, cenh3_INcen, cenh3_INcen_Coortrans)

        # Regions strictly within centromeres
        cenh3_OnlyINcen = f"{assembly}.{cenh3_sample}.CENH3.bdg.OnlyINcen"
        cenh3_OnlyINcen_Coortrans = f"{assembly}.{cenh3_sample}.CENH3.bdg.OnlyINcen.Coortrans"
        bedtools_intersect(cenh3, cenRegion, cenh3_OnlyINcen)
        set_negative_to_zero(cenh3_OnlyINcen)
        convert_element_relative(cenRegion_expand500k, cenh3_OnlyINcen, cenh3_OnlyINcen_Coortrans)


    # Sequencing depth
    seqs = ['ONT', 'hifi']
    for seq_type in seqs:
        seq = f"{assembly}.{seq_type}.bed"
        seq_INcen = f"{assembly}.{seq_type}.bed.INcen"
        seq_INcen_Coortrans = f"{assembly}.{seq_type}.bed.INcen.Coortrans"
        # Obtain records for the centromere (including 500kb flanks)
        bedtools_intersect(seq, cenRegion_expand500k)
        convert_element_relative(cenRegion_expand500k, seq_INcen, seq_INcen_Coortrans)
        # Obtain records strictly for the centromeric region
        shutil.copy(seq_INcen, f"{assembly}.{seq_type}.OnlyINcen")
        seq_onlyINcen = f"{assembly}.{seq_type}.OnlyINcen.INcen"
        seq_onlyINcen_Coortrans = f"{assembly}.{seq_type}.OnlyINcen.INcen.Coortrans"
        bedtools_intersect(f"{assembly}.{seq_type}.OnlyINcen", cenRegion)
        convert_element_relative(cenRegion_expand500k, seq_onlyINcen, seq_onlyINcen_Coortrans)


    # Hi-C contact matrix information
    hic = f"{assembly}.hic.bed"
    hic_INcen = f"{assembly}.hic.bed.INcen"
    hic_INcen_Coortrans = f"{assembly}.hic.bed.INcen.Coortrans"
    # Obtain records for the centromeric region (including 500kb flanks)
    bedtools_intersect(hic, cenRegion_expand500k)
    convert_element_relative(cenRegion_expand500k, hic_INcen, hic_INcen_Coortrans)
    # Obtain records strictly for the centromeric region
    hic_OnlyINcen = f"{assembly}.hic.bed.OnlyINcen"
    hic_OnlyINcen_Coortrans = f"{assembly}.hic.bed.OnlyINcen.Coortrans"
    bedtools_intersect(hic, cenRegion, hic_OnlyINcen)
    convert_element_relative(cenRegion_expand500k, hic_OnlyINcen, hic_OnlyINcen_Coortrans)


    # Centromere prediction results from the model
    predictCEN = f"{assembly}.predictCEN.bed"
    predictCEN_INcen = f"{assembly}.predictCEN.bed.INcen"
    predictCEN_INcen_Coortrans = f"{assembly}.predictCEN.bed.INcen.Coortrans"
    
    # Obtain records for the centromeric region (including 500kb flanks)
    bedtools_intersect(predictCEN, cenRegion_expand500k)
    convert_element_relative(cenRegion_expand500k, predictCEN_INcen, predictCEN_INcen_Coortrans)
    
    # Obtain records strictly for the centromeric region
    predictCEN_OnlyINcen = f"{assembly}.predictCEN.bed.OnlyINcen"
    predictCEN_OnlyINcen_Coortrans = f"{assembly}.predictCEN.bed.OnlyINcen.Coortrans"
    
    bedtools_intersect(predictCEN, cenRegion, predictCEN_OnlyINcen)
    convert_element_relative(cenRegion_expand500k, predictCEN_OnlyINcen, predictCEN_OnlyINcen_Coortrans)



    # Load various genomic elements
    element_types = ['satellite', 'intactLTR', 'NUMT', 'NUPT', 'rDNA', 'gene']
    
    for element_type in element_types:
        element = f"{assembly}.{element_type}.bed"
        element_INcen = f"{assembly}.{element_type}.bed.INcen"
        element_INcen_Coortrans = f"{assembly}.{element_type}.bed.INcen.Coortrans"
    
        # Verify if the element file exists and is not empty
        if not os.path.exists(element) or os.path.getsize(element) == 0:
            # Create element_INcen_Coortrans (empty file or initialization)
            open(element_INcen_Coortrans, 'w').close()
            print(f"{element} does not exist or is empty; created {element_INcen_Coortrans}")
        else:
            # File exists and is non-empty; proceed with standard workflow
            bedtools_intersect(element, cenRegion_expand500k)
            merge_bed_in_place(element_INcen)
            convert_element_relative(cenRegion_expand500k, element_INcen, element_INcen_Coortrans)

def compare_bed_regions(bed1, bed2, output_bed):
    """
    Compare interval sizes for each chromosome between two BED files and 
    output the larger interval per chromosome to output_bed.
    
    Args:
        bed1: Path to the first BED file.
        bed2: Path to the second BED file.
        output_bed: Path to the output BED file.
    """
    
    def read_bed(file_path):
        regions = {}
        with open(file_path) as f:
            for line in f:
                if line.strip() == "":
                    continue
                chrom, start, end = line.strip().split()[:3]
                start = int(start)
                end = int(end)
                regions[chrom] = (start, end)
        return regions

    bed1_regions = read_bed(bed1)
    bed2_regions = read_bed(bed2)

    with open(output_bed, "w") as out:
        all_chroms = set(bed1_regions.keys()).union(set(bed2_regions.keys()))
        
        for chrom in sorted(all_chroms):
            r1 = bed1_regions.get(chrom)
            r2 = bed2_regions.get(chrom)

            if r1 and r2:
                len1 = r1[1] - r1[0]
                len2 = r2[1] - r2[0]
                chosen = r1 if len1 >= len2 else r2
            elif r1:
                chosen = r1
            elif r2:
                chosen = r2
            else:
                continue

            out.write(f"{chrom}\t{chosen[0]}\t{chosen[1]}\n")

def main(hap1, hap2):
    ### 1. Data Processing
    process_element_file(hap1)
    process_element_file(hap2)

    ### 2. Configure Sector Sizes for Plotting -----------------------------------------------------------
    # Load CENH3 data
    
    hap1_cenRegion_expand500k = f"{hap1}.cenRegion.expand500k.bed"
    hap2_cenRegion_expand500k = f"{hap2}.cenRegion.expand500k.bed"
    bigger_cenRegion_expand500k = f"bigger.cenRegion.expand500k.bed"
    compare_bed_regions(hap1_cenRegion_expand500k, hap2_cenRegion_expand500k, bigger_cenRegion_expand500k)


    ### 3. Visualization / Plotting ----------------------------------------------------------------------
    # Read centromere region file and define sectors
    sector = centromere_size_in_mb(bigger_cenRegion_expand500k)
    print(sector)
    circos = Circos(sector, space=2, start=15)
    for sector in circos.sectors:
        # Plot sector name
        sector.text(f"{sector.name}", r=110, size=15)

        # Plot CENH3 signals (three files per CENH3)
        cenh3_samples = ['sample1', 'sample2', 'sample3']
        location_types = ['INcen', 'OnlyINcen']
        location_type_colors = {'INcen':'#E7E7E7', 'OnlyINcen':'#ac1f18'}
        assemblys = [hap1, hap2]
        assembly_locations = {hap1:(95,100), hap2:(90,95)}
        for cenh3_sample in cenh3_samples:
            for location_type in location_types:
                for assembly in assemblys:
                    assembly_loc = assembly_locations[assembly]
                    cenh3_INcen_Coortrans = f"{assembly}.{cenh3_sample}.CENH3.bdg.{location_type}.Coortrans"
                    # Load CENH3 data
                    cenh3_df = pd.read_csv(cenh3_INcen_Coortrans, sep = '\t', header=None, names=['chr', 'start', 'end', 'signal'])
                    cenh3_df_chr = cenh3_df[cenh3_df['chr'] == sector.name]
                    cenh3_x = (cenh3_df_chr['start'] + cenh3_df_chr['end']) / 2
                    cenh3_x = np.clip(cenh3_x, 0, sector.size)
                    cenh3_y = cenh3_df_chr['signal'].values
                    cenh3_y = np.clip(cenh3_y, 0, 2)
                    cenh3_track = sector.add_track(assembly_loc, r_pad_ratio=0.1)
                    cenh3_track.axis()
                    # Select color
                    location_type_color = location_type_colors[location_type]
                    cenh3_track.line(cenh3_x, cenh3_y, color=location_type_color, lw=0.2, vmin=0, vmax=2)
                    if assembly == 'hap1':
                        cenh3_track.xticks_by_interval(interval=5, outer=True, label_size=6, label_formatter=lambda v: f"{v/10:.1f}M")
            print(f'***** Completed plotting CENH3 signal track for {assembly}: {sector.name} *****')



        # Plot sequencing tracks
        assemblys = [hap1, hap2]
        seqs = ['ONT.bed', 'hifi.bed', 'ONT.OnlyINcen', 'hifi.OnlyINcen']
        seq_locations = {'ONT.bed':{hap1:(56,61), hap2:(51,56)}, 'hifi.bed':{hap1:(43, 48), hap2:(38, 43)}, 'ONT.OnlyINcen':{hap1:(56,61), hap2:(51,56)}, 'hifi.OnlyINcen':{hap1:(43, 48), hap2:(38, 43)}}
        seq_colors = {'ONT.bed':'#E7E7E7', 'hifi.bed':'#E7E7E7', 'ONT.OnlyINcen':'#849AB8', 'hifi.OnlyINcen':'#B5979A'}
        for seq in seqs:
            for assembly in assemblys:
                seq_INcen_Coortrans = f"{assembly}.{seq}.INcen.Coortrans"
                seq_df = pd.read_csv(seq_INcen_Coortrans, sep = '\t', header=None, names=['chr', 'start', 'end', 'depth'])
                seq_df_chr = seq_df[seq_df['chr'] == sector.name]
                seq_x = (seq_df_chr['start'] + seq_df_chr['end']) / 2
                seq_x = np.clip(seq_x, 0, sector.size)
                seq_y = seq_df_chr['depth'].values
                seq_y = np.clip(seq_y, 0, 100)
                # Create track
                seq_location = seq_locations[seq][assembly]
                seq_color = seq_colors[seq]
                seq_track = sector.add_track(seq_location, r_pad_ratio=0.1)
                seq_track.axis()
                seq_track.fill_between(seq_x, seq_y, color=seq_color, vmin=0, vmax=100)



        # Plot various genomic element tracks
        assemblys = [hap1, hap2]
        element_types = ['satellite', 'intactLTR', 'NUMT', 'NUPT', 'rDNA', 'gene']
        element_locs = {'satellite':{hap1:(30, 35), hap2:(25, 30)}, 'intactLTR':{hap1:(30, 35), hap2:(25, 30)}, 'NUMT':{hap1:(30, 35), hap2:(25, 30)}, 'NUPT':{hap1:(30, 35), hap2:(25, 30)}, 'rDNA':{hap1:(30, 35), hap2:(25, 30)}, 'gene':{hap1:(30, 35), hap2:(25, 30)}}
        element_colors = {'satellite':'red', 'intactLTR':'#4485c7', 'NUMT':'#D4562E', 'NUPT':'#84ba42', 'rDNA':'#682487', 'gene':'#dbb428'}
        for element_type in element_types:
            for assembly in assemblys:
                element_INcen_Coortrans = f"{assembly}.{element_type}.bed.INcen.Coortrans"
                element_df = pd.read_csv(element_INcen_Coortrans, sep = '\t', header=None, names=['chr', 'start', 'end'])
                # Check if file is empty
                if element_df.empty:
                    print(f"File {element_INcen_Coortrans} is empty; skipping track generation.")
                    continue
                element_df_chr = element_df[element_df['chr'] == sector.name]
                # Create track
                location = element_locs[element_type][assembly]
                element_color = element_colors[element_type]
                element_track = sector.add_track(location, r_pad_ratio=0.05)
                element_track.axis()
                for idx, row in element_df_chr.iterrows():
                    start = row['start']
                    end = row['end']
                    element_track.rect(start, end, color=element_color)


        # Plot Hi-C tracks
        assemblys = [hap1, hap2]
        assembly_locations = {hap1:(82, 87), hap2:(77, 82)}
        hic_locations = ['INcen', 'OnlyINcen']
        hic_colors = {'INcen':'#E7E7E7', 'OnlyINcen':'#D4562E'}
        for assembly in assemblys:
            for hic in hic_locations:
                assembly_loc = assembly_locations[assembly]
                hic_INcen_Coortrans = f"{assembly}.hic.bed.{hic}.Coortrans"
                hic_df = pd.read_csv(hic_INcen_Coortrans, sep = '\t', header=None, names=['chr', 'start', 'end', 'depth'])
                hic_df_chr = hic_df[hic_df['chr'] == sector.name]
                hic_x = (hic_df_chr['start'] + hic_df_chr['end']) / 2
                hic_x = np.clip(hic_x, 0, sector.size)
                hic_y = hic_df_chr['depth'].values
                hic_y = np.clip(hic_y, 0, 20)
                # Create track
                hic_color = hic_colors[hic]
                hic_track = sector.add_track(assembly_loc, r_pad_ratio=0.1)
                hic_track.axis()
                hic_track.line(hic_x, hic_y, color=hic_color, vmin=0, vmax=20)
                hic_track.fill_between(hic_x, hic_y, color=hic_color, vmin=0, vmax=20)


        # Plot model-predicted centromere probabilities
        predictCEN_locations = ['INcen', 'OnlyINcen']
        predictCEN_colors = {'INcen':'#E7E7E7', 'OnlyINcen':'#682487'}
        assemblys = [hap1, hap2]
        assembly_locations = {hap1:(69, 74), hap2:(64, 69)}
        
        for predictCEN in predictCEN_locations:
            for assembly in assemblys:
                assembly_loc = assembly_locations[assembly]
                predictCEN_INcen_Coortrans = f"{assembly}.predictCEN.bed.{predictCEN}.Coortrans"
                predictCEN_df = pd.read_csv(predictCEN_INcen_Coortrans, sep='\t', header=None,
                                            names=['chr', 'start', 'end', 'depth'])
                predictCEN_df_chr = predictCEN_df[predictCEN_df['chr'] == sector.name]
            
                predictCEN_x = (predictCEN_df_chr['start'] + predictCEN_df_chr['end']) / 2
                predictCEN_x = np.clip(predictCEN_x, 0, sector.size)
            
                predictCEN_y = predictCEN_df_chr['depth'].values
                predictCEN_y = np.clip(predictCEN_y, 0, 20)
            
                # Create track
                predictCEN_color = predictCEN_colors[predictCEN]
                predictCEN_track = sector.add_track(assembly_loc, r_pad_ratio=0.1)
                predictCEN_track.axis()
                predictCEN_track.line(predictCEN_x, predictCEN_y, color=predictCEN_color, vmin=0, vmax=1)
                predictCEN_track.fill_between(predictCEN_x, predictCEN_y, color=predictCEN_color, vmin=0, vmax=1)


        # Add annotation labels
        circos.text('A', r=95, deg=0, size=10)
        circos.text('B', r=82, deg=0, size=10)
        circos.text('C', r=69, deg=0, size=10)
        circos.text('D', r=56, deg=0, size=10)
        circos.text('E', r=43, deg=0, size=10)
        circos.text('F', r=30, deg=0, size=10)



        print(f'***** Completed plotting H3K9me2 track for {assembly} *****')
        print(f"Completed all tracks for {sector.name}")
        print(f"-----------------------------------")


        # Plot satellite track
        
    output_pdf = f"{assembly}.multi_omic.pdf"
    circos.savefig(output_pdf)


main(args.hap1, args.hap2)