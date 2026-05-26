import pandas as pd
import numpy as np
import os

# Configuration of file paths
hor_folder = r"D:\脚本\CEN155_HOR_all_species_results"
correct_hor_folder = r"D:\脚本\correct_HOR_results"
crm_file = r"D:\脚本\CRM_intersections_flanking_regions_full.tsv"
random_file = r"D:\脚本\random_1kb_intervals_no_CRM_full.tsv"

# Output directory configuration
output_dir = r"D:\脚本\HOR_score_results"
os.makedirs(output_dir, exist_ok=True)

crm_output = os.path.join(output_dir, "CRM_intervals_HOR_scores.tsv")
random_output = os.path.join(output_dir, "Random_intervals_HOR_scores.tsv")

print("=== 1. Data Acquisition: CRM and Random Interval Datasets ===")

# Load CRM dataset
crm_df = pd.read_csv(crm_file, sep='\t')
print(f"CRM dataset dimensions: {crm_df.shape}")

# Load random interval dataset
random_df = pd.read_csv(random_file, sep='\t')
print(f"Random interval dataset dimensions: {random_df.shape}")

print("\n=== 2. Interval Data Preparation ===")

# 2.1 Prepare CRM intervals (Consolidating Upstream and Downstream 1kb windows)
crm_upstream = crm_df[['sample_haplotype', 'chromosome', 'upstream_1kb_start', 'upstream_1kb_end']].copy()
crm_upstream['interval_type'] = 'upstream_1kb'
crm_upstream = crm_upstream.rename(columns={
    'sample_haplotype': 'haplotype',
    'chromosome': 'chrom',
    'upstream_1kb_start': 'interval_start',
    'upstream_1kb_end': 'interval_end'
})

crm_downstream = crm_df[['sample_haplotype', 'chromosome', 'downstream_1kb_start', 'downstream_1kb_end']].copy()
crm_downstream['interval_type'] = 'downstream_1kb'
crm_downstream = crm_downstream.rename(columns={
    'sample_haplotype': 'haplotype',
    'chromosome': 'chrom',
    'downstream_1kb_start': 'interval_start',
    'downstream_1kb_end': 'interval_end'
})

crm_intervals = pd.concat([crm_upstream, crm_downstream], ignore_index=True)
print(f"Total CRM intervals: {len(crm_intervals)}")

# 2.2 Prepare Random intervals
random_intervals = random_df[['sample', 'chrom', 'start', 'end']].copy()
random_intervals['interval_type'] = 'random_1kb'
random_intervals = random_intervals.rename(columns={
    'sample': 'haplotype',
    'start': 'interval_start',
    'end': 'interval_end'
})
print(f"Total random intervals: {len(random_intervals)}")

print("\n=== 3. Haplotype Mapping Configuration ===")

# Define mapping dictionary for haplotype standardization
haplotype_mapping = {
    'AA_Osat_hap1': 'AA_Osat_ind',
    'AA_Osat_hap2': 'AA_Osat_jap',
}

crm_intervals['hor_species'] = crm_intervals['haplotype'].map(lambda x: haplotype_mapping.get(x, x))
random_intervals['hor_species'] = random_intervals['haplotype'].map(lambda x: haplotype_mapping.get(x, x))

print(f"Mapped species in CRM intervals: {sorted(crm_intervals['hor_species'].unique())}")
print(f"Mapped species in random intervals: {sorted(random_intervals['hor_species'].unique())}")

print("\n=== 4. Loading HOR Metadata and Scores ===")

def load_hor_data():
    """Aggregates Higher-Order Repeat (HOR) data from multiple sources."""
    all_hor_dfs = []
    
    # 4.1 Load primary HOR data from specified results folder
    print("Parsing HOR results from the primary results directory...")
    
    # Identify species-specific subdirectories
    species_folders = [f for f in os.listdir(hor_folder) if os.path.isdir(os.path.join(hor_folder, f))]
    
    for species_folder in species_folders:
        species_path = os.path.join(hor_folder, species_folder)
        
        # Locate relevant CSV files within species directories
        csv_files = [f for f in os.listdir(species_path) if f.endswith('.csv') and 'CEN155_HOR_scores' in f]
        
        for csv_file in csv_files:
            file_path = os.path.join(species_path, csv_file)
            species_name = csv_file.replace('_CEN155_HOR_scores.csv', '')
            
            print(f"  Reading scores for {species_name}...")
            
            try:
                df = pd.read_csv(file_path, encoding='utf-8')
                
                # Select and standardize column headers
                df = df[['chromosome', 'physical_start', 'physical_end', 'HOR_score_percent']].copy()
                df = df.rename(columns={
                    'chromosome': 'chrom',
                    'physical_start': 'hor_position',
                    'physical_end': 'hor_end',
                    'HOR_score_percent': 'hor_score'
                })
                
                # Add species metadata
                df['species'] = species_name
                
                # Ensure numeric validity for coordinates and scores
                df['hor_position'] = pd.to_numeric(df['hor_position'], errors='coerce')
                df['hor_end'] = pd.to_numeric(df['hor_end'], errors='coerce')
                df['hor_score'] = pd.to_numeric(df['hor_score'], errors='coerce')
                
                # Drop records with invalid coordinates
                df = df.dropna(subset=['hor_position', 'hor_end'])
                
                all_hor_dfs.append(df)
                print(f"    Successfully loaded: {len(df)} records")
                
            except Exception as e:
                print(f"    Load failed: {e}")
    
    # 4.2 Load corrected HOR scores from specific Excel files
    print("\nParsing corrected HOR records...")
    
    correct_files = {
        'AA_Osat_ind': os.path.join(correct_hor_folder, "AA_Osat_ind_HOR_scores.xlsx"),
        'AA_Osat_jap': os.path.join(correct_hor_folder, "AA_Osat_jap_HOR_scores.xlsx")
    }
    
    for species, file_path in correct_files.items():
        if os.path.exists(file_path):
            print(f"  Reading corrected scores for {species}...")
            
            try:
                df = pd.read_excel(file_path)
                
                # Select and standardize column headers
                df = df[['chromosome', 'physical_start', 'physical_end', 'HOR_score_percent']].copy()
                df = df.rename(columns={
                    'chromosome': 'chrom',
                    'physical_start': 'hor_position',
                    'physical_end': 'hor_end',
                    'HOR_score_percent': 'hor_score'
                })
                
                df['species'] = species
                
                df['hor_position'] = pd.to_numeric(df['hor_position'], errors='coerce')
                df['hor_end'] = pd.to_numeric(df['hor_end'], errors='coerce')
                df['hor_score'] = pd.to_numeric(df['hor_score'], errors='coerce')
                
                df = df.dropna(subset=['hor_position', 'hor_end'])
                
                all_hor_dfs.append(df)
                print(f"    Successfully loaded: {len(df)} records")
                
            except Exception as e:
                print(f"    Load failed: {e}")
    
    # Aggregate all datasets
    if all_hor_dfs:
        all_hor_df = pd.concat(all_hor_dfs, ignore_index=True)
        print(f"\nAggregated HOR dataset: {len(all_hor_df)} total records")
        print(f"Unique species in dataset: {sorted(all_hor_df['species'].unique())}")
        return all_hor_df
    else:
        print("Error: No HOR data loaded successfully.")
        return None

# Execute data loading
all_hor_df = load_hor_data()
if all_hor_df is None:
    exit()

print("\n=== 5. Identifying Overlaps between Intervals and HOR Scores ===")

def find_overlaps(intervals_df, hor_df):
    """Determines spatial overlap between genomic intervals and HOR records."""
    results = []
    
    # Enforce integer types for coordinate logic
    intervals_df = intervals_df.copy()
    intervals_df['interval_start'] = intervals_df['interval_start'].astype(int)
    intervals_df['interval_end'] = intervals_df['interval_end'].astype(int)
    
    # Organize HOR data into nested dictionaries for efficient lookup (Species -> Chromosome)
    hor_dict = {}
    for species in hor_df['species'].unique():
        species_data = hor_df[hor_df['species'] == species]
        hor_dict[species] = {}
        
        for chrom in species_data['chrom'].unique():
            chrom_data = species_data[species_data['chrom'] == chrom]
            hor_dict[species][chrom] = chrom_data[['hor_position', 'hor_end', 'hor_score']].values.tolist()
    
    # Process each interval for overlap detection
    total_intervals = len(intervals_df)
    processed = 0
    
    for _, interval in intervals_df.iterrows():
        processed += 1
        if processed % 1000 == 0:
            print(f"  Processing progress: {processed}/{total_intervals}")
        
        haplotype = interval['haplotype']
        chrom = interval['chrom']
        interval_start = interval['interval_start']
        interval_end = interval['interval_end']
        interval_type = interval['interval_type']
        hor_species = interval['hor_species']
        
        # Check for matching species and chromosome metadata
        if hor_species in hor_dict and chrom in hor_dict[hor_species]:
            hor_entries = hor_dict[hor_species][chrom]
            overlaps = []
            
            # Evaluate each HOR record for intersection
            for hor_entry in hor_entries:
                hor_start = int(hor_entry[0])
                hor_end = int(hor_entry[1])
                hor_score = float(hor_entry[2])
                
                # Intersection validation logic
                if not (hor_end < interval_start or hor_start > interval_end):
                    overlaps.append({
                        'hor_position': hor_start,
                        'hor_score': hor_score
                    })
            
            if overlaps:
                # Record specific overlaps found
                for overlap in overlaps:
                    results.append({
                        'haplotype': haplotype,
                        'chrom': chrom,
                        'interval_start': interval_start,
                        'interval_end': interval_end,
                        'interval_type': interval_type,
                        'hor_position': overlap['hor_position'],
                        'hor_score': overlap['hor_score']
                    })
            else:
                # No overlap detected
                results.append({
                    'haplotype': haplotype,
                    'chrom': chrom,
                    'interval_start': interval_start,
                    'interval_end': interval_end,
                    'interval_type': interval_type,
                    'hor_position': 'no_overlap',
                    'hor_score': 'NA'
                })
        else:
            # Metadata mismatch or missing HOR data
            results.append({
                'haplotype': haplotype,
                'chrom': chrom,
                'interval_start': interval_start,
                'interval_end': interval_end,
                'interval_type': interval_type,
                'hor_position': 'no_overlap',
                'hor_score': 'NA'
            })
    
    return pd.DataFrame(results)

print("Calculating overlaps for CRM intervals...")
crm_results = find_overlaps(crm_intervals, all_hor_df)
print(f"CRM results processed: {len(crm_results)} records")

print("\nCalculating overlaps for random intervals...")
random_results = find_overlaps(random_intervals, all_hor_df)
print(f"Random results processed: {len(random_results)} records")

print("\n=== 6. Output Formatting and Storage ===")

# 6.1 Format CRM output headers
crm_final = crm_results.copy()
crm_final = crm_final.rename(columns={
    'haplotype': 'Species_Haplotype',
    'chrom': 'Chromosome',
    'interval_start': 'Interval_Start',
    'interval_end': 'Interval_End',
    'hor_position': 'HOR_Position',
    'hor_score': 'HOR_Score'
})

# Filter for relevant columns
crm_final = crm_final[['Species_Haplotype', 'Chromosome', 'Interval_Start', 'Interval_End', 'HOR_Position', 'HOR_Score']]

# 6.2 Format random output headers
random_final = random_results.copy()
random_final = random_final.rename(columns={
    'haplotype': 'Species_Haplotype',
    'chrom': 'Chromosome',
    'interval_start': 'Interval_Start',
    'interval_end': 'Interval_End',
    'hor_position': 'HOR_Position',
    'hor_score': 'HOR_Score'
})

# Filter for relevant columns
random_final = random_final[['Species_Haplotype', 'Chromosome', 'Interval_Start', 'Interval_End', 'HOR_Position', 'HOR_Score']]

# 6.3 Save results to disk
crm_final.to_csv(crm_output, sep='\t', index=False)
random_final.to_csv(random_output, sep='\t', index=False)

print(f"CRM interval results saved to: {crm_output}")
print(f"Random interval results saved to: {random_output}")

print("\n=== 7. Aggregate Statistics Summary ===")

crm_with_hor = crm_final[crm_final['HOR_Position'] != 'no_overlap']
random_with_hor = random_final[random_final['HOR_Position'] != 'no_overlap']

print(f"CRM Summary Statistics:")
print(f"  Total intervals analyzed: {len(crm_final)}")
print(f"  Intervals with HOR overlap: {len(crm_with_hor)} ({len(crm_with_hor)/len(crm_final)*100:.1f}%)")

print(f"\nRandom Background Summary Statistics:")
print(f"  Total intervals analyzed: {len(random_final)}")
print(f"  Intervals with HOR overlap: {len(random_with_hor)} ({len(random_with_hor)/len(random_final)*100:.1f}%)")

print("\n=== 8. Data Preview (Top 5 rows) ===")

print("\nCRM Result Preview:")
print(crm_final.head().to_string())
print("\nRandom Result Preview:")
print(random_final.head().to_string())

print("\nExecution finalized successfully!")