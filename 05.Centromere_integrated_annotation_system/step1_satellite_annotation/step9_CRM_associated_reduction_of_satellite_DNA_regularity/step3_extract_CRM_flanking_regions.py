import pandas as pd
import numpy as np

# Load datasets
print("Loading datasets...")
df_large_regions = pd.read_csv("/share/org/YZWL/yzwl_shahd/lxr/satellite_hor/new/CENH3_CEN155_merged_intersections.tsv", sep='\t')
df_crm = pd.read_csv("/share/org/YZWL/yzwl_shahd/lxr/satellite_hor/new/all_species_CRM_intervals.tsv", sep='\t')

print(f"Large regions file dimensions: {df_large_regions.shape}")
print(f"CRM file dimensions: {df_crm.shape}")

# Display sample metadata
print("\nSample types in large regions file (Top 10):")
print(df_large_regions['Sample'].unique()[:10])
print(f"Total unique samples in large regions file: {df_large_regions['Sample'].nunique()}")

print("\nSample types in CRM file (Top 10):")
print(df_crm['Sample'].unique()[:10])
print(f"Total unique samples in CRM file: {df_crm['Sample'].nunique()}")

# Rename columns for standardization
df_large_regions = df_large_regions.rename(columns={
    'Chromosome': 'chrom',
    'Merged_Intersection_Start': 'large_start',
    'Merged_Intersection_End': 'large_end'
})

df_crm = df_crm.rename(columns={
    'Chromosome': 'chrom',
    'Start': 'crm_start',
    'End': 'crm_end'
})

# Ensure numeric data types for coordinates
df_large_regions['large_start'] = pd.to_numeric(df_large_regions['large_start'], errors='coerce')
df_large_regions['large_end'] = pd.to_numeric(df_large_regions['large_end'], errors='coerce')
df_crm['crm_start'] = pd.to_numeric(df_crm['crm_start'], errors='coerce')
df_crm['crm_end'] = pd.to_numeric(df_crm['crm_end'], errors='coerce')

# Remove rows with missing coordinate values
df_large_regions = df_large_regions.dropna(subset=['large_start', 'large_end'])
df_crm = df_crm.dropna(subset=['crm_start', 'crm_end'])

print(f"\nCleaned large regions dimensions: {df_large_regions.shape}")
print(f"Cleaned CRM file dimensions: {df_crm.shape}")

# Initialize list for result storage (one result per intersection)
results = []

print("\nIdentifying intersections between CRM and large regions and calculating flanking windows...")
print(f"Total CRM records to process: {len(df_crm)}")

# Sample data verification for debugging purposes
print("\nDebugging Sample Data Check:")
print("CRM Samples:")
test_crm = df_crm.head(3).copy()
for idx, row in test_crm.iterrows():
    print(f"  {row['Sample']} {row['chrom']}: {row['crm_start']}-{row['crm_end']}")

print("\nCorresponding Large Region Matches:")
for idx, crm_row in test_crm.iterrows():
    chrom = crm_row['chrom']
    sample = crm_row['Sample']
    
    # Locate matching large regions based on chromosome and sample
    matching = df_large_regions[
        (df_large_regions['chrom'] == chrom) & 
        (df_large_regions['Sample'] == sample)
    ]
    
    if len(matching) > 0:
        for _, large_row in matching.iterrows():
            print(f"  Matched Large Region: {large_row['Sample']} {large_row['chrom']}: {large_row['large_start']}-{large_row['large_end']}")
    else:
        print(f"  No matching large region found for: {sample} {chrom}")

# Track intersection counts per species
species_counts = {}

# Iterate through each CRM to evaluate intersections with large regions
processed = 0
found_intersections = 0

for idx, crm_row in df_crm.iterrows():
    chrom = crm_row['chrom']
    crm_sample = crm_row['Sample']
    crm_start = crm_row['crm_start']
    crm_end = crm_row['crm_end']
    
    # Extract species info (assuming 'Species_Haplotype' format)
    species = '_'.join(crm_sample.split('_')[:2]) if '_' in crm_sample else crm_sample
    
    # Find matching large regions on the same sample and chromosome
    matching_large_regions = df_large_regions[
        (df_large_regions['chrom'] == chrom) & 
        (df_large_regions['Sample'] == crm_sample)
    ]
    
    # Fallback: Match by chromosome only if sample-specific match fails
    if len(matching_large_regions) == 0:
        matching_large_regions = df_large_regions[df_large_regions['chrom'] == chrom]
    
    intersection_found = False
    
    for _, large_row in matching_large_regions.iterrows():
        large_start = large_row['large_start']
        large_end = large_row['large_end']
        
        # Calculate overlap coordinates
        overlap_start = max(crm_start, large_start)
        overlap_end = min(crm_end, large_end)
        
        # Validate that a physical overlap exists
        if overlap_start < overlap_end:
            # Ensure the intersection is a subset of the current CRM
            if overlap_start >= crm_start and overlap_end <= crm_end:
                
                # Define flanking regions (1kb upstream/downstream)
                upstream_start = max(0, overlap_start - 1000)
                upstream_end = overlap_start
                
                downstream_start = overlap_end
                downstream_end = overlap_end + 1000
                
                # Construct result record
                result_row = {
                    'sample_haplotype': crm_sample,
                    'chromosome': chrom,
                    
                    # CRM intersection coordinates
                    'crm_intersect_start': int(overlap_start),
                    'crm_intersect_end': int(overlap_end),
                    
                    # Upstream 1kb flanking window
                    'upstream_1kb_start': int(upstream_start),
                    'upstream_1kb_end': int(upstream_end),
                    
                    # Downstream 1kb flanking window
                    'downstream_1kb_start': int(downstream_start),
                    'downstream_1kb_end': int(downstream_end),
                    
                    # Metadata and context
                    'original_crm_start': int(crm_start),
                    'original_crm_end': int(crm_end),
                    'large_region_sample': large_row['Sample'],
                    'large_region_start': int(large_start),
                    'large_region_end': int(large_end),
                    'intersection_length': int(overlap_end - overlap_start),
                    
                    # Species metadata
                    'species': species
                }
                
                results.append(result_row)
                found_intersections += 1
                intersection_found = True
                
                # Increment species-specific counter
                if species not in species_counts:
                    species_counts[species] = 0
                species_counts[species] += 1
                
                break  # Process only the first valid intersection per CRM
    
    processed += 1
    if processed % 1000 == 0:
        print(f"Processed {processed}/{len(df_crm)} CRMs ({processed/len(df_crm)*100:.1f}%)")
        print(f"  Total intersections identified: {found_intersections}")

# Conversion to DataFrame and result export
if results:
    results_df = pd.DataFrame(results)
    
    print(f"\nTotal CRM intersections identified: {len(results_df)}")
    
    # Define required column order for output
    required_columns = [
        'sample_haplotype',      # Col 1: Species/Haplotype
        'chromosome',            # Col 2: Chromosome
        'crm_intersect_start',   # Col 3: Intersection Start
        'crm_intersect_end',     # Col 4: Intersection End
        'upstream_1kb_start',    # Col 5: Upstream Start
        'upstream_1kb_end',      # Col 6: Upstream End
        'downstream_1kb_start',  # Col 7: Downstream Start
        'downstream_1kb_end'     # Col 8: Downstream End
    ]
    
    # Supplemental metadata columns
    optional_columns = [
        'original_crm_start', 'original_crm_end',
        'large_region_sample', 'large_region_start', 'large_region_end',
        'intersection_length', 'species'
    ]
    
    # Generate final standardized DataFrame
    final_df = results_df[required_columns]
    
    print("\nFormat Preview of Final Output:")
    print(final_df.head(10))
    
    # Export primary results
    main_output_path = "/share/org/YZWL/yzwl_shahd/lxr/satellite_hor/new/CRM_intersections_flanking_regions.tsv"
    final_df.to_csv(main_output_path, sep='\t', index=False)
    print(f"\nPrimary results saved to: {main_output_path}")
    
    # Export full dataset including metadata
    full_df = pd.concat([final_df, results_df[optional_columns]], axis=1)
    full_output_path = "/share/org/YZWL/yzwl_shahd/lxr/satellite_hor/new/CRM_intersections_flanking_regions_full.tsv"
    full_df.to_csv(full_output_path, sep='\t', index=False)
    print(f"Full results saved to: {full_output_path}")
    
    # Summary Statistics
    print("\n" + "="*80)
    print("Aggregate Statistics:")
    print(f"Total CRM input: {len(df_crm):,}")
    print(f"Total CRM intersections: {len(results_df):,}")
    print(f"Intersection ratio: {len(results_df)/len(df_crm)*100:.2f}%")
    
    # Species-specific statistics
    print("\n" + "="*80)
    print("Species-level Statistics Summary:")
    print("-"*80)
    
    species_stats = []
    
    for species, count in species_counts.items():
        # Calculate total CRMs for the specific species
        species_crms = df_crm[df_crm['Sample'].str.contains(species)].shape[0]
        percentage = (count / species_crms * 100) if species_crms > 0 else 0
        
        species_stats.append({
            'Species': species,
            'Total_CRMs': species_crms,
            'Intersection_CRMs': count,
            'Percentage': f"{percentage:.2f}%"
        })
        
        print(f"{species}:")
        print(f"  Total CRMs: {species_crms:,}")
        print(f"  Intersections found: {count:,}")
        print(f"  Intersection percentage: {percentage:.2f}%")
        print()
    
    # Save species statistics report
    species_df = pd.DataFrame(species_stats)
    species_output = "/share/org/YZWL/yzwl_shahd/lxr/satellite_hor/new/species_intersection_stats.tsv"
    species_df.to_csv(species_output, sep='\t', index=False)
    print(f"Species statistics report saved to: {species_output}")
    
    if len(results_df) > 0:
        # Intersection length metrics
        print("\n" + "="*80)
        print("Intersection Length Analysis:")
        print(f"  Minimum length: {results_df['intersection_length'].min():,} bp")
        print(f"  Maximum length: {results_df['intersection_length'].max():,} bp")
        print(f"  Mean length: {results_df['intersection_length'].mean():,.0f} bp")
        print(f"  Median length: {results_df['intersection_length'].median():,.0f} bp")
        
        # Distribution by chromosome
        print("\n" + "="*80)
        print("Distribution by Chromosome (Top 20):")
        chrom_counts = final_df['chromosome'].value_counts().head(20)
        for chrom, count in chrom_counts.items():
            print(f"  {chrom}: {count:,} intersections")
        
        # Distribution by sample
        print("\n" + "="*80)
        print("Distribution by Sample (Top 20):")
        sample_counts = final_df['sample_haplotype'].value_counts().head(20)
        for sample, count in sample_counts.items():
            print(f"  {sample}: {count:,} intersections")
    
    # BED file generation for downstream tools
    print("\n" + "="*80)
    print("Generating BED format files for visualization...")
    
    # Upstream BED
    upstream_bed = final_df[['chromosome', 'upstream_1kb_start', 'upstream_1kb_end']].copy()
    upstream_bed.columns = ['chrom', 'start', 'end']
    upstream_bed['name'] = final_df.apply(
        lambda row: f"{row['sample_haplotype']}_upstream_{row['crm_intersect_start']}-{row['crm_intersect_end']}", 
        axis=1
    )
    upstream_bed['score'] = 0
    upstream_bed['strand'] = '+'
    upstream_bed = upstream_bed[['chrom', 'start', 'end', 'name', 'score', 'strand']]
    
    upstream_bed_output = "/share/org/YZWL/yzwl_shahd/lxr/satellite_hor/new/CRM_upstream_1kb_regions.bed"
    upstream_bed.to_csv(upstream_bed_output, sep='\t', header=False, index=False)
    print(f"Upstream BED saved to: {upstream_bed_output}")
    
    # Downstream BED
    downstream_bed = final_df[['chromosome', 'downstream_1kb_start', 'downstream_1kb_end']].copy()
    downstream_bed.columns = ['chrom', 'start', 'end']
    downstream_bed['name'] = final_df.apply(
        lambda row: f"{row['sample_haplotype']}_downstream_{row['crm_intersect_start']}-{row['crm_intersect_end']}", 
        axis=1
    )
    downstream_bed['score'] = 0
    downstream_bed['strand'] = '+'
    downstream_bed = downstream_bed[['chrom', 'start', 'end', 'name', 'score', 'strand']]
    
    downstream_bed_output = "/share/org/YZWL/yzwl_shahd/lxr/satellite_hor/new/CRM_downstream_1kb_regions.bed"
    downstream_bed.to_csv(downstream_bed_output, sep='\t', header=False, index=False)
    print(f"Downstream BED saved to: {downstream_bed_output}")
    
    # Intersection BED
    intersect_bed = final_df[['chromosome', 'crm_intersect_start', 'crm_intersect_end']].copy()
    intersect_bed.columns = ['chrom', 'start', 'end']
    intersect_bed['name'] = final_df.apply(
        lambda row: f"{row['sample_haplotype']}_intersect_{row['crm_intersect_start']}-{row['crm_intersect_end']}", 
        axis=1
    )
    intersect_bed['score'] = 1000
    intersect_bed['strand'] = '+'
    intersect_bed = intersect_bed[['chrom', 'start', 'end', 'name', 'score', 'strand']]
    
    intersect_bed_output = "/share/org/YZWL/yzwl_shahd/lxr/satellite_hor/new/CRM_intersection_regions.bed"
    intersect_bed.to_csv(intersect_bed_output, sep='\t', header=False, index=False)
    print(f"Intersection BED saved to: {intersect_bed_output}")
    
    # Display final format example
    print("\n" + "="*80)
    print("Example Output for Publication/Analysis:")
    print("-"*80)
    print("sample_haplotype\tchromosome\tcrm_intersect_start\tcrm_intersect_end\tupstream_1kb_start\tupstream_1kb_end\tdownstream_1kb_start\tdownstream_1kb_end")
    print("-"*80)
    
    for i in range(min(5, len(final_df))):
        row = final_df.iloc[i]
        print(f"{row['sample_haplotype']}\t{row['chromosome']}\t{row['crm_intersect_start']}\t{row['crm_intersect_end']}\t{row['upstream_1kb_start']}\t{row['upstream_1kb_end']}\t{row['downstream_1kb_start']}\t{row['downstream_1kb_end']}")
    
    print("="*80)
    
else:
    print("No intersections were identified.")

print("\nProcessing complete.")