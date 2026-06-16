import pandas as pd
import numpy as np
import random

# Set random seeds to ensure reproducibility of results
random.seed(42)
np.random.seed(42)

# Data Acquisition
print("Loading datasets...")
df_large_regions = pd.read_csv("/share/org/YZWL/yzwl_shahd/lxr/satellite_hor/new/CENH3_CEN155_merged_intersections.tsv", sep='\t')
df_crm = pd.read_csv("/share/org/YZWL/yzwl_shahd/lxr/satellite_hor/new/all_species_CRM_intervals.tsv", sep='\t')

print(f"Large region file dimensions: {df_large_regions.shape}")
print(f"CRM interval file dimensions: {df_crm.shape}")

# Column Renaming for consistency
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

# Ensure coordinate data are numeric types
df_large_regions['large_start'] = pd.to_numeric(df_large_regions['large_start'], errors='coerce')
df_large_regions['large_end'] = pd.to_numeric(df_large_regions['large_end'], errors='coerce')
df_crm['crm_start'] = pd.to_numeric(df_crm['crm_start'], errors='coerce')
df_crm['crm_end'] = pd.to_numeric(df_crm['crm_end'], errors='coerce')

# Remove records with missing coordinate values
df_large_regions = df_large_regions.dropna(subset=['large_start', 'large_end'])
df_crm = df_crm.dropna(subset=['crm_start', 'crm_end'])

print(f"\nDimensions after cleaning (large regions): {df_large_regions.shape}")
print(f"Dimensions after cleaning (CRM intervals): {df_crm.shape}")

# Step 1: Identify all intersections between CRMs and large genomic regions, including flanking windows
print("\nStep 1: Identifying intersections between CRMs and target large regions...")

def find_all_intersections(df_large, df_crm):
    """Identifies all genomic intersections between CRM and defined large regions."""
    intersections = []
    
    for _, crm_row in df_crm.iterrows():
        chrom = crm_row['chrom']
        crm_sample = crm_row['Sample']
        crm_start = crm_row['crm_start']
        crm_end = crm_row['crm_end']
        
        # Extract species and haplotype metadata
        # Assuming sample format: AA_Osat_hap1, AA_Ogla_hap2, etc.
        species_parts = crm_sample.split('_')
        if len(species_parts) >= 3:
            species = f"{species_parts[0]}_{species_parts[1]}"  # e.g., AA_Osat
            haplotype = crm_sample  # Retain full haplotype designation
        else:
            species = crm_sample
            haplotype = crm_sample
        
        # Locate corresponding large regions based on sample and chromosome
        matching_large = df_large[
            (df_large['chrom'] == chrom) & 
            (df_large['Sample'] == crm_sample)
        ]
        
        if len(matching_large) == 0:
            matching_large = df_large[df_large['chrom'] == chrom]
        
        for _, large_row in matching_large.iterrows():
            large_start = large_row['large_start']
            large_end = large_row['large_end']
            
            # Calculate the intersection coordinates
            overlap_start = max(crm_start, large_start)
            overlap_end = min(crm_end, large_end)
            
            if overlap_start < overlap_end:  # Physical overlap confirmed
                intersections.append({
                    'species': species,
                    'haplotype': haplotype,
                    'sample': crm_sample,
                    'chrom': chrom,
                    'crm_start': crm_start,
                    'crm_end': crm_end,
                    'large_sample': large_row['Sample'],
                    'large_start': large_start,
                    'large_end': large_end,
                    'intersect_start': overlap_start,
                    'intersect_end': overlap_end,
                    # Define 1kb flanking regions
                    'upstream_start': max(0, overlap_start - 1000),
                    'upstream_end': overlap_start,
                    'downstream_start': overlap_end,
                    'downstream_end': overlap_end + 1000
                })
                break
    
    return pd.DataFrame(intersections)

# Execute Step 1
intersections_df = find_all_intersections(df_large_regions, df_crm)
print(f"Identified {len(intersections_df)} intersections between CRM intervals and large genomic regions.")

if len(intersections_df) == 0:
    print("No intersections detected; terminating subtraction process.")
    exit()

# Step 2: Subtract CRM-associated regions from large genomic regions to define available background segments
print("\nStep 2: Subtracting CRM-associated regions from large genomic regions...")

def subtract_regions_from_large(large_start, large_end, regions_to_remove):
    """
    Subtracts a list of specified intervals from a parent genomic region.
    
    Parameters:
    - large_start, large_end: Boundaries of the parent region.
    - regions_to_remove: List of tuples (start, end) representing regions to be excluded.
    
    Returns:
    - List of remaining genomic segments.
    """
    # Sort exclusion regions by start position
    regions_to_remove = sorted(regions_to_remove, key=lambda x: x[0])
    
    remaining_regions = []
    current_start = large_start
    
    for remove_start, remove_end in regions_to_remove:
        # Skip exclusion regions outside the parent boundary (upstream)
        if remove_end <= current_start:
            continue
        
        # Terminate if exclusion regions are outside parent boundary (downstream)
        if remove_start >= large_end:
            break
        
        # Standardize exclusion boundaries to fit within the parent region
        remove_start = max(remove_start, current_start)
        remove_end = min(remove_end, large_end)
        
        # If there is a gap before the exclusion region, append as a remaining segment
        if remove_start > current_start:
            remaining_regions.append((current_start, remove_start))
        
        # Advance current position beyond the excluded segment
        current_start = max(current_start, remove_end)
    
    # Process the final remaining segment
    if current_start < large_end:
        remaining_regions.append((current_start, large_end))
    
    return remaining_regions

# Store all available remaining segments across genomic regions
remaining_regions_all = []

# Process each large region by sample and chromosome
for _, large_row in df_large_regions.iterrows():
    chrom = large_row['chrom']
    sample = large_row['Sample']
    large_start = large_row['large_start']
    large_end = large_row['large_end']
    
    # Metadata extraction
    species_parts = sample.split('_')
    if len(species_parts) >= 3:
        species = f"{species_parts[0]}_{species_parts[1]}"
        haplotype = sample
    else:
        species = sample
        haplotype = sample
    
    # Identify CRM intersections associated with this specific large region
    relevant_intersections = intersections_df[
        (intersections_df['chrom'] == chrom) & 
        ((intersections_df['sample'] == sample) | (intersections_df['large_sample'] == sample))
    ]
    
    if len(relevant_intersections) == 0:
        # If no CRM intersections exist, the entire large region is available
        remaining_regions_all.append({
            'species': species,
            'haplotype': haplotype,
            'sample': sample,
            'chrom': chrom,
            'region_start': large_start,
            'region_end': large_end,
            'region_length': large_end - large_start,
            'original_large_start': large_start,
            'original_large_end': large_end,
            'contains_crm': 0
        })
        continue
    
    # Compile exclusion list: intersection core and 1kb flanking windows
    regions_to_remove = []
    
    for _, intersect_row in relevant_intersections.iterrows():
        # Core intersection
        regions_to_remove.append((intersect_row['intersect_start'], intersect_row['intersect_end']))
        # 1kb flanking segments
        regions_to_remove.append((intersect_row['upstream_start'], intersect_row['upstream_end']))
        regions_to_remove.append((intersect_row['downstream_start'], intersect_row['downstream_end']))
    
    # Merge overlapping or adjacent exclusion intervals
    regions_to_remove = sorted(regions_to_remove, key=lambda x: x[0])
    merged_regions = []
    
    for region in regions_to_remove:
        if not merged_regions:
            merged_regions.append(list(region))
        else:
            last_region = merged_regions[-1]
            if region[0] <= last_region[1]:
                last_region[1] = max(last_region[1], region[1])
            else:
                merged_regions.append(list(region))
    
    # Execute subtraction from the parent large region
    remaining_regions = subtract_regions_from_large(large_start, large_end, merged_regions)
    
    # Log valid background segments
    for rem_start, rem_end in remaining_regions:
        remaining_regions_all.append({
            'species': species,
            'haplotype': haplotype,
            'sample': sample,
            'chrom': chrom,
            'region_start': rem_start,
            'region_end': rem_end,
            'region_length': rem_end - rem_start,
            'original_large_start': large_start,
            'original_large_end': large_end,
            'contains_crm': 1
        })

# Conversion to DataFrame
remaining_df = pd.DataFrame(remaining_regions_all)

print(f"\nGenerated {len(remaining_df)} discrete background genomic segments.")

if len(remaining_df) == 0:
    print("No available background segments found for random sampling.")
    exit()

# Export background segments for reference
remaining_output = "/share/org/YZWL/yzwl_shahd/lxr/satellite_hor/new/CENH3_regions_without_CRM_fragments.tsv"
remaining_df.to_csv(remaining_output, sep='\t', index=False)
print(f"Background segments metadata saved to: {remaining_output}")

# Step 3: Generate all possible 1kb candidate intervals per haplotype
print("\n" + "="*80)
print("Step 3: Generating candidate 1kb intervals for each haplotype...")

def generate_all_possible_1kb_intervals_for_haplotype(remaining_df, haplotype, interval_length=1000):
    """
    Generates all possible sliding/discrete 1kb intervals for a specific haplotype.
    """
    haplotype_remaining = remaining_df[remaining_df['haplotype'] == haplotype].copy()
    
    if len(haplotype_remaining) == 0:
        print(f"  Haplotype {haplotype}: No background segments available.")
        return [], {
            'haplotype': haplotype,
            'total_regions': 0,
            'possible_intervals': 0,
            'total_length': 0
        }
    
    # Filter for segments with sufficient length (>= 1000bp)
    valid_regions = haplotype_remaining[haplotype_remaining['region_length'] >= interval_length].copy()
    
    if len(valid_regions) == 0:
        print(f"  Haplotype {haplotype}: No background segments of sufficient length (>=1kb) found.")
        return [], {
            'haplotype': haplotype,
            'total_regions': len(haplotype_remaining),
            'possible_intervals': 0,
            'total_length': haplotype_remaining['region_length'].sum()
        }
    
    # Enumerate all sliding window intervals of 1kb
    all_intervals = []
    
    for idx, region in valid_regions.iterrows():
        region_start = int(region['region_start'])
        region_end = int(region['region_end'])
        max_start = region_end - interval_length
        
        for start in range(region_start, max_start + 1):
            end = start + interval_length
            if end <= region_end:
                all_intervals.append({
                    'haplotype': haplotype,
                    'species': region['species'],
                    'sample': region['sample'],
                    'chrom': region['chrom'],
                    'start': start,
                    'end': end,
                    'length': interval_length,
                    'source_region_start': region_start,
                    'source_region_end': region_end
                })
    
    stats = {
        'haplotype': haplotype,
        'species': valid_regions['species'].iloc[0] if len(valid_regions) > 0 else '',
        'total_regions': len(haplotype_remaining),
        'valid_regions_ge_1kb': len(valid_regions),
        'total_length': haplotype_remaining['region_length'].sum(),
        'possible_intervals': len(all_intervals)
    }
    
    return all_intervals, stats

# Step 4: Random sampling of 1000 intervals (1kb) per haplotype
print("\n" + "="*80)
print("Step 4: Randomly sampling 1000 background intervals (1kb) per haplotype...")

all_haplotypes = sorted(remaining_df['haplotype'].unique())
print(f"\nProcessing {len(all_haplotypes)} unique haplotypes:")

all_random_intervals = []
all_stats = []

for haplotype in all_haplotypes:
    print(f"\n{'='*60}")
    print(f"Processing Haplotype: {haplotype}")
    
    # Generate population of possible intervals
    all_intervals, stats = generate_all_possible_1kb_intervals_for_haplotype(
        remaining_df, haplotype, interval_length=1000
    )
    
    print(f"  Total discrete segments: {stats['total_regions']}")
    print(f"  Segments >= 1kb: {stats['valid_regions_ge_1kb']}")
    print(f"  Aggregate available length: {stats['total_length']:,} bp")
    print(f"  Total candidate 1kb intervals: {stats['possible_intervals']}")
    
    target_num = 1000
    
    if len(all_intervals) == 0:
        print(f"  Failure: No 1kb intervals available for {haplotype}")
        stats.update({
            'requested': target_num,
            'selected': 0,
            'selection_rate': 0,
            'status': 'Failed: No possible intervals'
        })
    elif len(all_intervals) < target_num:
        print(f"  Warning: Only {len(all_intervals)} intervals available (requested {target_num})")
        selected_intervals = all_intervals.copy()
        random.shuffle(selected_intervals)
        
        stats.update({
            'requested': target_num,
            'selected': len(selected_intervals),
            'selection_rate': len(selected_intervals) / target_num * 100,
            'status': f'Partial Success: Only {len(selected_intervals)} available'
        })
    else:
        # Perform random sampling without replacement
        selected_intervals = random.sample(all_intervals, target_num)
        
        stats.update({
            'requested': target_num,
            'selected': target_num,
            'selection_rate': 100,
            'status': 'Success'
        })
    
    if 'selected_intervals' in locals() and selected_intervals:
        all_random_intervals.extend(selected_intervals)
        print(f"  Successfully sampled {len(selected_intervals)} intervals.")
    
    all_stats.append(stats)

# Result Consolidation
random_intervals_df = pd.DataFrame(all_random_intervals)
stats_df = pd.DataFrame(all_stats)

if len(random_intervals_df) > 0:
    print(f"\n{'='*80}")
    print(f"Aggregated Total: {len(random_intervals_df)} random 1kb intervals sampled.")
    
    # Validation of interval lengths
    random_intervals_df['actual_length'] = random_intervals_df['end'] - random_intervals_df['start']
    valid_count = (random_intervals_df['actual_length'] == 1000).sum()
    print(f"  Validation: {valid_count}/{len(random_intervals_df)} intervals match 1000bp length exactly.")
    
    # Detailed Sampling Statistics
    print("\nSampling Summary Report:")
    print("-" * 120)
    print("Haplotype\t\tSpecies\tRegions\t>=1kb Regions\tTotal Length\tCandidate Pools\tReq\tSelected\tRate\tStatus")
    print("-" * 120)
    for stats in all_stats:
        print(f"{stats['haplotype']}\t{stats.get('species', '')}\t{stats['total_regions']}\t\t"
              f"{stats['valid_regions_ge_1kb']}\t\t{stats['total_length']:,}\t"
              f"{stats['possible_intervals']}\t\t{stats['requested']}\t"
              f"{stats['selected']}\t{stats['selection_rate']:.1f}%\t{stats['status']}")
    
    # Data Export
    print("\n" + "="*80)
    print("Exporting results to files...")
    
    # Primary output: sampled 1kb intervals
    random_output = "/share/org/YZWL/yzwl_shahd/lxr/satellite_hor/new/random_1kb_intervals_no_CRM_by_haplotype.tsv"
    random_intervals_df[['haplotype', 'species', 'sample', 'chrom', 'start', 'end', 'length']].to_csv(
        random_output, sep='\t', index=False
    )
    print(f"Sampling results saved to: {random_output}")
    
    # Full metadata output
    random_full_output = "/share/org/YZWL/yzwl_shahd/lxr/satellite_hor/new/random_1kb_intervals_no_CRM_by_haplotype_full.tsv"
    random_intervals_df.to_csv(random_full_output, sep='\t', index=False)
    
    # BED Format Generation
    print("\nGenerating BED files for visualization...")
    bed_df = random_intervals_df[['chrom', 'start', 'end']].copy()
    bed_df['name'] = random_intervals_df.apply(
        lambda row: f"{row['haplotype']}_random_{row['start']}-{row['end']}", 
        axis=1
    )
    bed_df['score'] = 0
    bed_df['strand'] = '+'
    bed_df = bed_df[['chrom', 'start', 'end', 'name', 'score', 'strand']]
    
    bed_output = "/share/org/YZWL/yzwl_shahd/lxr/satellite_hor/new/random_1kb_intervals_no_CRM_by_haplotype.bed"
    bed_df.to_csv(bed_output, sep='\t', header=False, index=False)
    
    # Generate separate BED files per haplotype
    for haplotype in random_intervals_df['haplotype'].unique():
        haplotype_df = random_intervals_df[random_intervals_df['haplotype'] == haplotype]
        if len(haplotype_df) > 0:
            haplotype_bed = haplotype_df[['chrom', 'start', 'end']].copy()
            haplotype_bed['name'] = haplotype_df.apply(
                lambda row: f"random_{row['start']}-{row['end']}", 
                axis=1
            )
            haplotype_bed['score'] = 0
            haplotype_bed['strand'] = '+'
            haplotype_bed = haplotype_bed[['chrom', 'start', 'end', 'name', 'score', 'strand']]
            
            haplotype_simple = haplotype.replace('_', '')
            haplotype_bed_output = f"/share/org/YZWL/yzwl_shahd/lxr/satellite_hor/new/random_1kb_intervals_no_CRM_{haplotype_simple}.bed"
            haplotype_bed.to_csv(haplotype_bed_output, sep='\t', header=False, index=False)

    # Performance Analytics
    total_original_length = (df_large_regions['large_end'] - df_large_regions['large_start']).sum()
    total_remaining_length = remaining_df['region_length'].sum()
    total_removed_length = total_original_length - total_remaining_length
    removal_percentage = (total_removed_length / total_original_length * 100) if total_original_length > 0 else 0
    
    print(f"\nGenomic Region Subtraction Analysis:")
    print(f"Aggregate original target length: {total_original_length:,} bp")
    print(f"Aggregate remaining background length: {total_remaining_length:,} bp")
    print(f"Aggregate excluded (CRM-associated) length: {total_removed_length:,} bp")
    print(f"Exclusion percentage: {removal_percentage:.2f}%")

else:
    print("\nSampling failed: No background intervals were generated.")

print("\nPipeline execution complete.")