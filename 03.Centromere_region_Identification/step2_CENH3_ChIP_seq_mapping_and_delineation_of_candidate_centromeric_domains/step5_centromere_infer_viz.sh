/share/org/YZWL/yzwl_liubx/wangdan/pycharm_script/Oryza_Centromere_project/run_AllChrSegmentationCENH3_CBS-v6.R

$Rscript ~/wangdan/pycharm_script/Oryza_Centromere_project/run_AllChrSegmentationCENH3_CBS-v6.R
Usage: 
Options:
        --sample_file=FILE Segment file matched to original coordinates after the first round of Circular Binary Segmentation (CBS). Required.
        --raw_log2_file=FILE Raw log2-transformed ratio signal intensity values. Required.
        --known_cen_file=FILE File containing manually curated centromere coordinates. Required.
        --prefix=STRING  Prefix for the output directory and files. Required.
        --known_cen_filename_filter  Filter to specify which haplotype to read from the manually curated centromere file; must match the entries in the file. Required.
        --haplotype_suffix  Haplotype suffix used for plotting raw scatter plots. Default: hap1.
        --plot_segmented_points Whether to visualize segmented data points (intermediate files). Default: TRUE.
        --undo_sd_vec=VALUES Comma-separated list of undo.SD values. Default: 1.
        --cen_logratio_threshold  Log2-ratio threshold for centromere inference. Default: 0.5.
        --max_gap_points=INTEGER Maximum number of allowed gap points in low-signal regions during centromere inference. Default: 20.
        --min_peak_points=INTEGER Minimum number of consecutive data points (after gap filling) required for a region to be considered a valid peak for centromere inference. Default: 1.
        --chromosomes=CHR_LIST Comma-separated list of target chromosomes. Use 'ALL' to process all chromosomes found in the input file. Default: "ALL" (e.g., "Chr01,Chr02").
        --multi_cen_min_raw_block_points=INTEGER Enrichment parameter: Minimum number of consecutive data points in a raw block (all points exceeding cen_logratio_threshold) to be considered a candidate for multi-centromere detection. Set to 0 for single centromere mode; other values define the local window of scatter points to consider. Default: 50.
        --help
		
		
Execution Command:
Rscript ~/wangdan/pycharm_script/Oryza_Centromere_project/run_AllChrSegmentationCENH3_CBS-v6.R --sample_file GG_Omey_hap1.sample1.CENH3.segmentation.bdg --raw_log2_file ~/Oryza_Centromere/02.ChIP-seq_CENH3/02.CENH3_bdg_analysis/GG_Omey/GG_Omey_Sample1/GG_Omey_1_cenH3_log2ratio_2k.bdg --known_cen_file ~/Oryza_Centromere/tmmp/00.genome_for_centromere/00.Centromere\ regions.xlsx --prefix GG_Omey_hap1 --known_cen_filename_filter GG_Omey_hap1