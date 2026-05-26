--- START OF FILE text/x-sh ---

# Navigate to the working directory
cd /share/org/YZWL/yzwl_shahd/wangdan/Oryza_Centromere/centromere_distance/mash_distance_noXX

# Step 1: Generate a Mash sketch file. 
# Comparative tests were performed using -s 10,000 and 100,000; results were highly similar.
/share/org/YZWL/yzwl_shahd/wangdan/software/mash-Linux64-v2.2/mash sketch -k 17 -s 10000 -o 02.mash_all_centromeres_k17_s10000 01.mash_input_fasta/*.fasta

# Step 2: Calculate pairwise distances using Mash.
/share/org/YZWL/yzwl_shahd/wangdan/software/mash-Linux64-v2.2/mash dist 02.mash_all_centromeres_k17_s10000.msh 02.mash_all_centromeres_k17_s10000.msh > 03.mash_distance_k17_s10000.tab

# Step 3: Convert the output tab-separated table into a square distance matrix.
python ../mash_distance/mash_to_matrix.py 03.mash_distance_k17_s10000.tab 04.mash_distance_k17_s10000.tsv

# Step 4: Perform dimensionality reduction and visualize results (MDS, etc.).
python plot_dimension_reduction_in_main.py 04.mash_distance_k17_s10000.tsv plot_config.tsv 05.mash_distance_k17_s10000_DR_plot --method mds --width 8 --height 2.2

# Step 5: Conduct PERMANOVA analysis utilizing the scikit-bio Python library.
python analyze_permanova_structure_addR2.py 04.mash_distance_k17_s10000.tsv 06.mash_distance_k17_s10000_Pseudo-F_R2 --width 3 --height 3

# Step 6: Perform PERMDISP (Permutational Analysis of Multivariate Dispersions) analysis.
# Rscript run_permdisp_comprehensive.R -i 04.mash_distance_k17_s10000.tsv -o 06.mash_distance_k17_s10000_permdisp