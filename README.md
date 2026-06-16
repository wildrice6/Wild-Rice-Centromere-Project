Wild‑Rice‑Centromere‑Project
Oryza Genus Centromere Comparative Genomics Analysis Pipeline
This repository contains a complete set of bioinformatics analysis scripts for centromere‑focused research across 26 Oryza species. It supports the characterization of centromeric satellite repeats, high‑order repeat (HOR) structures, transposable elements, epigenetic associations and evolutionary dynamics based on T2T gap‑free genome assemblies, specifically designed for manuscript‑scale data analysis.

Project Overview
All pipelines are dedicated to comparative genomic investigation of centromeres within the genus Oryza. Major analytical modules include:
- Centromere boundary identification and integration with CENH3 ChIP‑seq data
- Discovery and classification of centromeric satellite monomers (CEN126 / CEN155)
- Detection, scoring and structural characterization of high‑order repeat (HOR) units
- Enrichment analysis of CRM‑LTR retrotransposons and their interplay with centromeric regions
- Evaluation of monomer similarity, AT‑content, copy number and centromeric coverage
- Species‑wide comparison of HOR block distribution and adjacent‑block distance dynamics
- Pan‑centromere gene density profiling and functional enrichment analysis
- Batch statistical computing and visualization for publication‑ready figures

Full Pipeline Directory Structure
- 01.Quality_assessment (https://github.com/wildrice6/Wild-Rice-Centromere-Project/tree/main/01.Quality_assessment) — Genome quality evaluation and sequence quality control
- 02.SV_between_haplotypes (https://github.com/wildrice6/Wild-Rice-Centromere-Project/tree/main/02.SV_between_haplotypes) — Structural variation analysis between haplotypes
- 03.Centromere_region_Identification (https://github.com/wildrice6/Wild-Rice-Centromere-Project/tree/main/03.Centromere_region_Identification) — Centromere region localization and boundary annotation
- 04.Centromere_position_normalization_and_haplotype_similarity_analysis (https://github.com/wildrice6/Wild-Rice-Centromere-Project/tree/main/04.Centromere_position_normalization_and_haplotype_similarity_analysis) — Centromere position normalization and haplotype similarity comparison
- 05.Centromere_integrated_annotation_system (https://github.com/wildrice6/Wild-Rice-Centromere-Project/tree/main/05.Centromere_integrated_annotation_system) — Integrated annotation workflow for centromeric regions
- 06.Centromeric_sequences_and_evolutionary_analysis (https://github.com/wildrice6/Wild-Rice-Centromere-Project/tree/main/06.Centromeric_sequences_and_evolutionary_analysis) — Centromeric sequence mining and evolutionary dynamics analysis
- 07.Repeat_monomer_and_HOR_inference (https://github.com/wildrice6/Wild-Rice-Centromere-Project/tree/main/07.Repeat_monomer_and_HOR_inference) — Centromeric repeat monomer identification and HOR‑structure inference

Software & Dependencies
- Python ≥ 3.9: Biopython, Pandas, Matplotlib, Seaborn, Scipy
- R ≥ 4.2: ggplot2, tidyverse, ggpubr
- Bioinformatics tools: ModDotPlot, Cent‑Mind, Jellyfish
- Shell utilities: GNU parallel for batch‑mode execution

Usage
1. Prepare input datasets including T2T genome assemblies, CENH3‑ChIP‑seq BAM files and repeat annotation results
2. Execute analysis scripts sequentially following numerical step order
3. Customize file paths and species lists for batch processing
4. Export statistical outputs and visualization figures to designated directories

Citation
If you utilize this pipeline in your research, please cite our related publications.
