# Wild‑Rice‑Centromere‑Project
**Oryza Genus Centromere Comparative Genomics Analysis Pipeline**
稻属野生稻着丝粒比较基因组学全套分析脚本，用于26个稻属物种着丝粒区域解析、卫星重复、HOR结构、转座子、表观关联、进化动力学研究，适配T2T无间隙基因组，用于学术论文数据分析。

## Project Overview
This repository contains all bioinformatics scripts for centromere‑focused comparative genomic analysis across **26 Oryza species (AA/BB/CC/EE etc. genomes)**.
Main research topics:
- Centromere boundary identification & CENH3‑ChIP‑seq integration
- Centromeric satellite monomer (CEN126/CEN155) identification & classification
- HOR (High‑Order Repeat) structure detection, scoring and structural characterization
- CRM‑LTR retrotransposon enrichment and centromere‑transposon interplay
- Monomer similarity, AT‑content, copy number and centromeric coverage analysis
- HOR block distribution, adjacent distance dynamics, species‑level comparison
- Pan‑centromere gene density, functional enrichment and evolutionary patterns
- Batch plotting and statistical visualization for manuscript figures

## Full Pipeline Directory Structure
- [01.Quality_assessment](https://github.com/wildrice6/Wild-Rice-Centromere-Project/tree/main/01.Quality_assessment) — 基因组质量评估、序列质控
- [02.SV_between_haplotypes](https://github.com/wildrice6/Wild-Rice-Centromere-Project/tree/main/02.SV_between_haplotypes) — 单倍型间结构变异分析
- [03.Centromere_region_Identification](https://github.com/wildrice6/Wild-Rice-Centromere-Project/tree/main/03.Centromere_region_Identification) — 着丝粒区域定位与边界识别
- [04.Centromere_position_normalization_and_haplotype_similarity_analysis](https://github.com/wildrice6/Wild-Rice-Centromere-Project/tree/main/04.Centromere_position_normalization_and_haplotype_similarity_analysis) — 着丝粒位置标准化、单倍型相似性分析
- [05.Centromere_integrated_annotation_system](https://github.com/wildrice6/Wild-Rice-Centromere-Project/tree/main/05.Centromere_integrated_annotation_system) — 着丝粒整合注释体系构建
- [06.Centromeric_sequences_and_evolutionary_analysis](https://github.com/wildrice6/Wild-Rice-Centromere-Project/tree/main/06.Centromeric_sequences_and_evolutionary_analysis) — 着丝粒序列挖掘与进化动力学分析
- [07.Repeat_monomer_and_HOR_inference](https://github.com/wildrice6/Wild-Rice-Centromere-Project/tree/main/07.Repeat_monomer_and_HOR_inference)

## Software & Dependencies
- Python ≥3.9: Biopython, Pandas, Matplotlib, Seaborn, Scipy
- R ≥4.2: ggplot2, tidyverse, ggpubr
- Bioinformatics tools: ModDotPlot, Cent‑Mind, Jellyfish, MMseqs2, Syri, Mummer, eggNOG‑mapper
- Shell: GNU parallel for batch execution

## Usage
1. Organize T2T genome assemblies, CENH3‑ChIP bam files and repeat annotation results
2. Run scripts sequentially by step number
3. Modify file paths and species list for batch analysis
4. Statistical results and figures are output in separate folders

## Citation
If you use this pipeline in your research, please cite our related work.
