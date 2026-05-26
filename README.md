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

## Full Pipeline Directory Structure（点击文件夹直接跳转）
- [01.Genome preprocessing](https://github.com/wildrice6/Wild-Rice-Centromere-Project/tree/main/01.Genome%20preprocessing) — 基因组预处理、序列清洗、格式标准化
- [02.Centromere boundary prediction](https://github.com/wildrice6/Wild-Rice-Centromere-Project/tree/main/02.Centromere%20boundary%20prediction) — 着丝粒区间预测(Cent‑Mind)、CENH3富集区整合
- [03.Repeat annotation](https://github.com/wildrice6/Wild-Rice-Centromere-Project/tree/main/03.Repeat%20annotation) — 全基因组重复序列注释、卫星序列挖掘
- [04.CENH3‑ChIP‑seq analysis](https://github.com/wildrice6/Wild-Rice-Centromere-Project/tree/main/04.CENH3-ChIP-seq%20analysis) — ChIP‑seq mapping、peak calling、功能着丝粒定位
- [05.CRM‑LTR transposon analysis](https://github.com/wildrice6/Wild-Rice-Centromere-Project/tree/main/05.CRM-LTR%20transposon%20analysis) — CRM元件、LTR转座子注释、着丝粒富集分析
- [06.Monomer identification](https://github.com/wildrice6/Wild-Rice-Centromere-Project/tree/main/06.Monomer%20identification) — ModDotPlot调用、k‑mer分析、单体聚类
- [07.Repeat monomer and HOR inference](https://github.com/wildrice6/Wild-Rice-Centromere-Project/tree/main/07.Repeat%20monomer%20and%20HOR%20inference)
  - [step3_centromeric_monomer_identification](https://github.com/wildrice6/Wild-Rice-Centromere-Project/tree/main/07.Repeat%20monomer%20and%20HOR%20inference/step3_centromeric_monomer_identification) — 着丝粒单体鉴定、脚本批量生成
  - [step4_monomer_characteristics](https://github.com/wildrice6/Wild-Rice-Centromere-Project/tree/main/07.Repeat%20monomer%20and%20HOR%20inference/step4_monomer_characteristics) — 拷贝数、覆盖度、AT含量、种内/种间相似性
  - [step6_HORs_characteristics](https://github.com/wildrice6/Wild-Rice-Centromere-Project/tree/main/07.Repeat%20monomer%20and%20HOR%20inference/step6_HORs_characteristics) — HOR块大小、分布、距离、CENH3关联
- [08.Pan‑centromere gene analysis](https://github.com/wildrice6/Wild-Rice-Centromere-Project/tree/main/08.Pan-centromere%20gene%20analysis) — 着丝粒基因注释、eggNOG功能富集
- [09.Statistics & visualization](https://github.com/wildrice6/Wild-Rice-Centromere-Project/tree/main/09.Statistics%20%26%20visualization) — 批量R/Python绘图、统计检验、结果汇总

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

## Author
wildrice6
Research: Comparative genomics of centromeres in Oryza genus
Lab: Qian Qian / Xiaoming Zheng Group

## Citation
If you use this pipeline in your research, please cite our related work.
