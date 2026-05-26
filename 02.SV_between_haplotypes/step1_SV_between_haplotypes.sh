#!/bin/bash
#CSUB -J AA_Oglu_SV
#CSUB -q c01
#CSUB -o AA_Oglu_SV
#CSUB -e AA_Oglu_SV
#CSUB -n 32
#CSUB -cwd /share/org/YZWL/yzwl_liubx/bj/7.3jjhsj/SV_output/AA_Oglu

source /share/apps/anaconda3/bin/activate
conda activate /share/org/YZWL/yzwl_liubx/miniconda3/envs/SV

minimap2 -ax asm5 -t 16 --eqx /share/org/YZWL/yzwl_liubx/bj/7.3jjhsj/Ref/00.genome_centromere_new/AA_Oglu/hap1/AA_Oglu_hap1.fasta \
    /share/org/YZWL/yzwl_liubx/bj/7.3jjhsj/Ref/00.genome_centromere_new/AA_Oglu/hap2/AA_Oglu_hap2.fasta \
    | samtools sort -O BAM - > /share/org/YZWL/yzwl_liubx/bj/7.3jjhsj/SV_output/AA_Oglu/AA_Oglu_Lhex.bam

syri -c  AA_Oglu_Lhex.bam -r /share/org/YZWL/yzwl_liubx/bj/7.3jjhsj/Ref/00.genome_centromere_new/AA_Oglu/hap1/AA_Oglu_hap1.fasta -q \
    /share/org/YZWL/yzwl_liubx/bj/7.3jjhsj/Ref/00.genome_centromere_new/AA_Oglu/hap2/AA_Oglu_hap2.fasta -F B --prefix AA_Oglu_Lhex.

plotsr --sr AA_Oglu_Lhex.syri.out --genomes \
    /share/org/YZWL/yzwl_liubx/bj/7.3jjhsj/SV_output/AA_Oglu/AA_Oglu.txt \
    -H 8 -W 5 -o AA_Oglu_Lhex_Centromere_syntenic_plot.pdf

awk '$11=="DEL" && $3-$2>49'  AA_Oglu_Lhex.syri.out  > AA_Oglu.syri.50.bed
awk '$11=="INS" && $8-$7>49'  AA_Oglu_Lhex.syri.out  >> AA_Oglu.syri.50.bed
awk '$11=="INV" && $3-$2>49'  AA_Oglu_Lhex.syri.out  >> AA_Oglu.syri.50.bed
awk '$11=="DUP" && $3-$2>49'  AA_Oglu_Lhex.syri.out  >> AA_Oglu.syri.50.bed
