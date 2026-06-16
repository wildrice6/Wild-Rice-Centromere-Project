#!/bin/bash
#CSUB -J bowtie2_samtools
#CSUB -q c01
#CSUB -o %J.out
#CSUB -e %J.error
#CSUB -n 88
#CSUB -I
#CSUB -R span[hosts=1]
#CSUB -cwd /share/home/yzbsl_chengch/Run_out

cd /share/org/YZWL/yzbsl_chengch/ChIP_seq_CENH3_2/Analysis/YZ002/YZ002_merged

source /share/apps/anaconda3/bin/activate
conda activate myenv

bowtie2-build -f YZ002_merge.fasta --threads 35 YZ002_merge

bowtie2 -p 35 -x YZ002_merge -1 /share/org/YZWL/yzbsl_chengch/ChIP_seq_CENH3_2/Analysis/YZ002/YZ002/YZ002_1.clean.paired.R1.fq -2 /share/org/YZWL/yzbsl_chengch/ChIP_seq_CENH3_2/Analysis/YZ002/YZ002/YZ002_1.clean.paired.R2.fq --very-sensitive --no-mixed --no-discordant --maxins 800 --rg-id YZ002_merge  --rg "PL:ILLUMINA" --rg "SM:YZ002_merge" -S YZ002_merge_1_ChIP_input.sam

bowtie2 -p 35 -x YZ002_merge -1 /share/org/YZWL/yzbsl_chengch/ChIP_seq_CENH3_2/Analysis/YZ002/YZ002/YZ002_1_cenH3.clean.paired.R1.fq -2 /share/org/YZWL/yzbsl_chengch/ChIP_seq_CENH3_2/Analysis/YZ002/YZ002/YZ002_1_cenH3.clean.paired.R2.fq --very-sensitive --no-mixed --no-discordant --maxins 800 --rg-id YZ002_merge  --rg "PL:ILLUMINA" --rg "SM:YZ002_merge" -S YZ002_merge_1_ChIP_IP.sam

conda activate samtools

samtools view -@ 35 -bS YZ002_merge_1_ChIP_input.sam > YZ002_merge_1_ChIP_input.bam

samtools view -@ 35 -bS YZ002_merge_1_ChIP_IP.sam > YZ002_merge_1_ChIP_IP.bam

samtools sort -@ 35 YZ002_merge_1_ChIP_input.bam -o YZ002_merge_1_ChIP_input.sorted.bam

samtools sort -@ 35 YZ002_merge_1_ChIP_IP.bam -o YZ002_merge_1_ChIP_IP.sorted.bam

samtools index -@ 35  YZ002_merge_1_ChIP_input.sorted.bam

samtools index -@ 35  YZ002_merge_1_ChIP_IP.sorted.bam

samtools flagstat -@ 35  YZ002_merge_1_ChIP_input.sorted.bam > YZ002_merge_1_ChIP_input.sorted.flagstat

samtools flagstat -@ 35  YZ002_merge_1_ChIP_IP.sorted.bam > YZ002_merge_1_ChIP_IP.sorted.flagstat


samtools view YZ002_merge_1_ChIP_input.sam -F 256 -F 2048 -F4 > YZ002_merge_1_ChIP_input.best_alignment.sam

head -n 100 YZ002_merge_1_ChIP_input.sam | grep "^@" > head_YZ002_merge_1_ChIP_input

cat head_YZ002_merge_1_ChIP_input YZ002_merge_1_ChIP_input.best_alignment.sam >> YZ002_merge_1_ChIP_input.best_alignment_head.sam

samtools view -bS YZ002_merge_1_ChIP_input.best_alignment_head.sam | samtools sort -o YZ002_merge_1_ChIP_input.filtered_sorted.bam

samtools index -@ 35 YZ002_merge_1_ChIP_input.filtered_sorted.bam


samtools view YZ002_merge_1_ChIP_IP.sam -F 256 -F 2048 -F4 > YZ002_merge_1_ChIP_IP.best_alignment.sam

head -n 100 YZ002_merge_1_ChIP_IP.sam | grep "^@" > head_YZ002_merge_1_ChIP_IP

cat head_YZ002_merge_1_ChIP_IP YZ002_merge_1_ChIP_IP.best_alignment.sam >> YZ002_merge_1_ChIP_IP.best_alignment_head.sam

samtools view -bS YZ002_merge_1_ChIP_IP.best_alignment_head.sam | samtools sort -o YZ002_merge_1_ChIP_IP.filtered_sorted.bam

samtools index -@ 35 YZ002_merge_1_ChIP_IP.filtered_sorted.bam


conda activate deeptools

bamCompare -b1 YZ002_merge_1_ChIP_IP.filtered_sorted.bam -b2 YZ002_merge_1_ChIP_input.filtered_sorted.bam --outFileFormat bedgraph -o YZ002_merge_1_log2ratio_2k.bdg --binSize 2000 --operation log2 -p 5 --extendReads

bamCompare -b1 YZ002_merge_1_ChIP_IP.filtered_sorted.bam -b2 YZ002_merge_1_ChIP_input.filtered_sorted.bam --outFileFormat bedgraph -o YZ002_merge_1_log2ratio_10b.bdg --binSize 10 --operation log2 -p 5 --extendReads
