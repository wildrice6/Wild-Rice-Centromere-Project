# BUSCO assessment
busco -i ${ASSEMBLY} \
    -l embryophyta_odb10 \
    -o ${SAMPLE}_busco \
    -m genome \
    -c ${THREADS}

# LAI
/usr/bin/perl LTR_FINDER_parallel -seq ${SAMPLE}.fa -harvest_out -threads ${THREADS}
gt suffixerator -db ${SAMPLE}.fa -indexname ${SAMPLE} -tis -suf -lcp -des -ssp -sds -dna
gt ltrharvest -index ${SAMPLE} -maxlenltr 7000 >${SAMPLE}.harvest.scn
cat ${SAMPLE}.harvest.scn ${SAMPLE}.finder.combine.scn >${SAMPLE}_rawLTR_7k.scn
perl LTR_retriever -genome ${SAMPLE}.fa -inharvest ${SAMPLE}_rawLTR_7k.scn -threads ${THREADS}


# GCI
minimap2 -t ${THREADS} -ax map-hifi ${SAMPLE}.fa  ${HIFI_READS} > ${SAMPLE}.minimap2.hifi.sam
samtools view -@ ${THREADS} -Sb  ${SAMPLE}.minimap2.hifi.sam | samtools sort -@ ${THREADS} -m 2g -o ${SAMPLE}.minimap2.hifi.bam
samtools index ${SAMPLE}.minimap2.hifi.bam
paftools.js sam2paf ${SAMPLE}.minimap2.hifi.sam | sort -k6,6V -k8,8n > ${SAMPLE}.minimap2.hifi.paf

meryl count k=${N} output ${SAMPLE}.merylDB ${SAMPLE}.fa
meryl print greater-than distinct=0.9998 ${SAMPLE}.merylDB > ${SAMPLE}.repetitive_k15.txt
winnowmap -W ${SAMPLE}.repetitive_k15.txt -ax map-pb -t ${THREADS}  ${SAMPLE}.fa ${HIFI_READS} > ${SAMPLE}.winnowmap.hifi.sam
samtools view -@ ${THREADS} -Sb ${SAMPLE}.winnowmap.hifi.sam | samtools sort -@ ${THREADS} -m 2g -o ${SAMPLE}.winnowmap.hifi.bam
samtools index ${SAMPLE}.winnowmap.hifi.bam
paftools.js sam2paf ${SAMPLE}.winnowmap.hifi.sam | sort -k6,6V -k8,8n > ${SAMPLE}.winnowmap.hifi.paf

minimap2 -t ${THREADS} -ax map-ont ${SAMPLE}.fa ${ONT_READS} > minimap2.ont.sam
samtools view -@ ${THREADS} -Sb ${SAMPLE}.minimap2.ont.sam | samtools sort -@ ${THREADS} -m 2g -o ${SAMPLE}.minimap2.ont.bam
samtools index ${SAMPLE}.minimap2.ont.bam
paftools.js sam2paf ${SAMPLE}.minimap2.ont.sam | sort -k6,6V -k8,8n > ${SAMPLE}.minimap2.ont.paf

winnowmap -W ${SAMPLE}.repetitive_k15.txt -ax map-ont -t ${THREADS} ${SAMPLE}.fa ${ONT_READS} > ${SAMPLE}.winnowmap.ont.sam
samtools view -@ ${THREADS} -Sb ${SAMPLE}.winnowmap.ont.sam | samtools sort -@ ${THREADS} -m 2g -o ${SAMPLE}.winnowmap.ont.bam
samtools index ${SAMPLE}.winnowmap.ont.bam
paftools.js sam2paf ${SAMPLE}.winnowmap.ont.sam | sort -k6,6V -k8,8n > ${SAMPLE}.winnowmap.ont.paf

ulimit -n 10240  # 提高文件打开数限制
python ../GCI.py -r ${SAMPLE}.fa \
  --hifi ${SAMPLE}.winnowmap.hifi.bam ${SAMPLE}.minimap2.hifi.paf \
  --nano ${SAMPLE}.winnowmap.ont.bam ${SAMPLE}.minimap2.ont.paf \
  -t ${THREADS} -p -it pdf
  
# CRAQ
craq -g ${SAMPLE}.fa -sms ${HIFI_READS} -ngs ${SAMPLE}_R1.fq.gz,${SAMPLE}_R2.fq.gz -x map-hifi -t ${THREADS}
