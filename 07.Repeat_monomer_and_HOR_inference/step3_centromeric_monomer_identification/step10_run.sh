cd /share/org/YZWL/yzwl_shahd/wangdan/Oryza_Centromere/moddotplot/all_batch

# 1. Extract centromere region fasta files
cat Oryza13_centromere_pos.xls | grep -v Species | perl -ne 'chomp;@tmp=split/\t/,$_;print "samtools faidx genome/$tmp[0]\_$tmp[1].fasta $tmp[2]:$tmp[3]-$tmp[4] >centromere_fasta/$tmp[0]\_$tmp[1].$tmp[2].centromere.fa\n"' | sh

# 2. moddotplot analysis: Iterate through various window sizes. 
# Note: Visualization for w50 often fails, so batch analysis is performed without plotting to generate BED files only.
ls centromere_fasta/*.fa | perl -ne 'chomp;@i=(50,100,150,200,250,300,400,500,600,700,800,900,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000);$o=(split/\//,$_)[-1];$o=~s/.fa//;foreach $i (@i){print "moddotplot static -f $_ -w $i -t 2 -o moddotplot_analysis/$o.w$i --no-plot --no-hist\n"}' >moddotplot_analysis.sh

# Generate tasks for csub submission
python generate_csub_scripts.py moddotplot_analysis.sh "moddotplot;2;/share/org/YZWL/yzwl_shahd/wangdan/Oryza_Centromere/moddotplot/all_batch" 140 moddotplot_analysis_csub

# Submit tasks
ls moddotplot_analysis_csub/*.sh | sed 's/^/csub </' | sh

# which may be incompatible with large window size analysis.

# 3. Statistics of repeats
ls moddotplot_analysis/*/*.bed >moddotplot.bed.list
cat moddotplot.bed.list | perl -ne 'chomp;@tmp=split/\./,(split/\//,$_)[1];$pre="$tmp[0]\.$tmp[1]";$w=$tmp[-1];print "python analyze_moddotplot.py $_ 20 repeat_stat/$pre.$w\n";' >repeat_stat.sh
python generate_csub_scripts.py repeat_stat.sh "repeat_stat;2;/share/org/YZWL/yzwl_shahd/wangdan/Oryza_Centromere/moddotplot/all_batch" 60 repeat_stat_csub

ls repeat_stat/*repeats.txt | awk -F "_top" '{print $1}' >repeats.list
ls repeat_stat_csub_remain1468/*.sh | sed 's/^/csub </' | sh

# 4. Merge repeat results and perform dynamic window analysis to select the optimal window size
# 4.1. Merge repeat results
python aggregate_repeats_results.py repeat_stat repeat_summary/merge_all_repeats.txt

# 4.2. Identify the optimal window size
cat Oryza13_centromere_pos.xls | grep -v "Hap" | cut -f 1,2 | sort | uniq | sed 's/\t/_/' | perl -ne 'chomp;print "python analyze_centromere_dynamics.py repeat_summary_new/merge_all_repeats.txt $_ repeat_summary_new/$_\n"' >repeat_summary_new.sh

# 5. Annotate Top-1 repeats in the optimal window and analyze correspondence with other genomic elements.

# 6. Extract Top-1 sequences: Extract fasta sequences corresponding to the identified optimal window size.
# Processing all materials
cat repeat_summary_new/*/*summary_metrics_classified.tsv | cut -f 1,2,4 | grep -v "Material" | perl -ne 'chomp;($s,$chr,$w)=split/\t/,$_;print "/share/org/YZWL/yzwl_shahd/wangdan/Oryza_Centromere/moddotplot/all_batch/repeat_stat/$s.$chr.w$w\_top20_regions_list.txt\n"' >06.repeat_top1_fa.list
python extract_top1_seqs.py --pos_list 06.repeat_top1_fa.list --ref_dir /share/org/YZWL/yzwl_shahd/alltest/0823PAN_CEN/00.rawdata/genome/ --out_dir 06.repeat_top1_fa_new

# 7. Generate supplementary tables
ls repeat_summary_new/*/*contribution_details.tsv | perl -ne 'chomp;$pre=(split/_contribution/,(split/\//,$_)[-1])[0];open(IN,"<$_");while(<IN>){print "$pre\t$_";}close IN;' >repeat_summary_new.all.contribution_details.tsv
ls repeat_stat/*top20_repeats.txt | perl -ne 'chomp;@pre=split/\./,(split/\//,$_)[-1];$w=(split/\_/,$pre[2])[0];open(IN,"<$_");while(<IN>){print "$pre[0]\t$pre[1]\t$w\t$_";}close IN;' >repeat_stat.all.top20_repeats.txt

# 8. Plot the dominance ratio of Top-1 repeat sequences
python plot_top1_dominance.py repeat_summary_new.all.summary_metrics_classified.tsv plot_species_col.tsv repeat_summary_plot/Top1_repeat_ratio.h6w8.pdf --height 6 --width 8