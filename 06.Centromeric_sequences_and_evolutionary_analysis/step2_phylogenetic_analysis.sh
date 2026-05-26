cd /share/org/YZWL/yzbsl_chengch/Centromere_analysis_250915/Satellite/Run4/CEN155/Tree

source /share/apps/anaconda3/bin/activate

conda activate cd-hit
cd-hit -i merge.fa -o merge.nonredundant.fasta -c 0.96

conda activate Phylogenetic

sed '/^>/s/:/_/g' merge.nonredundant.fasta > merge.nonredundant.M.fasta

sed '/^>/s/(-)//g; s/(+)//g' merge.nonredundant.M.fasta > merge.nonredundant.M2.fasta

mafft --auto --thread 4 --anysymbol --maxiterate 1000 merge.nonredundant.M2.fasta > merge.nonredundant.M2.mafft

trimal -in merge.nonredundant.M2.mafft -out merge.nonredundant.M2.trimal -automated1

fasttree -nt  -gtr merge.nonredundant.M2.trimal > merge.nonredundant.M2_CEN155.tree