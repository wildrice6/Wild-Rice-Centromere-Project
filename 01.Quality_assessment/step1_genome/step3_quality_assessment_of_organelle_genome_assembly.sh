conda activate himt
for file in *.fasta; do himt assess -i ${file} -o ${file/.fasta/}; done