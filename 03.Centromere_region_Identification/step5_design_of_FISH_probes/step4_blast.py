import argparse
import os

# Initialize argument parser
parser = argparse.ArgumentParser(description='A command-line wrapper for various BLAST sequence alignment modes.')
parser.add_argument('-m', '--model', help='Alignment mode selection. Options: blastn, short_blastn, matrix_blastn, blastp', default='blastn')
parser.add_argument('-t', '--target', help='Target sequence file (Database)')
parser.add_argument('-q', '--query', help='Query sequence file')

# Argument parsing
args = parser.parse_args()
model = args.model
input_target = args.target
input_query = args.query

class blast:
    def __init__(self, target, query):
        self.target = target
        self.query = query

    def mvdb(self):
        """Organizes generated database files into a dedicated directory."""
        target_pre = self.target.rsplit('.', 1)[0]
        query_pre = self.query.rsplit('.', 1)[0]
        command1 = f"mkdir {target_pre}.{query_pre}_db"
        command2 = f"mv {target_pre}.{query_pre}_db.* {target_pre}.{query_pre}_db"
        os.system(command1)
        os.system(command2)


class blastn(blast):
    def __init__(self, target, query):
        super().__init__(target, query)

    def dobidui(self):
        """Executes standard nucleotide BLAST (blastn)."""
        target_pre = self.target.rsplit('.', 1)[0]
        query_pre = self.query.rsplit('.', 1)[0]
        command1 = f"makeblastdb -in {self.target} -dbtype nucl -out {target_pre}.{query_pre}_db"
        command2 = f"blastn -query {self.query} -db {target_pre}.{query_pre}_db -outfmt 6 -out {target_pre}.{query_pre}.txt"
        os.system(command1)
        os.system(command2)

class short_blastn(blast):
    def __init__(self, target, query):
        super().__init__(target, query)

    def dobidui(self):
        """Executes BLAST optimized for short nucleotide sequences."""
        target_pre = self.target.rsplit('.', 1)[0]
        query_pre = self.query.rsplit('.', 1)[0]
        command1 = f"makeblastdb -in {self.target} -dbtype nucl -out {target_pre}.{query_pre}_db"
        command2 = f"blastn -query {self.query} -db {target_pre}.{query_pre}_db -outfmt 6 -out {target_pre}.{query_pre}.txt -task blastn-short -word_size 7 -evalue 1"
        os.system(command1)
        os.system(command2)

class matrix_blastn(blast):
    def __init__(self, target, query):
        super().__init__(target, query)

    def dobidui(self):
        """Executes nucleotide BLAST restricted to a single High-scoring Segment Pair (HSP)."""
        target_pre = self.target.rsplit('.', 1)[0]
        query_pre = self.query.rsplit('.', 1)[0]
        command1 = f"makeblastdb -in {self.target} -dbtype nucl -out {target_pre}.{query_pre}_db"
        command2 = f"blastn -query {self.query} -db {target_pre}.{query_pre}_db -outfmt 6 -out {target_pre}.{query_pre}.txt -max_hsps 1"
        os.system(command1)
        os.system(command2)

class blastp(blast):
    def __init__(self, target, query):
        super().__init__(target, query)

    def dobidui(self):
        """Executes protein BLAST (blastp)."""
        target_pre = self.target.rsplit('.', 1)[0]
        query_pre = self.query.rsplit('.', 1)[0]
        command1 = f"makeblastdb -in {self.target} -dbtype prot -out {target_pre}.{query_pre}_db"
        command2 = f"blastp -query {self.query} -db {target_pre}.{query_pre}_db -outfmt 6 -out {target_pre}.{query_pre}.txt"
        os.system(command1)
        os.system(command2)


# Main execution logic
if model == 'blastn':
    my_blast = blastn(input_target, input_query)
    my_blast.dobidui()
    my_blast.mvdb()
elif model == 'short_blastn':
    my_blast = short_blastn(input_target, input_query)
    my_blast.dobidui()
    my_blast.mvdb()
elif model == 'matrix_blastn':
    my_blast = matrix_blastn(input_target, input_query)
    my_blast.dobidui()
    my_blast.mvdb()
elif model == 'blastp':
    my_blast = blastp(input_target, input_query)
    my_blast.dobidui()
    my_blast.mvdb()
else:
    print('Error: The specified mode is not supported.')