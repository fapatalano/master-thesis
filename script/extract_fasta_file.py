import os
import sys
from Bio import SeqIO
import pandas as pd

path = './fasta_files'
codes_tsv = open('code_dimers.tsv')
tsv_read = pd.read_csv(codes_tsv, sep='\t')
df_codes = tsv_read['code']
codes = [code[:4] for code in df_codes]
seqiter = SeqIO.parse(open('allseqs.fa'), 'fasta')
SeqIO.write((seq for seq in seqiter if seq.id[:4] in codes), 'output.fasta', "fasta")
# for seq in seqiter:
#     if seq.id[:4] not in codes:
#         continue
#     path = os.path.join('./fasta_files', seq.id + ".fasta")
#     sys.stdout = open(path, 'w')
#     SeqIO.write(seq, sys.stdout, "fasta")

