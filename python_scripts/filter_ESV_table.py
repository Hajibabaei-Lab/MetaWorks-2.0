# Teresita M. Porter, August 18, 2022
# first arg is ESV.table.tmp filename
# second arg is for longest.orfs.fasta or orfs.fasta.nt.filtered filename2

import numpy as np
import pandas as pd
import sys
from Bio import SeqIO

# read in ESV.table.tmp
filename = sys.argv[1]
df = pd.read_csv(filename, sep='\t')

# read in longest.orfs.fasta to get good ESV ids
filename2 = sys.argv[2]
headers = []
for record in SeqIO.parse(open(filename2), 'fasta'):
	headers.append(record.id)

# filter ESV.table by GlobalESVs in taxonomy file
df_filtered = df[df["#OTU ID"].isin(headers)]
print(df_filtered.to_csv(sep='\t', index=False))

