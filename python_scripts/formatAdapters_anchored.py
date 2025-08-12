# Teresita M. Porter, May 5, 2021

import numpy as np
import pandas as pd
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# vars
file_out = 'adapters.fasta'

# read in csv
filename = sys.argv[-1]
df = pd.read_csv(filename)

# create fasta header for linked adapters
df['SampleIDamplicon'] = df['SampleID'].str.cat(df['Amplicon'],sep="_")
df['SampleIDamplicon'] = df['SampleIDamplicon'].astype(str) + ";"

# replace I with N where needed
df['Forward'] = df['Forward'].str.replace('I', 'N')
df['Reverse'] = df['Reverse'].str.replace('I', 'N')

# reverse complement reverse primer
i = 0
rc = [None]*len(df['Reverse'])
for x in df['Reverse']:
	seq = Seq(x)
	rc[i] = seq.reverse_complement()
	rc[i] = str(rc[i])
	i = i + 1

df['ReverseRC'] = rc

# Linked adapters that are NOT anchored
df['LinkedAdapters'] = df['Forward'].str.cat(df['ReverseRC'], sep="...")

# Linked adapters anchored at the 5' and 3' ends
df['LinkedAdapters'] = df['LinkedAdapters'].map('^{}$'.format)

record_list=[]
with open(file_out, 'w') as f_out:
	for index, row in df.iterrows():
		record = SeqRecord(Seq(df['LinkedAdapters'].iloc[index]), 
		description = "",
		id=df['SampleIDamplicon'].iloc[index])
		record_list.append(record)
	SeqIO.write(record_list, f_out, 'fasta')

