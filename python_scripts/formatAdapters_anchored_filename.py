# Teresita M. Porter, May 5, 2021

#import numpy as np
import pandas as pd
import sys
from Bio import SeqIO
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord

# read in unique sample names into a list of strings
lines=[]
with open(sys.argv[2]) as f:
	lines = f.readlines()
#	print(lines)

# for each unique sample, print a new adapters.fasta file
for sample in lines:
	sample=sample.strip()
	filename=sample + '_' + sys.argv[1]
#	print(filename)

	# read in adapters.fasta
	record_list=[]
	for record in SeqIO.parse(sys.argv[1], "fasta"):
		record.id=sample + '_' + record.id
		record_list.append(record)
	SeqIO.write(record_list, filename, 'fasta')
