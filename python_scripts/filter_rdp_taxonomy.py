# Teresita M. Porter, August 18, 2022
# first arg is ITSx_out.ITS1/2.fasta.2 or longest.orfs.fasta or orfs.fasta.nt.filtered filename
# second arg is for rdp.out.tmp filename
# third arg is config["marker"] to pick right set of headers

import numpy as np
import pandas as pd
import sys
from Bio import SeqIO

# read in longest.orfs.fasta
filename = sys.argv[1]
headers = []
for record in SeqIO.parse(open(filename), 'fasta'):
	headers.append(record.id)


# read in rdp.out.tmp
filename2 = sys.argv[2]
df = pd.read_csv(filename2, sep='\t', header=None)

# read in marker
marker = sys.argv[3]

orf3_tax4_abund12 = ['COI', 'rbcL_landPlant', 'rbcL_eukaryota']
orf3_tax4_abund11 = ['rbcL_diatom']

# filter rdp.out.tmp by ESVs in longest.orfs.fasta headers
df_filtered = df[df[0].isin(headers)]

if marker in orf3_tax4_abund12:

	df_filtered.columns = ['GlobalESV','Strand','Root','RootRank','rBP','SuperKingdom','SuperKingdomRank','skBP','Kingdom','KingdomRank','kBP','Phylum','PhylumRank','pBP','Class','ClassRank','cBP','Order','OrderRank','oBP','Family','FamilyRank','fBP','Genus','GenusRank','gBP','Species','SpeciesRank','sBP']
	print(df_filtered.to_csv(index=False, header=True))

elif marker in orf3_tax4_abund11:

	df_filtered.columns = ['GlobalESV','Strand','Root','RootRank','rBP','Domain','DomainRank','dBP','Kingdom','KingdomRank','kBP','SubKingdom','SubKingdomRank','skBP','Phylum','PhylumRank','pBP','Class','ClassRank','cBP','Order','OrderRank','oBP','Family','FamilyRank','fBP','Genus','GenusRank','gBP','Species','SpeciesRank','sBP']
	print(df_filtered.to_csv(index=False, header=True))

elif marker == 'ITS_fungi':

	df_filtered.columns = ['GlobalESV','Strand','Root','RootRank','rBP','Kingdom','KingdomRank','kBP','Phylum','PhylumRank','pBP','Class','ClassRank','cBP','Order','OrderRank','oBP','Family','FamilyRank','fBP','Genus','GenusRank','gBP','Species','SpeciesRank','sBP']
	print(df_filtered.to_csv(index=False, header=True))



