# Teresita M. Porter, August 20, 2022
# first arg is for rdp.out.tmp filename
# second arg is config["marker"] to pick right set of headers

import numpy as np
import pandas as pd
import sys
from Bio import SeqIO

# read in rdp.out.tmp
filename = sys.argv[1]
df = pd.read_csv(filename, sep='\t', header=None)

# read in marker
marker = sys.argv[2]

tax3_abund6 = ['COI', 'rbcL_eukaryota', 'rbcL_landPlant', '12S_fish', '12S_vertebrate']
tax3_abund7 = ['16S', '28S_fungi']
tax3_abund9 = ['18S_eukaryota']

# add correct headers
if marker in tax3_abund6:

	df.columns = ['GlobalESV','Strand','Root','RootRank','rBP','SuperKingdom','SuperKingdomRank','skBP','Kingdom','KingdomRank','kBP','Phylum','PhylumRank','pBP','Class','ClassRank','cBP','Order','OrderRank','oBP','Family','FamilyRank','fBP','Genus','GenusRank','gBP','Species','SpeciesRank','sBP']
	print(df.to_csv(index=False, header=True))

elif marker in tax3_abund7:

	df.columns = ['GlobalESV','Strand','Domain','DomainRank','dBP','Phylum','PhylumRank','pBP','Class','ClassRank','cBP','Order','OrderRank','oBP','Family','FamilyRank','fBP','Genus','GenusRank','gBP']
	print(df.to_csv(index=False, header=True))

elif marker == tax3_abund9:

	df.columns = ['GlobalESV','Strand','Root','RootRank','rBP','Domain','DomainRank','dBP','Kingdom','KingdomRank','kBP','Phylum','PhylumRank','pBP','Class','ClassRank','cBP','Order','OrderRank','oBP','Family','FamilyRank','fBP','Genus','GenusRank','gBP']
	print(df.to_csv(index=False, header=True))



