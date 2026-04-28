import sys

import pandas as pd
from Bio import SeqIO
from marker_defs import get_rdp_csv_header, MARKER_TO_CONDITION

filename = sys.argv[1]
headers = []
for record in SeqIO.parse(open(filename), 'fasta'):
    headers.append(record.id)

filename2 = sys.argv[2]
df = pd.read_csv(filename2, sep='\t', header=None)

marker = sys.argv[3]

if marker not in MARKER_TO_CONDITION:
    print(f"Warning: Unknown marker '{marker}', outputting without column headers", file=sys.stderr)
    print(df.to_csv(index=False, header=False))
    sys.exit(0)

prefix = ['GlobalESV', 'Strand']
taxonomy_cols = [c for c in get_rdp_csv_header(marker).split(',') if c]
df_filtered = df[df[0].isin(headers)]
df_filtered.columns = prefix + taxonomy_cols
print(df_filtered.to_csv(index=False, header=True))
