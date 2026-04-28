import sys

import pandas as pd
from marker_defs import get_rdp_csv_header, MARKER_TO_CONDITION

filename = sys.argv[1]
df = pd.read_csv(filename, sep='\t', header=None)

marker = sys.argv[2]

if marker not in MARKER_TO_CONDITION:
    print(f"Error: Unknown marker '{marker}'", file=sys.stderr)
    sys.exit(1)

prefix = ['GlobalESV', 'Strand']
taxonomy_cols = [c for c in get_rdp_csv_header(marker).split(',') if c]
df.columns = prefix + taxonomy_cols
print(df.to_csv(index=False, header=True))
