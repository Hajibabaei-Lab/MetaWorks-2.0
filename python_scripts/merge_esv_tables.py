import sys
import pandas as pd

tables = [pd.read_csv(f, sep='\t', index_col=0) for f in sys.argv[1:]]
merged = pd.concat(tables, axis=1).fillna(0).astype(int)
merged.to_csv(sys.stdout, sep='\t')