import argparse
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent))
from marker_defs import get_rdp_csv_header, MARKER_TO_CONDITION


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Filter RDP taxonomy output by FASTA headers and add column names."
    )
    parser.add_argument("filename", help="FASTA file with sequence headers to keep")
    parser.add_argument("filename2", help="RDP taxonomy output (tab-delimited)")
    parser.add_argument("marker", help="Marker name (e.g. COI, 16S)")
    args = parser.parse_args(argv)

    headers = []
    for record in SeqIO.parse(open(args.filename), 'fasta'):
        headers.append(record.id)

    max_cols = max(len(line.split('\t')) for line in open(args.filename2) if line.strip())
    df = pd.read_csv(args.filename2, sep='\t', header=None, names=range(max_cols))

    if args.marker not in MARKER_TO_CONDITION:
        print(f"Warning: Unknown marker '{args.marker}', outputting without column headers", file=sys.stderr)
        print(df.to_csv(index=False, header=False))
        sys.exit(0)

    prefix = ['GlobalESV', 'Strand']
    taxonomy_cols = [c for c in get_rdp_csv_header(args.marker).split(',') if c]
    expected = len(prefix) + len(taxonomy_cols)
    actual = len(df.columns)
    if actual > expected:
        taxonomy_cols += [f'extra_rank{i}' for i in range(actual - expected)]
    elif actual < expected:
        taxonomy_cols = taxonomy_cols[: actual - len(prefix)]
    df_filtered = df[df[0].isin(headers)]
    df_filtered.columns = prefix + taxonomy_cols
    print(df_filtered.to_csv(index=False, header=True))


if __name__ == "__main__":
    main()
