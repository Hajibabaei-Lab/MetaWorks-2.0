#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Filter ESV table by ESV IDs present in ORF FASTA file."
    )
    parser.add_argument("esv_table_path", help="ESV table in TSV format")
    parser.add_argument("orfs_path", help="FASTA file with ESV IDs to keep")
    args = parser.parse_args(argv)

    esv_table_path = Path(args.esv_table_path)
    orfs_path = Path(args.orfs_path)

    if not esv_table_path.exists():
        print(f"ESV table file not found: {esv_table_path}", file=sys.stderr)
        sys.exit(1)

    if not orfs_path.exists():
        print(f"ORF FASTA file not found: {orfs_path}", file=sys.stderr)
        sys.exit(1)

    try:
        df = pd.read_csv(esv_table_path, sep="\t")
    except pd.errors.EmptyDataError:
        print(f"ESV table is empty: {esv_table_path}", file=sys.stderr)
        sys.exit(1)
    except pd.errors.ParserError:
        print(f"Failed to parse ESV table: {esv_table_path}", file=sys.stderr)
        sys.exit(1)

    try:
        headers = [record.id for record in SeqIO.parse(open(orfs_path), "fasta")]
    except Exception as exc:
        print(f"Failed to read ORF FASTA: {orfs_path}: {exc}", file=sys.stderr)
        sys.exit(1)

    if "#OTU ID" not in df.columns:
        print("ESV table missing required column: #OTU ID", file=sys.stderr)
        sys.exit(1)

    df_filtered = df[df["#OTU ID"].isin(headers)]
    print(df_filtered.to_csv(sep="\t", index=False))


if __name__ == "__main__":
    main()
