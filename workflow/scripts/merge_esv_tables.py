#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path

import pandas as pd


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Merge multiple ESV tables into a single combined table."
    )
    parser.add_argument("table_paths", nargs="+", help="ESV table files in TSV format")
    args = parser.parse_args(argv)

    table_paths = [Path(p) for p in args.table_paths]

    for table_path in table_paths:
        if not table_path.exists():
            print(f"ESV table file not found: {table_path}", file=sys.stderr)
            sys.exit(1)

    tables = []
    for table_path in table_paths:
        try:
            df = pd.read_csv(table_path, sep="\t", index_col=0)
            tables.append(df)
        except pd.errors.EmptyDataError:
            print(f"ESV table is empty: {table_path}", file=sys.stderr)
            sys.exit(1)
        except pd.errors.ParserError:
            print(f"Failed to parse ESV table: {table_path}", file=sys.stderr)
            sys.exit(1)
        except Exception as exc:
            print(f"Failed to read ESV table: {table_path}: {exc}", file=sys.stderr)
            sys.exit(1)

    if not tables:
        print("No valid ESV tables to merge", file=sys.stderr)
        sys.exit(1)

    try:
        merged = pd.concat(tables, axis=1).fillna(0).astype(int)
    except Exception as exc:
        print(f"Failed to merge ESV tables: {exc}", file=sys.stderr)
        sys.exit(1)

    merged.to_csv(sys.stdout, sep="\t")


if __name__ == "__main__":
    main()
