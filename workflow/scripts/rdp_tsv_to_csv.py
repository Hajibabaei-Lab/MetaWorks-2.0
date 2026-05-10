import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from marker_defs import get_rdp_csv_header, MARKER_TO_CONDITION


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Convert RDP TSV output to CSV with proper column headers."
    )
    parser.add_argument("filename", help="RDP output file (tab-delimited)")
    parser.add_argument("marker", help="Marker name (e.g. COI, 16S)")
    args = parser.parse_args(argv)

    df = pd.read_csv(args.filename, sep='\t', header=None)

    if args.marker not in MARKER_TO_CONDITION:
        print(f"Error: Unknown marker '{args.marker}'", file=sys.stderr)
        sys.exit(1)

    prefix = ['GlobalESV', 'Strand']
    taxonomy_cols = [c for c in get_rdp_csv_header(args.marker).split(',') if c]
    df.columns = prefix + taxonomy_cols
    print(df.to_csv(index=False, header=True))


if __name__ == "__main__":
    main()
