#!/usr/bin/env python3
import argparse
import gzip
import sys


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Rename FASTA headers by prepending the filename base."
    )
    parser.add_argument("filename", help="Input FASTA file (.fasta or .fasta.gz)")
    args = parser.parse_args(argv)

    filename = args.filename

    if filename.endswith("gz"):
        try:
            with gzip.open(filename, "rt") as f:
                lines = f.readlines()
        except Exception as e:
            sys.exit(f"Cannot open gzipped infile: {e}")
    else:
        try:
            with open(filename, "r") as f:
                lines = f.readlines()
        except Exception as e:
            sys.exit(f"Cannot open infile: {e}")

    base = filename.split('.')[0]
    base2 = base.split('/')[-1]

    for line in lines:
        line = line.rstrip("\n")
        if line.startswith(">"):
            new_header = f">{base2}_{line[1:]}"
            print(new_header)
        else:
            print(line)


if __name__ == "__main__":
    main()
