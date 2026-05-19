#!/usr/bin/env python3
import argparse
import gzip
import statistics
import sys


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Compute read length statistics from a gzipped FASTQ file."
    )
    parser.add_argument("input_file", help="Input FASTQ file (.fastq.gz)")
    args = parser.parse_args(argv)

    lengths = []

    try:
        with gzip.open(args.input_file, "rt") as f:
            for i, line in enumerate(f):
                if i % 4 == 1:
                    lengths.append(len(line.strip()))
    except Exception as e:
        sys.exit(f"Error reading {args.input_file}: {e}")

    if not lengths:
        sys.exit("No sequence lines found in the file.")

    count = len(lengths)
    minimum = min(lengths)
    maximum = max(lengths)
    mean_val = statistics.mean(lengths)
    median_val = statistics.median(lengths)
    try:
        mode_val = statistics.mode(lengths)
    except statistics.StatisticsError:
        mode_val = "NA"

    print(f"{args.input_file}\t{count}\t{minimum}\t{maximum}\t{mean_val}\t{median_val}\t{mode_val}")


if __name__ == "__main__":
    main()
