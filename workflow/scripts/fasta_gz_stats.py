#!/usr/bin/env python3
import argparse
import gzip
import statistics
import sys


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Compute sequence length statistics from a gzipped FASTA file."
    )
    parser.add_argument("input_file", help="Input FASTA file (.fasta.gz)")
    args = parser.parse_args(argv)

    lengths = []

    try:
        with gzip.open(args.input_file, "rt") as f:
            in_record = False
            current_length = 0

            for line in f:
                if line.startswith(">"):
                    if in_record:
                        lengths.append(current_length)
                    in_record = True
                    current_length = 0
                    continue

                if not in_record:
                    continue

                stripped = line.strip()
                if stripped:
                    current_length += len(stripped)

            if in_record:
                lengths.append(current_length)
    except Exception as e:
        sys.exit(f"Error reading {args.input_file}: {e}")

    count = len(lengths)
    if not lengths:
        minimum = 0
        maximum = 0
        mean_val = 0
        median_val = 0
        mode_val = "NA"
    else:
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
