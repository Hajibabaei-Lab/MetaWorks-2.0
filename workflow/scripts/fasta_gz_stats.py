#!/usr/bin/env python3
import gzip
import statistics
import sys

if len(sys.argv) < 2:
    sys.exit("Usage: python3 fasta_gz_stats.py infile.fasta.gz > outfile")

input_file = sys.argv[1]
lengths: list[int] = []

try:
    with gzip.open(input_file, "rt") as f:
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
    sys.exit(f"Error reading {input_file}: {e}")

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

# Output: input filename, count, min, max, mean, median, and mode (tab-seperated)
print(f"{input_file}\t{count}\t{minimum}\t{maximum}\t{mean_val}\t{median_val}\t{mode_val}")
