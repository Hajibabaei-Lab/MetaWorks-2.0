#!/usr/bin/env python3
import gzip
import statistics
import sys

input_file = sys.argv[1]

lengths = []

try:
    with gzip.open(input_file, "rt") as f:
        for i, line in enumerate(f):
            # For a FASTQ file, every 2nd line (0-indexed line 1, 5, 9, etc.) is a sequence.
            if i % 4 == 1:
                lengths.append(len(line.strip()))
except Exception as e:
    sys.exit(f"Error reading {input_file}: {e}")

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

# Output the statistics as a tab-separated line: input_file, count, min, max, mean, median, mode
print(f"{input_file}\t{count}\t{minimum}\t{maximum}\t{mean_val}\t{median_val}\t{mode_val}")
