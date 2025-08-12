#!/usr/bin/env python3
import gzip
import statistics
import sys

# Use sys.argv to get the input file.
if len(sys.argv) < 2:
    sys.exit("Usage: python3 fasta_gz_stats.py infile.fasta.gz > outfile")

input_file = sys.argv[1]
lengths = []

# Open the gzipped FASTA file in text mode.
try:
    with gzip.open(input_file, "rt") as f:
        # Read all lines from the file.
        lines = f.readlines()
except Exception as e:
    sys.exit(f"Error reading {input_file}: {e}")

allseqs = [line.rstrip("\n") for i, line in enumerate(lines) if (i + 1) % 2 == 0]

num = len(allseqs)

# For each sequence line, calculate its length.
for seq in allseqs:
    lengths.append(len(seq))

num2 = len(lengths)
if num != num2:
    print("Possible error", file=sys.stderr)

count = len(lengths)
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
