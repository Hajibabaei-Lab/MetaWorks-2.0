#!/usr/bin/env python3
import sys
import os
import gzip


filename = sys.argv[1]

# Open the file, using gzip if the filename ends with 'gz'
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

# First, split on dots and take the first element; then split on '/' and take the last element.
base = filename.split('.')[0]
base2 = base.split('/')[-1]

# Process each line: modify header lines, leave sequence lines unchanged.
for line in lines:
    line = line.rstrip("\n")
    if line.startswith(">"):
        # Remove the '>' and prepend base2 and an underscore, then add '>' back.
        new_header = f">{base2}_{line[1:]}"
        print(new_header)
    else:
        print(line)
