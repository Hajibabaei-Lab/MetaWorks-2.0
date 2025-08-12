#!/usr/bin/env python3
# Alex Song, July 2025
# This script filters the rdp.out.tmp file to include only those records that have sequences in the LongestOrfs file.
# It outputs a tab-delimited file with the record ID,
import sys

# Usage: python3 filter_rdp.py LongestOrfs rdp.out.tmp > rdp.csv.tmp
if len(sys.argv) < 3:
    sys.exit("Usage: python3 filter_rdp.py LongestOrfs rdp.out.tmp > rdp.csv.tmp")

infile1 = sys.argv[1]  # LongestOrfs file (FASTA format)
infile2 = sys.argv[2]  # rdp.out.tmp file (tab-delimited taxonomy info)

# Parse the LongestOrfs FASTA file into a dictionary.
esvs = {}
current_id = None
current_seq_lines = []

with open(infile1, "r") as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if current_id is not None:
                # Save the accumulated sequence for the previous record.
                esvs[current_id] = "".join(current_seq_lines)
            # Remove the '>' and set current_id.
            current_id = line[1:]
            current_seq_lines = []
        else:
            current_seq_lines.append(line)
# Don't forget to add the last record.
if current_id is not None:
    esvs[current_id] = "".join(current_seq_lines)

# Parse the rdp.out.tmp file and, for each record, print the id, sequence (if available), and remaining fields.
with open(infile2, "r") as f:
    for line in f:
        line = line.rstrip("\n")
        parts = line.split("\t")
        if not parts:
            continue
        record_id = parts[0]
        record = "\t".join(parts[1:])
        if record_id in esvs:
            seq = esvs[record_id]
            print(f"{record_id}\t{seq}\t{record}")
        # If record_id is not found in esvs, do nothing (don't print an assignment).
