#!/usr/bin/env python3
# Alex Song, July 2025
# Script to grab ITSx output or chimera-free output and add to taxonomy file
# Account for a non-strictly formatted FASTA file with linebreaks in seq
# Usage: python3 python_scripts/add_seqs_to_tax3.py ITSxOutput/ChimeraFreeOut rdp.out.tmp > rdp.csv.tmp
import sys

if len(sys.argv) < 3:
    sys.exit("Usage: python3 python_scripts/add_seqs_to_tax3.py ITSxOutput/ChimeraFreeOut rdp.out.tmp > rdp.csv.tmp")

infile1 = sys.argv[1]
infile2 = sys.argv[2]

# Parse the first file: a FASTA file that may have sequences split across multiple lines.
esvs = {}         # dictionary to hold sequences keyed by id
current_id = None
current_seq_lines = []

with open(infile1, "r") as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith(">"):
            # If we already have an id, save its accumulated sequence.
            if current_id is not None:
                esvs[current_id] = "".join(current_seq_lines)
            # Remove the leading '>' and set the new id.
            current_id = line[1:]
            current_seq_lines = []
        else:
            # Accumulate sequence lines.
            current_seq_lines.append(line)
# Save the last record.
if current_id is not None:
    esvs[current_id] = "".join(current_seq_lines)

# Parse the second file (rdp.out.tmp) and print combined output.
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
        # If you want to print even when the sequence is missing, uncomment the following:
        # else:
        #     print(f"{record_id}\t\t{record}")
