#!/usr/bin/env python3
import sys

if len(sys.argv) < 3:
    sys.exit("Usage: {} taxon.zotus cat.denoised.nonchimeras".format(sys.argv[0]))

taxon_file = sys.argv[1]
fasta_file = sys.argv[2]

# -------------------------------
# Read the list of taxon IDs (zotus)
# -------------------------------
with open(taxon_file, "r") as f:
    tax_ids = [line.strip() for line in f if line.strip()]

# -------------------------------
# Parse the FASTA file into a dictionary mapping zotu -> sequence
# -------------------------------
fasta_dict = {}
current_zotu = None
current_seq = []

with open(fasta_file, "r") as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith(">"):
            # If a record is already in progress, save it
            if current_zotu is not None:
                fasta_dict[current_zotu] = "".join(current_seq)
            # Start a new record: remove the ">" and strip whitespace
            current_zotu = line[1:].strip()
            current_seq = []
        else:
            # Append the line (sequence part) to the current sequence list
            current_seq.append(line.strip())
    # After the loop, add the last record (if any)
    if current_zotu is not None:
        fasta_dict[current_zotu] = "".join(current_seq)

# -------------------------------
# For each taxon ID, if it exists in the FASTA dictionary, print the record
# -------------------------------
for zotu in tax_ids:
    if zotu in fasta_dict:
        print(f">{zotu}")
        print(fasta_dict[zotu])
