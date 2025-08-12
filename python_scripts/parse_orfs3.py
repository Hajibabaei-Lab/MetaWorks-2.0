#!/usr/bin/env python3
import sys


input_file = sys.argv[1]
output_file = sys.argv[2]

# Read all lines from the input FASTA file (cds.fasta.tmp)
with open(input_file, 'r') as f:
    lines = [line.rstrip("\n") for line in f]

# Data structures: 
# nt_length: { otu: { orf: length, ... }, ... }
# nt_seq: { otu: { orf: sequence, ... }, ... }
nt_length = {}
nt_seq = {}

# Process the FASTA records
# Assumes headers are like:
#   >lcl|Otu1:2:310 ORF1_Otu1:1:309
# We use the second field (splitting by whitespace) which is "ORF1_Otu1:1:309"
# Then, split on ":" to get [ "ORF1_Otu1", "1", "309" ].
# From that, start = 1, stop = 309, so length = 309 - 1 + 1 = 309.
# Also, split "ORF1_Otu1" on "_" to extract ORF and OTU (here: "ORF1", "Otu1").
i = 0
while i < len(lines):
    line = lines[i]
    if line.startswith(">"):
        parts = line.split()
        if len(parts) < 2:
            i += 1
            continue
        orf_info = parts[1]  # e.g., "ORF1_Otu1:1:309"
        orf_parts = orf_info.split(":")
        if len(orf_parts) < 3:
            i += 1
            continue
        header_part = orf_parts[0]  # "ORF1_Otu1"
        try:
            start = int(orf_parts[1])
            stop = int(orf_parts[2])
        except ValueError:
            i += 1
            continue
        length_val = stop - start + 1

        # Split header_part by underscore to get ORF and OTU.
        header_split = header_part.split("_")
        if len(header_split) < 2:
            i += 1
            continue
        orf_id = header_split[0]  # e.g., "ORF1"
        otu_id = header_split[1]  # e.g., "Otu1"

        # Accumulate the sequence lines (which may span multiple lines) until the next header.
        seq = ""
        j = i + 1
        while j < len(lines) and not lines[j].startswith(">"):
            seq += lines[j].strip()
            j += 1

        # Save length and sequence in dictionaries (organized by OTU, then ORF)
        if otu_id not in nt_length:
            nt_length[otu_id] = {}
            nt_seq[otu_id] = {}
        nt_length[otu_id][orf_id] = length_val
        nt_seq[otu_id][orf_id] = seq

        i = j  # jump to next header or EOF
    else:
        i += 1

# For each OTU, pick the longest ORF.
nt_length_longest = {}
nt_seq_longest = {}
longest_lengths = []  # list to hold lengths of each longest ORF

for otu in nt_length:
    if nt_length[otu]:
        # Get the ORF with the maximum length for this OTU.
        longest_orf, longest_length = max(nt_length[otu].items(), key=lambda x: x[1])
        nt_length_longest[otu] = {longest_orf: longest_length}
        nt_seq_longest[otu] = {longest_orf: nt_seq[otu][longest_orf]}
        longest_lengths.append(longest_length)

# Calculate the 25th and 75th percentiles (lower and upper quartiles) of the longest ORF lengths.
sorted_lengths = sorted(longest_lengths)
n = len(sorted_lengths)
p25_index = int(n * 0.25)
p75_index = int(n * 0.75)
if p25_index >= n:
    p25_index = n - 1
if p75_index >= n:
    p75_index = n - 1
percentile25 = sorted_lengths[p25_index]
percentile75 = sorted_lengths[p75_index]
iqr = percentile75 - percentile25

# Calculate cutoff values: lower cutoff = 25th percentile - 1.5*IQR, upper cutoff = 75th percentile + 1.5*IQR.
lower_cutoff = percentile25 - (1.5 * iqr)
upper_cutoff = percentile75 + (1.5 * iqr)

# Append filtered longest ORFs (those within cutoff range) to the output file (cds.fasta).
with open(output_file, "a") as out_f:
    for otu in nt_length_longest:
        for orf in nt_length_longest[otu]:
            length_val = nt_length_longest[otu][orf]
            if lower_cutoff <= length_val <= upper_cutoff:
                seq = nt_seq_longest[otu][orf]
                out_f.write(f">{otu}\n{seq}\n")
