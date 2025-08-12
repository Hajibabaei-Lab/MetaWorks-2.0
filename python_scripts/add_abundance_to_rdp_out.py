#!/usr/bin/env python3
# Yaye (Alex) Song, July 2025
# Script to add read abundance from ESV.table to rdp.out
# Prints to STDOUT so redirect to a file in snakemake
# Usage: python3 python_scripts/add_abundance_to_rdp_out.py ESV.table rdp.csv.tmp [header] > rdp.csv.tmp2

import sys

if len(sys.argv) < 3:
    sys.exit("Usage: python3 python_scripts/add_abundance_to_rdp_out.py ESV.table rdp.csv.tmp [header] > rdp.csv.tmp2")

table_file = sys.argv[1]
rdp_file = sys.argv[2]
# Optional header parameter (if provided, print as the first line of output)
header = sys.argv[3] if len(sys.argv) >= 4 else None

# --------------------------
# Parse the RDP taxonomic assignment file.
# Build a dictionary:
#   key: global OTU (first part of the header, split by semicolon)
#   value: taxonomic assignment fields joined by commas.
# --------------------------
assignment = {}

with open(rdp_file, "r") as f:
    for line in f:
        line = line.rstrip("\n")
        if not line:
            continue
        parts = line.split("\t")
        if not parts:
            continue
        idline = parts.pop(0)
        id_parts = idline.split(";")
        global_otu = id_parts[0]
        assignment_str = ",".join(parts)
        assignment[global_otu] = assignment_str

# --------------------------
# Parse the ESV table file.
# The first line is a header row (split on whitespace, dropping the first two columns).
# Each subsequent line has a global OTU followed by abundance values per sample.
# Build a nested dictionary (table_data):
#    key: global OTU,
#    value: dict mapping sample name to abundance.
# --------------------------
table_data = {}
headers = []

with open(table_file, "r") as f:
    lines = f.readlines()

for i, line in enumerate(lines):
    line = line.rstrip("\n")
    if i == 0:
        # Header row: split on whitespace and drop the first two columns.
        headers = line.split()
        if len(headers) >= 2:
            headers = headers[2:]
    else:
        parts = line.split()
        if not parts:
            continue
        global_otu = parts.pop(0).strip()
        if len(parts) != len(headers):
            print("Check FAILED", file=sys.stderr)
        sample_abundances = {}
        for j, abund_str in enumerate(parts):
            sample = headers[j]
            try:
                abund = int(abund_str)
            except ValueError:
                try:
                    abund = float(abund_str)
                except ValueError:
                    abund = 0
            sample_abundances[sample] = abund
        table_data[global_otu] = sample_abundances

# --------------------------
# Print header (if provided) and then loop through the assignments.
# Only output records with abundance >= 3.
# --------------------------
if header:
    print(header)

for global_otu, tax_assignment in assignment.items():
    if global_otu in table_data:
        sample_abundances = table_data[global_otu]
        for sample, abund in sample_abundances.items():
            if abund >= 3:
                print(f"{global_otu},{sample},{abund},{tax_assignment}")
    else:
        print(f"Cannot find global_otu {global_otu} in table", file=sys.stderr)
