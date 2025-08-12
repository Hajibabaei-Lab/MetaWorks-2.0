#!/usr/bin/env python3
import sys

nt_file = sys.argv[1]
aa_file = sys.argv[2]

outfile1 = nt_file.replace(".nt.hmm", ".nt.filtered.hmm")
outfile2 = aa_file.replace(".aa.hmm", ".aa.filtered.hmm")


# --- Helper functions to parse FASTA header lines ---

def parse_nt_header(line):
    """
    Parse an nt header line.
    Expected formats:
      >lcl|KR389058:1-588 ORF1_KR389058:0:588
      >lcl|... ORF1_BOLD:AAE3122:0:86  (for BOLD records)
    Returns a tuple: (orf, otu, length) or None if parsing fails.
    """
    parts = line.strip().split()
    if len(parts) < 2:
        return None
    orfline = parts[1]  # e.g. "ORF1_KR389058:0:588" or "ORF1_BOLD:AAE3122:0:86"
    orf_parts = orfline.split(":")
    if "BOLD" in orfline:
        # Expected parts: [ "ORF1_BOLD", "AAE3122", "0", "86" ]
        if len(orf_parts) < 4:
            return None
        orfbold = orf_parts[0]  # e.g. "ORF1_BOLD"
        accession = orf_parts[1]  # e.g. "AAE3122"
        try:
            start = int(orf_parts[2])
            stop  = int(orf_parts[3])
        except ValueError:
            return None
        length_val = stop - start + 1
        # Create a composite string and split by underscore.
        orfaccession = orfbold + ":" + accession  # e.g. "ORF1_BOLD:AAE3122"
        subparts = orfaccession.split("_")
        if len(subparts) < 2:
            return None
        orf = subparts[0]           # e.g. "ORF1"
        otu = subparts[1]           # e.g. "BOLD:AAE3122" (the OTU/accession)
    else:
        # Expected parts: [ "ORF1_KR389058", "0", "588" ]
        if len(orf_parts) < 3:
            return None
        orfaccession = orf_parts[0]  # e.g. "ORF1_KR389058"
        try:
            start = int(orf_parts[1])
            stop  = int(orf_parts[2])
        except ValueError:
            return None
        length_val = stop - start + 1
        subparts = orfaccession.split("_")
        if len(subparts) < 2:
            return None
        orf = subparts[0]           # e.g. "ORF1"
        otu = subparts[1]           # e.g. "KR389058"
    return orf, otu, length_val

def parse_aa_header(line):
    """
    Parse an aa header line.
    Expected formats:
      >lcl|ORF1_KR389058:0:587
      >lcl|ORF1_BOLD:AAE3122:0:86
    Returns a tuple: (orf, otu, length) or None if parsing fails.
    """
    # For aa, take the first field and remove the prefix ">lcl|"
    parts = line.strip().split()
    if len(parts) < 1:
        return None
    header = parts[0]
    if header.startswith(">lcl|"):
        header = header[5:]
    orf_parts = header.split(":")
    if "BOLD" in header:
        if len(orf_parts) < 4:
            return None
        orfbold = orf_parts[0]   # e.g. "ORF1_BOLD"
        accession = orf_parts[1]  # e.g. "AAE3122"
        try:
            start = int(orf_parts[2])
            stop  = int(orf_parts[3])
        except ValueError:
            return None
        length_val = stop - start + 1
        orfaccession = orfbold + ":" + accession
        subparts = orfaccession.split("_")
        if len(subparts) < 2:
            return None
        orf = subparts[0]
        otu = subparts[1]
    else:
        if len(orf_parts) < 3:
            return None
        orfaccession = orf_parts[0]  # e.g. "ORF1_KR389058"
        try:
            start = int(orf_parts[1])
            stop  = int(orf_parts[2])
        except ValueError:
            return None
        length_val = stop - start + 1
        subparts = orfaccession.split("_")
        if len(subparts) < 2:
            return None
        orf = subparts[0]
        otu = subparts[1]
    return orf, otu, length_val

# --- Parse the nt FASTA file ---
with open(nt_file, "r") as f:
    nt_lines = [line.rstrip("\n") for line in f]

nt_length = {}  # {otu: {orf: length}}
nt_seq = {}     # {otu: {orf: sequence}}

i = 0
while i < len(nt_lines):
    line = nt_lines[i]
    if line.startswith(">"):
        header_data = parse_nt_header(line)
        if header_data is None:
            i += 1
            continue
        orf, otu, length_val = header_data
        # Accumulate sequence lines until next header
        seq = ""
        i += 1
        while i < len(nt_lines) and not nt_lines[i].startswith(">"):
            seq += nt_lines[i].strip()
            i += 1
        if otu not in nt_length:
            nt_length[otu] = {}
            nt_seq[otu] = {}
        nt_length[otu][orf] = length_val
        nt_seq[otu][orf] = seq
    else:
        i += 1

# --- Parse the aa FASTA file ---
with open(aa_file, "r") as f:
    aa_lines = [line.rstrip("\n") for line in f]

aa_length = {}  # {otu: {orf: length}}
aa_seq = {}     # {otu: {orf: sequence}}

i = 0
while i < len(aa_lines):
    line = aa_lines[i]
    if line.startswith(">"):
        header_data = parse_aa_header(line)
        if header_data is None:
            i += 1
            continue
        orf, otu, length_val = header_data
        seq = ""
        i += 1
        while i < len(aa_lines) and not aa_lines[i].startswith(">"):
            seq += aa_lines[i].strip()
            i += 1
        if otu not in aa_length:
            aa_length[otu] = {}
            aa_seq[otu] = {}
        aa_length[otu][orf] = length_val
        aa_seq[otu][orf] = seq
    else:
        i += 1

# --- Keep records that match between nt and aa ---
# For each otu and each orf in nt, if a matching otu/orf exists in aa, record it.
match = {}  # {otu: {orf: length}}
lengths = []  # list of nt lengths for matched records (not used further)

for otu in nt_length:
    for orf in nt_length[otu]:
        if otu in aa_length and orf in aa_length[otu]:
            if otu not in match:
                match[otu] = {}
            match[otu][orf] = nt_length[otu][orf]
            lengths.append(nt_length[otu][orf])

# --- For each otu in the match set, choose the longest ORF and output its nt and aa sequences ---
with open(outfile1, "w") as out_nt, open(outfile2, "w") as out_aa:
    for otu in match:
        # Sort the ORFs for this OTU by ascending length and choose the longest one.
        orfs_sorted = sorted(match[otu], key=lambda x: match[otu][x])
        if not orfs_sorted:
            continue
        longest_orf = orfs_sorted[-1]
        nt_sequence = nt_seq[otu][longest_orf]
        aa_sequence = aa_seq[otu][longest_orf]
        out_nt.write(f">{otu}\n{nt_sequence}\n")
        out_aa.write(f">{otu}\n{aa_sequence}\n")
