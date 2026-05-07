#!/usr/bin/env python3
import argparse


def parse_nt_header(line):
    parts = line.strip().split()
    if len(parts) < 2:
        return None
    orfline = parts[1]
    orf_parts = orfline.split(":")
    if "BOLD" in orfline:
        if len(orf_parts) < 4:
            return None
        orfbold = orf_parts[0]
        accession = orf_parts[1]
        try:
            start = int(orf_parts[2])
            stop = int(orf_parts[3])
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
        orfaccession = orf_parts[0]
        try:
            start = int(orf_parts[1])
            stop = int(orf_parts[2])
        except ValueError:
            return None
        length_val = stop - start + 1
        subparts = orfaccession.split("_")
        if len(subparts) < 2:
            return None
        orf = subparts[0]
        otu = subparts[1]
    return orf, otu, length_val


def parse_aa_header(line):
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
        orfbold = orf_parts[0]
        accession = orf_parts[1]
        try:
            start = int(orf_parts[2])
            stop = int(orf_parts[3])
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
        orfaccession = orf_parts[0]
        try:
            start = int(orf_parts[1])
            stop = int(orf_parts[2])
        except ValueError:
            return None
        length_val = stop - start + 1
        subparts = orfaccession.split("_")
        if len(subparts) < 2:
            return None
        orf = subparts[0]
        otu = subparts[1]
    return orf, otu, length_val


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Consolidate nt/aa ORF files, matching records and keeping longest ORF per OTU."
    )
    parser.add_argument("nt_file", help="Nucleotide ORF FASTA file")
    parser.add_argument("aa_file", help="Amino acid ORF FASTA file")
    parser.add_argument("output_nt", help="Output nucleotide FASTA file")
    parser.add_argument("output_aa", help="Output amino acid FASTA file")
    args = parser.parse_args(argv)

    with open(args.nt_file, "r") as f:
        nt_lines = [line.rstrip("\n") for line in f]

    nt_length = {}
    nt_seq = {}

    i = 0
    while i < len(nt_lines):
        line = nt_lines[i]
        if line.startswith(">"):
            header_data = parse_nt_header(line)
            if header_data is None:
                i += 1
                continue
            orf, otu, length_val = header_data
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

    with open(args.aa_file, "r") as f:
        aa_lines = [line.rstrip("\n") for line in f]

    aa_length = {}
    aa_seq = {}

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

    match = {}

    for otu in nt_length:
        for orf in nt_length[otu]:
            if otu in aa_length and orf in aa_length[otu]:
                if otu not in match:
                    match[otu] = {}
                match[otu][orf] = nt_length[otu][orf]

    with open(args.output_nt, "w") as out_nt, open(args.output_aa, "w") as out_aa:
        for otu in match:
            orfs_sorted = sorted(match[otu], key=lambda x: match[otu][x])
            if not orfs_sorted:
                continue
            longest_orf = orfs_sorted[-1]
            nt_sequence = nt_seq[otu][longest_orf]
            aa_sequence = aa_seq[otu][longest_orf]
            out_nt.write(f">{otu}\n{nt_sequence}\n")
            out_aa.write(f">{otu}\n{aa_sequence}\n")


if __name__ == "__main__":
    main()
