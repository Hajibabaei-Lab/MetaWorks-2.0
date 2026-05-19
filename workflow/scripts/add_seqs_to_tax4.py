#!/usr/bin/env python3
import argparse


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Filter RDP taxonomy output to include only records with sequences in a FASTA file."
    )
    parser.add_argument("infile1", help="LongestOrfs FASTA file")
    parser.add_argument("infile2", help="RDP taxonomy file (tab-delimited)")
    args = parser.parse_args(argv)

    esvs = {}
    current_id = None
    current_seq_lines = []

    with open(args.infile1, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_id is not None:
                    esvs[current_id] = "".join(current_seq_lines)
                current_id = line[1:]
                current_seq_lines = []
            else:
                current_seq_lines.append(line)
    if current_id is not None:
        esvs[current_id] = "".join(current_seq_lines)

    with open(args.infile2, "r") as f:
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


if __name__ == "__main__":
    main()
