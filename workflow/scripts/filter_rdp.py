#!/usr/bin/env python3
import argparse

import numpy as np


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Filter RDP output by removing sequences without matching ORFs."
    )
    parser.add_argument("hmm_file", help="HMMER output file")
    parser.add_argument("orfs_file", help="Filtered ORFs FASTA file")
    parser.add_argument("rdp_file", help="RDP output file to filter")
    args = parser.parse_args(argv)

    hmm_dict = {}
    scores = []

    with open(args.hmm_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            id_val = parts[2]
            try:
                score_val = float(parts[5])
            except ValueError:
                continue
            hmm_dict[id_val] = score_val
            scores.append(score_val)

    p25 = np.percentile(scores, 25)
    p75 = np.percentile(scores, 75)
    iqr = p75 - p25
    lc = p25 - (1.5 * iqr)
    uc = p75 + (1.5 * iqr)

    good_ids = [id_val for id_val, score in hmm_dict.items() if not (score < lc or score > uc)]

    orfs_dict = {}

    with open(args.orfs_file, "r") as f:
        orfs_lines = f.readlines()

    i = 0
    while i < len(orfs_lines):
        line = orfs_lines[i].rstrip("\n")
        if line.startswith(">"):
            id_line = line[1:].strip()
            if i + 1 < len(orfs_lines):
                seq_line = orfs_lines[i+1].rstrip("\n")
                orfs_dict[id_line] = seq_line
            i += 2
        else:
            i += 1

    with open(args.rdp_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if not parts:
                continue
            id_val = parts[0]
            record = "\t".join(parts[1:])
            if id_val in orfs_dict and id_val in good_ids:
                seq = orfs_dict[id_val]
                print(f"{id_val}\t{seq}\t{record}")


if __name__ == "__main__":
    main()
