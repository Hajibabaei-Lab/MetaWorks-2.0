#!/usr/bin/env python3
import sys
import numpy as np

if len(sys.argv) < 4:
        sys.exit("Usage: {} hmm.txt orfs.fasta.nt.filtered rdp.out.tmp".format(sys.argv[0]))

        hmm_file = sys.argv[1]
        orfs_file = sys.argv[2]
        rdp_file = sys.argv[3]

        # -------------------------------
        # Parse HMMER output and hash scores
        # -------------------------------
        hmm_dict = {}   # key: id, value: score
        scores = []     # list of scores for percentile calculation

        with open(hmm_file, "r") as f:
            for line in f:
                line = line.rstrip("\n")
                if line.startswith("#"):
                    continue
                parts = line.split()  # split on any whitespace
                if len(parts) < 6:
                    continue  # skip incomplete lines
                id_val = parts[2]
                try:
                    score_val = float(parts[5])
                except ValueError:
                    continue
                hmm_dict[id_val] = score_val
                scores.append(score_val)

        # -------------------------------
        # Compute cutoff values using percentiles (25th and 75th)
        # -------------------------------
        p25 = np.percentile(scores, 25)
        p75 = np.percentile(scores, 75)
        iqr = p75 - p25
        lc = p25 - (1.5 * iqr)
        uc = p75 + (1.5 * iqr)

        # Get list of good ids (skip over ids with outlier scores)
        good_ids = [id_val for id_val, score in hmm_dict.items() if not (score < lc or score > uc)]
        # Note: In this script we calculate good_ids but do not use them further.
        # You could use them to filter the ORFs if desired.

        # -------------------------------
        # Parse filtered nt orfs FASTA and hash sequences
        # -------------------------------
        orfs_dict = {}  # key: id, value: sequence

        with open(orfs_file, "r") as f:
            orfs_lines = f.readlines()

        i = 0
        while i < len(orfs_lines):
            line = orfs_lines[i].rstrip("\n")
            if line.startswith(">"):
                id_line = line[1:].strip()  # remove '>' and strip whitespace
                # Next line is the sequence (assuming one-line sequences)
                if i + 1 < len(orfs_lines):
                    seq_line = orfs_lines[i+1].rstrip("\n")
                    orfs_dict[id_line] = seq_line
                i += 2
            else:
                i += 1

        # -------------------------------
        # Parse rdp.out.tmp and add ORF sequence if available
        # -------------------------------
        with open(rdp_file, "r") as f:
            for line in f:
                line = line.rstrip("\n")
                if not line:
                    continue
                parts = line.split("\t")
                if not parts:
                    continue
                id_val = parts[0]
                record = "\t".join(parts[1:])
                if id_val in orfs_dict:
                    seq = orfs_dict[id_val]
                    print(f"{id_val}\t{seq}\t{record}")
                # Else: do not print rows that were screened out
