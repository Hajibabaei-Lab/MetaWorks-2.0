#!/usr/bin/env python3
import argparse


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Parse ORF FASTA and select longest ORF per OTU with IQR filtering."
    )
    parser.add_argument("input_file", help="Input FASTA file (cds.fasta.tmp)")
    parser.add_argument("output_file", help="Output FASTA file")
    args = parser.parse_args(argv)

    with open(args.input_file, 'r') as f:
        lines = [line.rstrip("\n") for line in f]

    nt_length = {}
    nt_seq = {}

    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith(">"):
            parts = line.split()
            if len(parts) < 2:
                i += 1
                continue
            orf_info = parts[1]
            orf_parts = orf_info.split(":")
            if len(orf_parts) < 3:
                i += 1
                continue
            header_part = orf_parts[0]
            try:
                start = int(orf_parts[1])
                stop = int(orf_parts[2])
            except ValueError:
                i += 1
                continue
            length_val = stop - start + 1

            header_split = header_part.split("_")
            if len(header_split) < 2:
                i += 1
                continue
            orf_id = header_split[0]
            otu_id = header_split[1]

            seq = ""
            j = i + 1
            while j < len(lines) and not lines[j].startswith(">"):
                seq += lines[j].strip()
                j += 1

            if otu_id not in nt_length:
                nt_length[otu_id] = {}
                nt_seq[otu_id] = {}
            nt_length[otu_id][orf_id] = length_val
            nt_seq[otu_id][orf_id] = seq

            i = j
        else:
            i += 1

    nt_length_longest = {}
    nt_seq_longest = {}
    longest_lengths = []

    for otu in nt_length:
        if nt_length[otu]:
            longest_orf, longest_length = max(nt_length[otu].items(), key=lambda x: x[1])
            nt_length_longest[otu] = {longest_orf: longest_length}
            nt_seq_longest[otu] = {longest_orf: nt_seq[otu][longest_orf]}
            longest_lengths.append(longest_length)

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

    lower_cutoff = percentile25 - (1.5 * iqr)
    upper_cutoff = percentile75 + (1.5 * iqr)

    with open(args.output_file, "w") as out_f:
        for otu in nt_length_longest:
            for orf in nt_length_longest[otu]:
                length_val = nt_length_longest[otu][orf]
                if lower_cutoff <= length_val <= upper_cutoff:
                    seq = nt_seq_longest[otu][orf]
                    out_f.write(f">{otu}\n{seq}\n")


if __name__ == "__main__":
    main()
