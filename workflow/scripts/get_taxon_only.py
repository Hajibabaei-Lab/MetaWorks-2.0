#!/usr/bin/env python3
import argparse


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Extract FASTA sequences matching a list of taxon IDs."
    )
    parser.add_argument("taxon_file", help="File with taxon IDs (one per line)")
    parser.add_argument("fasta_file", help="FASTA file to extract from")
    args = parser.parse_args(argv)

    with open(args.taxon_file, "r") as f:
        tax_ids = [line.strip() for line in f if line.strip()]

    fasta_dict = {}
    current_zotu = None
    current_seq = []

    with open(args.fasta_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_zotu is not None:
                    fasta_dict[current_zotu] = "".join(current_seq)
                current_zotu = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_zotu is not None:
            fasta_dict[current_zotu] = "".join(current_seq)

    for zotu in tax_ids:
        if zotu in fasta_dict:
            print(f">{zotu}")
            print(fasta_dict[zotu])


if __name__ == "__main__":
    main()
