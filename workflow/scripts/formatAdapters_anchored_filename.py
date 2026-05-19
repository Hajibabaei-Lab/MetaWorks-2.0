import argparse
import os

from Bio import SeqIO


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Create per-sample adapter FASTA files from a base adapters file and sample list."
    )
    parser.add_argument("adapters_file", help="Base adapters FASTA file")
    parser.add_argument("samples_file", help="File with sample names (one per line)")
    parser.add_argument("--output-dir", default=".", help="Output directory for per-sample files")
    args = parser.parse_args(argv)

    os.makedirs(args.output_dir, exist_ok=True)

    with open(args.samples_file) as f:
        samples = [line.strip() for line in f if line.strip()]

    for sample in samples:
        filename = os.path.join(args.output_dir, f"{sample}_adapters.fasta")
        record_list = []
        for record in SeqIO.parse(args.adapters_file, "fasta"):
            record.id = sample + "_" + record.id
            record.description = ""
            record_list.append(record)
        SeqIO.write(record_list, filename, "fasta")


if __name__ == "__main__":
    main()
