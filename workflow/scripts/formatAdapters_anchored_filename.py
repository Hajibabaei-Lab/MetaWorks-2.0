import argparse

from Bio import SeqIO


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Create per-sample adapter FASTA files from a base adapters file and sample list."
    )
    parser.add_argument("adapters_file", help="Base adapters FASTA file")
    parser.add_argument("samples_file", help="File with sample names (one per line)")
    args = parser.parse_args(argv)

    lines = []
    with open(args.samples_file) as f:
        lines = f.readlines()

    for sample in lines:
        sample = sample.strip()
        filename = sample + '_' + args.adapters_file

        record_list = []
        for record in SeqIO.parse(args.adapters_file, "fasta"):
            record.id = sample + '_' + record.id
            record_list.append(record)
        SeqIO.write(record_list, filename, 'fasta')


if __name__ == "__main__":
    main()
