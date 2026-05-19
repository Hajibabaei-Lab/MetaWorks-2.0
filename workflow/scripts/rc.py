#!/usr/bin/env python3
import argparse
import gzip

from Bio import SeqIO


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Reverse-complement all reads in a gzipped FASTQ file."
    )
    parser.add_argument("inputfile", help="Input FASTQ file (.fastq.gz)")
    parser.add_argument("outputfile", help="Output FASTQ file (.fastq.gz)")
    args = parser.parse_args(argv)

    handle_in = gzip.open(args.inputfile, "rt")
    handle_out = gzip.open(args.outputfile, "wt")

    fq = SeqIO.parse(handle_in, "fastq")
    for read in fq:
        reverse = read.reverse_complement()
        read.seq = reverse.seq
        read.letter_annotations = reverse.letter_annotations
        handle_out.write(read.format("fastq"))

    handle_in.close()
    handle_out.close()


if __name__ == "__main__":
    main()
