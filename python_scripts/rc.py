# Teresita M. Porter, May 6/21

from Bio import SeqIO
import gzip
import sys

inputfile = sys.argv[1]
outputfile = sys.argv[2]
handle_in = gzip.open(inputfile, "rt")
handle_out = gzip.open(outputfile, "wt")

fq = SeqIO.parse(handle_in, "fastq")
for read in fq:
	reverse = read.reverse_complement()
	read.seq = reverse.seq
	read.letter_annotations = reverse.letter_annotations
	handle_out.write(read.format("fastq"))

handle_in.close()
handle_out.close()

