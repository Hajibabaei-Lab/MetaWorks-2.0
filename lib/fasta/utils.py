"""
FASTA/FASTQ utility functions.

This module provides utility functions for manipulating FASTA/FASTQ files.
"""

import gzip
from contextlib import ExitStack
from pathlib import Path
from typing import IO, Dict, Union

from Bio import SeqIO

from ..exceptions import FileProcessingError


def reverse_complement(input_file: Union[str, Path], output_file: Union[str, Path]) -> None:
    """
    Reverse complement all sequences in a FASTQ file.

    Args:
        input_file: Path to input FASTQ file (gzipped)
        output_file: Path to output FASTQ file (will be gzipped)

    Raises:
        FileProcessingError: If files cannot be read or written
        ValidationError: If file is invalid

    Example:
        >>> from lib.fasta import reverse_complement
        >>> reverse_complement('reads.fastq.gz', 'rc_reads.fastq.gz')

    Author: Teresita M. Porter, May 6/21
    """
    input_path = Path(input_file)
    output_path = Path(output_file)

    if not input_path.exists():
        raise FileProcessingError(
            f"Input file not found: {input_path}",
            filepath=str(input_path),
            suggestion="Check that input file exists and is readable",
        )

    try:
        # Create output directory if needed
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Open files (assumes gzipped input and output)
        with gzip.open(input_path, "rt") as handle_in, gzip.open(output_path, "wt") as handle_out:
            # Parse FASTQ file
            fq = SeqIO.parse(handle_in, "fastq")
            for read in fq:
                reverse = read.reverse_complement()
                read.seq = reverse.seq
                read.letter_annotations = reverse.letter_annotations
                handle_out.write(read.format("fastq"))

    except Exception as exc:
        raise FileProcessingError(
            f"Error processing {input_path}: {exc}",
            filepath=str(input_path),
            suggestion="Check file format and permissions",
        ) from exc


def rename_fasta_gzip(
    input_file: Union[str, Path],
    output_file: Union[str, Path],
    name_mapping: Dict[str, str],
) -> None:
    """
    Rename sequences in a FASTA file based on mapping dictionary.

    Args:
        input_file: Path to input FASTA file
        output_file: Path to output FASTA file
        name_mapping: Dictionary mapping old names to new names

    Raises:
        FileProcessingError: If files cannot be read or written

    Example:
        >>> from lib.fasta import rename_fasta_gzip
        >>> mapping = {'seq1': 'renamed_seq1', 'seq2': 'renamed_seq2'}
        >>> rename_fasta_gzip('input.fasta', 'output.fasta', mapping)
    """
    input_path = Path(input_file)
    output_path = Path(output_file)

    if not input_path.exists():
        raise FileProcessingError(
            f"Input file not found: {input_path}",
            filepath=str(input_path),
            suggestion="Check that input file exists and is readable",
        )

    try:
        # Create output directory if needed
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with ExitStack() as stack:
            handle_in: IO[str]
            handle_out: IO[str]
            if input_path.suffix == ".gz":
                handle_in = stack.enter_context(gzip.open(input_path, "rt"))
            else:
                handle_in = stack.enter_context(input_path.open("r"))

            if output_path.suffix == ".gz":
                handle_out = stack.enter_context(gzip.open(output_path, "wt"))
            else:
                handle_out = stack.enter_context(output_path.open("w"))

            # Parse FASTA file
            for record in SeqIO.parse(handle_in, "fasta"):
                # Rename sequence if mapping exists
                if record.id in name_mapping:
                    record.id = name_mapping[record.id]
                    record.description = name_mapping[record.id]

                handle_out.write(record.format("fasta"))

    except Exception as exc:
        raise FileProcessingError(
            f"Error processing {input_path}: {exc}",
            filepath=str(input_path),
            suggestion="Check file format and permissions",
        ) from exc
