"""
FASTA/FASTQ file statistics utilities.

This module provides functions to calculate statistics for sequence files.
"""

import gzip
import statistics
from contextlib import ExitStack
from pathlib import Path
from typing import IO, Dict, Union

from ..exceptions import FileProcessingError


def get_fasta_stats(
    input_file: Union[str, Path], is_gzipped: bool = True
) -> Dict[str, Union[int, float, str]]:
    """
    Calculate statistics for a FASTA file.

    Args:
        input_file: Path to FASTA file (can be gzipped)
        is_gzipped: Whether file is gzipped (default: True)

    Returns:
        Dictionary with statistics: input_file, count, min, max, mean, median, mode

    Raises:
        FileProcessingError: If file cannot be read
        ValidationError: If file is invalid

    Example:
        >>> from lib.fasta import get_fasta_stats
        >>> stats = get_fasta_stats('sequences.fasta.gz')
        >>> print(stats)
    """
    input_path = Path(input_file)

    if not input_path.exists():
        raise FileProcessingError(
            f"FASTA file not found: {input_path}",
            filepath=str(input_path),
            suggestion="Check that FASTA file exists and is readable",
        )

    lengths = []

    try:
        with ExitStack() as stack:
            handle: IO[str]
            if is_gzipped or input_path.suffix == ".gz":
                handle = stack.enter_context(gzip.open(input_path, "rt"))
            else:
                handle = stack.enter_context(input_path.open("r"))

            for i, line in enumerate(handle):
                if (i + 1) % 2 == 0:  # Sequence lines (every other line in FASTA)
                    lengths.append(len(line.strip()))

    except Exception as exc:
        raise FileProcessingError(
            f"Error reading {input_path}: {exc}",
            filepath=str(input_path),
            suggestion="Check file permissions and format",
        ) from exc

    if not lengths:
        raise FileProcessingError(
            f"No sequences found in {input_path}",
            filepath=str(input_path),
            suggestion="Ensure file contains valid FASTA sequences",
        )

    # Calculate statistics
    count = len(lengths)
    minimum = min(lengths)
    maximum = max(lengths)
    mean_val = statistics.mean(lengths)
    median_val = statistics.median(lengths)

    mode_val: Union[int, str]
    try:
        mode_val = statistics.mode(lengths)
    except statistics.StatisticsError:
        mode_val = "NA"

    return {
        "input_file": str(input_path),
        "count": count,
        "min": minimum,
        "max": maximum,
        "mean": mean_val,
        "median": median_val,
        "mode": mode_val,
    }


def get_fastq_stats(
    input_file: Union[str, Path], is_gzipped: bool = True
) -> Dict[str, Union[int, float, str]]:
    """
    Calculate statistics for a FASTQ file.

    Args:
        input_file: Path to FASTQ file (can be gzipped)
        is_gzipped: Whether file is gzipped (default: True)

    Returns:
        Dictionary with statistics: input_file, count, min, max, mean, median, mode

    Raises:
        FileProcessingError: If file cannot be read
        ValidationError: If file is invalid

    Example:
        >>> from lib.fasta import get_fastq_stats
        >>> stats = get_fastq_stats('reads.fastq.gz')
        >>> print(stats)
    """
    input_path = Path(input_file)

    if not input_path.exists():
        raise FileProcessingError(
            f"FASTQ file not found: {input_path}",
            filepath=str(input_path),
            suggestion="Check that FASTQ file exists and is readable",
        )

    lengths = []

    try:
        with ExitStack() as stack:
            handle: IO[str]
            if is_gzipped or input_path.suffix == ".gz":
                handle = stack.enter_context(gzip.open(input_path, "rt"))
            else:
                handle = stack.enter_context(input_path.open("r"))

            for i, line in enumerate(handle):
                if (i + 1) % 4 == 2:  # Sequence lines in FASTQ (lines 2, 6, 10, ...)
                    lengths.append(len(line.strip()))

    except Exception as exc:
        raise FileProcessingError(
            f"Error reading {input_path}: {exc}",
            filepath=str(input_path),
            suggestion="Check file permissions and format",
        ) from exc

    if not lengths:
        raise FileProcessingError(
            f"No sequences found in {input_path}",
            filepath=str(input_path),
            suggestion="Ensure file contains valid FASTQ sequences",
        )

    # Calculate statistics
    count = len(lengths)
    minimum = min(lengths)
    maximum = max(lengths)
    mean_val = statistics.mean(lengths)
    median_val = statistics.median(lengths)

    mode_val: Union[int, str]
    try:
        mode_val = statistics.mode(lengths)
    except statistics.StatisticsError:
        mode_val = "NA"

    return {
        "input_file": str(input_path),
        "count": count,
        "min": minimum,
        "max": maximum,
        "mean": mean_val,
        "median": median_val,
        "mode": mode_val,
    }
