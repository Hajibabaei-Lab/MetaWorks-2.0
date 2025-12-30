"""
FASTA file processing utilities.

This module provides functionality for manipulating FASTA files.
"""

from .stats import get_fasta_stats, get_fastq_stats
from .utils import rename_fasta_gzip, reverse_complement

__all__ = ["get_fasta_stats", "get_fastq_stats", "reverse_complement", "rename_fasta_gzip"]
