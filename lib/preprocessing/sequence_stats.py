"""
Unified sequence statistics module for FASTQ/FASTA files.

This consolidates the functionality from fasta_gz_stats.py and fastq_gz_stats.py
into a single, more maintainable module.
"""

import gzip
import statistics
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Union


class SequenceStats:
    """Calculate statistics for sequence files (FASTQ/FASTA, gzipped or not)"""

    def __init__(self, filepath: Union[str, Path]):
        """
        Initialize SequenceStats.

        Args:
            filepath: Path to sequence file
        """
        self.filepath = Path(filepath)
        self.lengths: List[int] = []
        self.file_type: Optional[str] = None

    def calculate(self) -> Dict[str, Any]:
        """
        Calculate all statistics for the file.

        Returns:
            Dictionary with statistics
        """
        if not self.filepath.exists():
            raise FileNotFoundError(f"File not found: {self.filepath}")

        # Determine if file is gzipped
        if self.filepath.suffix == ".gz" or self._is_gzipped():
            return self._process_gzip()
        return self._process_text()

    def _is_gzipped(self) -> bool:
        """Check if file is gzipped by reading magic number"""
        try:
            with open(self.filepath, "rb") as f:
                return f.read(2) == b"\x1f\x8b"
        except (OSError, IOError):
            return False

    def _process_gzip(self) -> Dict[str, Any]:
        """Process gzipped file"""
        try:
            with gzip.open(self.filepath, "rt") as f:
                self._extract_lengths(f)
        except Exception as e:
            raise RuntimeError(f"Error reading {self.filepath}: {e}")

        return self._compute_stats()

    def _process_text(self) -> Dict[str, Any]:
        """Process uncompressed text file"""
        try:
            with open(self.filepath, "r") as f:
                self._extract_lengths(f)
        except Exception as e:
            raise RuntimeError(f"Error reading {self.filepath}: {e}")

        return self._compute_stats()

    def _extract_lengths(self, file_handle):
        """
        Extract sequence lengths from file.

        Automatically detects FASTQ vs FASTA format.
        """
        first_line = file_handle.readline()
        file_handle.seek(0)

        # Detect file type
        if first_line.startswith("@"):
            self.file_type = "fastq"
            self._extract_fastq_lengths(file_handle)
        elif first_line.startswith(">"):
            self.file_type = "fasta"
            self._extract_fasta_lengths(file_handle)
        else:
            raise ValueError(f"Unknown file format: {self.filepath}")

    def _extract_fastq_lengths(self, file_handle):
        """Extract sequence lengths from FASTQ file"""
        for i, line in enumerate(file_handle):
            # FASTQ: sequence is every 4th line, starting at index 1
            if (i + 1) % 4 == 2:
                self.lengths.append(len(line.strip()))

    def _extract_fasta_lengths(self, file_handle):
        """Extract sequence lengths from FASTA file"""
        for i, line in enumerate(file_handle):
            # FASTA: sequence lines don't start with '>'
            if not line.startswith(">") and line.strip():
                self.lengths.append(len(line.strip()))

    def _compute_stats(self) -> Dict[str, Any]:
        """Compute statistics from collected lengths"""
        if not self.lengths:
            raise ValueError(f"No sequences found in {self.filepath}")

        return {
            "file": str(self.filepath),
            "count": len(self.lengths),
            "min": min(self.lengths),
            "max": max(self.lengths),
            "mean": statistics.mean(self.lengths),
            "median": statistics.median(self.lengths),
            "mode": self._safe_mode(self.lengths),
            "stdev": statistics.stdev(self.lengths) if len(self.lengths) > 1 else 0,
            "file_type": self.file_type,
        }

    @staticmethod
    def _safe_mode(values: List[int]) -> Union[int, str]:
        """Calculate mode with error handling"""
        try:
            return statistics.mode(values)
        except statistics.StatisticsError:
            return "NA"


def calculate_stats(
    filepath: Union[str, Path], output_format: str = "tsv"
) -> Union[str, Dict[str, Any]]:
    """
    Calculate statistics for a sequence file.

    Args:
        filepath: Path to sequence file
        output_format: Output format ('tsv', 'json', or 'dict')

    Returns:
        Formatted statistics string or dict
    """
    stats = SequenceStats(filepath)
    result = stats.calculate()

    if output_format == "dict":
        return result
    elif output_format == "json":
        import json

        return json.dumps(result, indent=2)
    else:  # tsv (default)
        return (
            f"{result['file']}\t{result['count']}\t{result['min']}\t"
            f"{result['max']}\t{result['mean']:.2f}\t{result['median']}\t"
            f"{result['mode']}"
        )


def main():
    """CLI interface for sequence statistics"""
    if len(sys.argv) < 2:
        sys.exit("Usage: python -m lib.preprocessing.sequence_stats <file> [format]")

    filepath = sys.argv[1]
    output_format = sys.argv[2] if len(sys.argv) > 2 else "tsv"

    try:
        result = calculate_stats(filepath, output_format)
        print(result)
    except Exception as e:
        sys.exit(f"Error: {e}")


if __name__ == "__main__":
    main()
