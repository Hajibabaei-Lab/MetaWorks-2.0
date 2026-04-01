#!/usr/bin/env python3
"""
Filter ESV table by ESV IDs present in ORF FASTA file.

Usage:
    python filter_ESV_table.py <esv_table.tsv> <orfs.fasta>

Arguments:
    esv_table.tsv: ESV table in TSV format (input)
    orfs.fasta: FASTA file with ESV IDs (input)
    Output: Filtered ESV table printed to stdout

Author: Teresita M. Porter, August 18, 2022
"""

import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

# Add lib to path for custom exceptions
sys.path.insert(0, str(Path(__file__).parent.parent))
from lib.exceptions import FileProcessingError, ValidationError


def validate_arguments():
    """Validate command line arguments."""
    if len(sys.argv) != 3:
        raise ValidationError(
            "Incorrect number of arguments",
            suggestion="Usage: python filter_ESV_table.py <esv_table.tsv> <orfs.fasta>",
            details=f"Expected 2 arguments, got {len(sys.argv) - 1}"
        )

    esv_table_path = Path(sys.argv[1])
    orfs_path = Path(sys.argv[2])

    if not esv_table_path.exists():
        raise FileProcessingError(
            f"ESV table file not found: {esv_table_path}",
            filepath=str(esv_table_path),
            suggestion="Check that ESV table file exists and is readable"
        )

    if not orfs_path.exists():
        raise FileProcessingError(
            f"ORF FASTA file not found: {orfs_path}",
            filepath=str(orfs_path),
            suggestion="Check that ORF FASTA file exists and is readable"
        )

    return esv_table_path, orfs_path


def main():
    """Main function to filter ESV table by ORF IDs."""
    try:
        esv_table_path, orfs_path = validate_arguments()
    except ValidationError as exc:
        print(f"Validation Error: {exc.message}", file=sys.stderr)
        sys.exit(1)
    except FileProcessingError as exc:
        print(f"File Error: {exc.message}", file=sys.stderr)
        sys.exit(1)

    # Read ESV table
    try:
        df = pd.read_csv(esv_table_path, sep='\t')
    except pd.errors.EmptyDataError as exc:
        raise FileProcessingError(
            f"ESV table is empty: {esv_table_path}",
            filepath=str(esv_table_path),
            suggestion="Check that ESV table contains data"
        ) from exc
    except pd.errors.ParserError as exc:
        raise FileProcessingError(
            f"Failed to parse ESV table: {esv_table_path}",
            filepath=str(esv_table_path),
            suggestion="Ensure file is valid TSV format with proper headers"
        ) from exc

    # Read ORF FASTA to get good ESV IDs
    try:
        headers = []
        for record in SeqIO.parse(open(orfs_path), 'fasta'):
            headers.append(record.id)
    except IOError as exc:
        raise FileProcessingError(
            f"Failed to read ORF FASTA: {orfs_path}",
            filepath=str(orfs_path),
            suggestion="Check file permissions and ensure file exists"
        ) from exc
    except Exception as exc:
        raise FileProcessingError(
            f"Failed to parse ORF FASTA: {orfs_path}",
            filepath=str(orfs_path),
            suggestion="Ensure file is valid FASTA format"
        ) from exc

    # Filter ESV table by ESV IDs present in ORF file
    if "#OTU ID" not in df.columns:
        raise FileProcessingError(
            "ESV table missing required column: #OTU ID",
            filepath=str(esv_table_path),
            suggestion="Ensure ESV table has '#OTU ID' column"
        )

    df_filtered = df[df["#OTU ID"].isin(headers)]
    print(df_filtered.to_csv(sep='\t', index=False))


if __name__ == "__main__":
    try:
        main()
    except FileProcessingError as exc:
        print(f"File Processing Error: {exc.message}", file=sys.stderr)
        sys.exit(1)
    except Exception as exc:
        print(f"Unexpected Error: {exc}", file=sys.stderr)
        sys.exit(1)
