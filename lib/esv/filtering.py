"""
ESV table filtering utilities.

This module provides functions to filter ESV tables based on ORF FASTA files.
"""

from pathlib import Path
from typing import Union

import pandas as pd
from Bio import SeqIO

from ..exceptions import FileProcessingError


def filter_esv_table(
    esv_table_path: Union[str, Path],
    orf_fasta_path: Union[str, Path],
    output_path: Union[str, Path] = None,
) -> pd.DataFrame:
    """
    Filter ESV table by ESV IDs present in ORF FASTA file.

    Args:
        esv_table_path: Path to ESV table in TSV format
        orf_fasta_path: Path to ORF FASTA file with ESV IDs
        output_path: Optional path to write filtered table (if not provided, returns DataFrame)

    Returns:
        Filtered ESV table as DataFrame (if output_path is None)

    Raises:
        FileProcessingError: If files cannot be read or are invalid
        ValidationError: If required columns are missing

    Example:
        >>> from lib.esv import filter_esv_table
        >>> df = filter_esv_table('esv_table.tsv', 'orfs.fasta')
        >>> df.to_csv('filtered.tsv', sep='\\t', index=False)
    """
    esv_table_path = Path(esv_table_path)
    orf_fasta_path = Path(orf_fasta_path)

    # Validate input files exist
    if not esv_table_path.exists():
        raise FileProcessingError(
            f"ESV table file not found: {esv_table_path}",
            filepath=str(esv_table_path),
            suggestion="Check that ESV table file exists and is readable",
        )

    if not orf_fasta_path.exists():
        raise FileProcessingError(
            f"ORF FASTA file not found: {orf_fasta_path}",
            filepath=str(orf_fasta_path),
            suggestion="Check that ORF FASTA file exists and is readable",
        )

    # Read ESV table
    try:
        df = pd.read_csv(esv_table_path, sep="\t")
    except pd.errors.EmptyDataError as exc:
        raise FileProcessingError(
            f"ESV table is empty: {esv_table_path}",
            filepath=str(esv_table_path),
            suggestion="Check that ESV table contains data",
        ) from exc
    except pd.errors.ParserError as exc:
        raise FileProcessingError(
            f"Failed to parse ESV table: {esv_table_path}",
            filepath=str(esv_table_path),
            suggestion="Ensure file is valid TSV format with proper headers",
        ) from exc

    # Read ORF FASTA to get ESV IDs
    try:
        esv_ids = []
        for record in SeqIO.parse(orf_fasta_path, "fasta"):
            esv_ids.append(record.id)
    except IOError as exc:
        raise FileProcessingError(
            f"Failed to read ORF FASTA: {orf_fasta_path}",
            filepath=str(orf_fasta_path),
            suggestion="Check file permissions and ensure file exists",
        ) from exc
    except Exception as exc:
        raise FileProcessingError(
            f"Failed to parse ORF FASTA: {orf_fasta_path}",
            filepath=str(orf_fasta_path),
            suggestion="Ensure file is valid FASTA format",
        ) from exc

    # Validate required column exists
    if "#OTU ID" not in df.columns:
        raise FileProcessingError(
            "ESV table missing required column: #OTU ID",
            filepath=str(esv_table_path),
            suggestion="Ensure ESV table has '#OTU ID' column",
        )

    # Filter ESV table by ESV IDs present in ORF file
    df_filtered = df[df["#OTU ID"].isin(esv_ids)]

    # Write to output path if provided
    if output_path:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        df_filtered.to_csv(output_path, sep="\t", index=False)

    return df_filtered
