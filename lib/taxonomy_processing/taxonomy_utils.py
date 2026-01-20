"""
Taxonomy utility functions.

This module provides utility functions for manipulating taxonomy data.
"""

from pathlib import Path
from typing import Optional, Union

import pandas as pd

from ..exceptions import FileProcessingError, ValidationError


def add_abundance(
    taxonomy_file: Union[str, Path],
    esv_table: Union[str, Path],
    output_file: Optional[Union[str, Path]] = None,
) -> pd.DataFrame:
    """
    Add abundance data from ESV table to taxonomy output.

    Args:
        taxonomy_file: Path to RDP taxonomy output
        esv_table: Path to ESV table with abundance data
        output_file: Optional path to write output (if not provided, returns DataFrame)

    Returns:
        DataFrame with abundance data added (if output_file is None)

    Raises:
        FileProcessingError: If files cannot be read
        ValidationError: If files are invalid

    Example:
        >>> from lib.taxonomy_processing import add_abundance
        >>> df = add_abundance('rdp_out.tsv', 'esv_table.tsv')
        >>> df.to_csv('rdp_with_abundance.tsv', sep='\\t', index=False)
    """
    taxonomy_path = Path(taxonomy_file)
    esv_path = Path(esv_table)

    # Validate input files exist
    if not taxonomy_path.exists():
        raise FileProcessingError(
            f"Taxonomy file not found: {taxonomy_path}",
            filepath=str(taxonomy_path),
            suggestion="Check that taxonomy file exists and is readable",
        )

    if not esv_path.exists():
        raise FileProcessingError(
            f"ESV table not found: {esv_path}",
            filepath=str(esv_path),
            suggestion="Check that ESV table exists and is readable",
        )

    # Read files
    try:
        taxonomy_df = pd.read_csv(taxonomy_path, sep="\t")
        esv_df = pd.read_csv(esv_path, sep="\t")
    except Exception as exc:
        raise FileProcessingError(
            "Failed to read taxonomy or ESV table",
            filepath=str(taxonomy_path),
            suggestion="Ensure files are valid TSV format",
        ) from exc

    # Check for required columns
    if "GlobalESV" not in taxonomy_df.columns:
        raise FileProcessingError(
            "Taxonomy file missing required column: GlobalESV",
            filepath=str(taxonomy_path),
            suggestion="Ensure taxonomy file has 'GlobalESV' column",
        )

    if "#OTU ID" not in esv_df.columns:
        raise FileProcessingError(
            "ESV table missing required column: #OTU ID",
            filepath=str(esv_path),
            suggestion="Ensure ESV table has '#OTU ID' column",
        )

    # Merge abundance data
    try:
        merged = pd.merge(
            taxonomy_df,
            esv_df[
                ["#OTU ID"] + [col for col in esv_df.columns if col not in ["#OTU ID", "taxonomy"]]
            ],
            left_on="GlobalESV",
            right_on="#OTU ID",
            how="left",
        )
    except Exception as exc:
        raise FileProcessingError(
            "Failed to merge taxonomy with ESV table",
            suggestion="Ensure GlobalESV and #OTU ID columns have matching values",
            details=str(exc),
        ) from exc

    # Write to output path if provided
    if output_file:
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        merged.to_csv(output_path, sep="\t", index=False)

    return merged


def get_taxon_only(
    taxonomy_file: Union[str, Path],
    target_rank: str,
    output_file: Optional[Union[str, Path]] = None,
) -> pd.DataFrame:
    """
    Extract only taxonomic information for a specific rank.

    Args:
        taxonomy_file: Path to RDP taxonomy output
        target_rank: Taxonomic rank to extract (e.g., 'Genus', 'Family', 'Order')
        output_file: Optional path to write output (if not provided, returns DataFrame)

    Returns:
        DataFrame with only target rank information (if output_file is None)

    Raises:
        FileProcessingError: If file cannot be read
        ValidationError: If target rank is not found

    Example:
        >>> from lib.taxonomy_processing import get_taxon_only
        >>> df = get_taxon_only('rdp_out.tsv', 'Genus')
        >>> df.to_csv('genus_only.tsv', sep='\\t', index=False)
    """
    taxonomy_path = Path(taxonomy_file)

    if not taxonomy_path.exists():
        raise FileProcessingError(
            f"Taxonomy file not found: {taxonomy_path}",
            filepath=str(taxonomy_path),
            suggestion="Check that taxonomy file exists and is readable",
        )

    # Read taxonomy file
    try:
        df = pd.read_csv(taxonomy_path, sep="\t")
    except Exception as exc:
        raise FileProcessingError(
            f"Failed to read taxonomy file: {taxonomy_path}",
            filepath=str(taxonomy_path),
            suggestion="Ensure file is valid TSV format",
        ) from exc

    # Check if target rank exists
    if target_rank not in df.columns:
        available_ranks = [col for col in df.columns if col not in ["GlobalESV", "Strand"]]
        raise ValidationError(
            f"Target rank not found: {target_rank}",
            suggestion=f"Available ranks: {', '.join(available_ranks)}",
            details=f"Columns: {list(df.columns)}",
        )

    # Extract GlobalESV and target rank columns
    rank_col = target_rank
    rank_bp_col = f"{target_rank}BP"

    columns_to_keep = ["GlobalESV"]
    if rank_col in df.columns:
        columns_to_keep.append(rank_col)
    if rank_bp_col in df.columns:
        columns_to_keep.append(rank_bp_col)

    result = df[columns_to_keep]

    # Write to output path if provided
    if output_file:
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(output_path, sep="\t", index=False)

    return result
