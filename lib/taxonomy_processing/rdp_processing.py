"""
RDP classifier processing utilities.

This module provides functions for processing RDP classifier outputs.
"""

from pathlib import Path
from typing import Union

import pandas as pd

from ..exceptions import FileProcessingError, ValidationError


def tsv_to_csv(
    input_file: Union[str, Path], marker: str, output_file: Union[str, Path] = None
) -> pd.DataFrame:
    """
    Convert RDP classifier TSV output to CSV with proper headers.

    Args:
        input_file: Path to RDP classifier output (TSV format)
        marker: Marker type (determines column headers)
        output_file: Optional path to write CSV (if not provided, returns DataFrame)

    Returns:
        DataFrame with proper column headers (if output_file is None)

    Raises:
        FileProcessingError: If file cannot be read
        ValidationError: If marker is invalid

    Example:
        >>> from lib.taxonomy_processing import tsv_to_csv
        >>> df = tsv_to_csv('rdp.out.tmp', '16S')
        >>> df.to_csv('rdp.csv', index=False)

    Author: Teresita M. Porter, August 20, 2022
    """
    input_path = Path(input_file)

    if not input_path.exists():
        raise FileProcessingError(
            f"RDP output file not found: {input_path}",
            filepath=str(input_path),
            suggestion="Check that RDP output file exists and is readable",
        )

    # Read RDP output
    try:
        df = pd.read_csv(input_path, sep="\t", header=None)
    except Exception as exc:
        raise FileProcessingError(
            f"Failed to read RDP output: {input_path}",
            filepath=str(input_path),
            suggestion="Ensure file is valid TSV format",
        ) from exc

    # Define marker types and their column headers
    tax3_abund6 = ["COI", "rbcL_eukaryota", "rbcL_landPlant", "12S_fish", "12S_vertebrate"]
    tax3_abund7 = ["16S", "28S_fungi"]
    tax3_abund9 = ["18S_eukaryota"]

    # Add correct headers based on marker
    if marker in tax3_abund6:
        df.columns = [
            "GlobalESV",
            "Strand",
            "Root",
            "RootRank",
            "rBP",
            "SuperKingdom",
            "SuperKingdomRank",
            "skBP",
            "Kingdom",
            "KingdomRank",
            "kBP",
            "Phylum",
            "PhylumRank",
            "pBP",
            "Class",
            "ClassRank",
            "cBP",
            "Order",
            "OrderRank",
            "oBP",
            "Family",
            "FamilyRank",
            "fBP",
            "Genus",
            "GenusRank",
            "gBP",
            "Species",
            "SpeciesRank",
            "sBP",
        ]
    elif marker in tax3_abund7:
        df.columns = [
            "GlobalESV",
            "Strand",
            "Domain",
            "DomainRank",
            "dBP",
            "Phylum",
            "PhylumRank",
            "pBP",
            "Class",
            "ClassRank",
            "cBP",
            "Order",
            "OrderRank",
            "oBP",
            "Family",
            "FamilyRank",
            "fBP",
            "Genus",
            "GenusRank",
            "gBP",
        ]
    elif marker in tax3_abund9:
        df.columns = [
            "GlobalESV",
            "Strand",
            "Root",
            "RootRank",
            "rBP",
            "Domain",
            "DomainRank",
            "dBP",
            "Kingdom",
            "KingdomRank",
            "kBP",
            "Phylum",
            "PhylumRank",
            "pBP",
            "Class",
            "ClassRank",
            "cBP",
            "Order",
            "OrderRank",
            "oBP",
            "Family",
            "FamilyRank",
            "fBP",
            "Genus",
            "GenusRank",
            "gBP",
        ]
    else:
        raise ValidationError(
            f"Invalid marker type: {marker}",
            suggestion=f"Valid markers are: {', '.join(tax3_abund6 + tax3_abund7 + tax3_abund9)}",
            details=f"Got marker: {marker}",
        )

    # Write to output path if provided
    if output_file:
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_path, index=False, header=True)

    return df


def filter_rdp_taxonomy(
    input_file: Union[str, Path], min_confidence: float = 0.8, output_file: Union[str, Path] = None
) -> pd.DataFrame:
    """
    Filter RDP taxonomy output by confidence threshold.

    Args:
        input_file: Path to RDP taxonomy output
        min_confidence: Minimum confidence threshold (0.0-1.0)
        output_file: Optional path to write filtered output (if not provided, returns DataFrame)

    Returns:
        Filtered DataFrame (if output_file is None)

    Raises:
        FileProcessingError: If file cannot be read
        ValidationError: If confidence threshold is invalid

    Example:
        >>> from lib.taxonomy_processing import filter_rdp_taxonomy
        >>> df = filter_rdp_taxonomy('rdp_out.tsv', 0.8)
    """
    input_path = Path(input_file)

    if not input_path.exists():
        raise FileProcessingError(
            f"RDP taxonomy file not found: {input_path}",
            filepath=str(input_path),
            suggestion="Check that RDP taxonomy file exists and is readable",
        )

    if not 0.0 <= min_confidence <= 1.0:
        raise ValidationError(
            f"Invalid confidence threshold: {min_confidence}",
            suggestion="Confidence threshold must be between 0.0 and 1.0",
            details=f"Got confidence: {min_confidence}",
        )

    # Read RDP output
    try:
        df = pd.read_csv(input_path, sep="\t")
    except Exception as exc:
        raise FileProcessingError(
            f"Failed to read RDP taxonomy: {input_path}",
            filepath=str(input_path),
            suggestion="Ensure file is valid TSV format",
        ) from exc

    # Filter by confidence (look for columns ending in 'BP' for confidence scores)
    bp_columns = [col for col in df.columns if col.endswith("BP")]
    if bp_columns:
        # Keep rows where all BP columns meet threshold
        for col in bp_columns:
            df = df[df[col] >= min_confidence]

    # Write to output path if provided
    if output_file:
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_path, sep="\t", index=False)

    return df


def parallel_rdp(
    input_file: Union[str, Path], output_dir: Union[str, Path], num_threads: int = 4
) -> None:
    """
    Run RDP classifier in parallel mode.

    Args:
        input_file: Path to input sequences (FASTA format)
        output_dir: Directory to write output files
        num_threads: Number of threads for parallel processing

    Raises:
        FileProcessingError: If input file not found
        ValidationError: If parameters are invalid

    Example:
        >>> from lib.taxonomy_processing import parallel_rdp
        >>> parallel_rdp('sequences.fasta', 'rdp_output/', num_threads=8)
    """
    input_path = Path(input_file)
    output_path = Path(output_dir)

    if not input_path.exists():
        raise FileProcessingError(
            f"Input file not found: {input_path}",
            filepath=str(input_path),
            suggestion="Check that input file exists and is readable",
        )

    if num_threads < 1:
        raise ValidationError(
            f"Invalid number of threads: {num_threads}",
            suggestion="Number of threads must be >= 1",
            details=f"Got threads: {num_threads}",
        )

    # Create output directory
    output_path.mkdir(parents=True, exist_ok=True)

    # Note: Actual RDP classifier execution would go here
    # This is a placeholder for the implementation
    raise NotImplementedError(
        "parallel_rdp is not yet implemented",
        suggestion="Use the RDP classifier directly or implement this function",
    )
