"""
ESV table merging utilities.

This module provides functions to merge multiple ESV tables.
"""

from pathlib import Path
from typing import List, Optional, Union

import pandas as pd

from ..exceptions import FileProcessingError, ValidationError


def merge_esv_tables(
    table_paths: List[Union[str, Path]], output_path: Optional[Union[str, Path]] = None
) -> pd.DataFrame:
    """
    Merge multiple ESV tables into a single combined table.

    Args:
        table_paths: List of paths to ESV tables in TSV format
        output_path: Optional path to write merged table (if not provided, returns DataFrame)

    Returns:
        Merged ESV table as DataFrame (if output_path is None)

    Raises:
        FileProcessingError: If files cannot be read or are invalid
        ValidationError: If no valid tables are provided

    Example:
        >>> from lib.esv import merge_esv_tables
        >>> df = merge_esv_tables(['table1.tsv', 'table2.tsv'])
        >>> df.to_csv('merged.tsv', sep='\\t')
    """
    # Validate input paths
    if not table_paths:
        raise ValidationError(
            "No ESV tables provided", suggestion="Provide at least one ESV table path"
        )

    path_list: List[Path] = [Path(p) for p in table_paths]

    for table_path in path_list:
        if not table_path.exists():
            raise FileProcessingError(
                f"ESV table file not found: {table_path}",
                filepath=str(table_path),
                suggestion="Check that ESV table file exists and is readable",
            )

    # Read all ESV tables
    tables = []
    for table_path in path_list:
        try:
            df = pd.read_csv(table_path, sep="\t", index_col=0)
            tables.append(df)
        except pd.errors.EmptyDataError as exc:
            raise FileProcessingError(
                f"ESV table is empty: {table_path}",
                filepath=str(table_path),
                suggestion="Check that ESV table contains data",
            ) from exc
        except pd.errors.ParserError as exc:
            raise FileProcessingError(
                f"Failed to parse ESV table: {table_path}",
                filepath=str(table_path),
                suggestion="Ensure file is valid TSV format with proper headers",
            ) from exc
        except Exception as exc:
            raise FileProcessingError(
                f"Failed to read ESV table: {table_path}",
                filepath=str(table_path),
                suggestion="Check file permissions and format",
            ) from exc

    if not tables:
        raise FileProcessingError(
            "No valid ESV tables to merge", suggestion="Check that input files are valid ESV tables"
        )

    # Merge tables along columns, fill NA with 0, convert to int
    try:
        merged = pd.concat(tables, axis=1).fillna(0).astype(int)
    except Exception as exc:
        raise FileProcessingError(
            "Failed to merge ESV tables",
            suggestion="Ensure tables have compatible formats and indices",
            details=str(exc),
        ) from exc

    # Write to output path if provided
    if output_path:
        out_path = Path(output_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        merged.to_csv(out_path, sep="\t")

    return merged
