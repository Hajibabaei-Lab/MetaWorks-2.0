#!/usr/bin/env python3
"""
Merge multiple ESV tables into a single combined table.

Usage:
    python merge_esv_tables.py <esv_table1.tsv> <esv_table2.tsv> ...

Arguments:
    One or more ESV tables in TSV format (input)
    Output: Merged ESV table printed to stdout

Author: Alex Song, July 2025
"""

import sys
from pathlib import Path

import pandas as pd

# Add lib to path for custom exceptions
sys.path.insert(0, str(Path(__file__).parent.parent))
from lib.exceptions import FileProcessingError, ValidationError


def validate_arguments():
    """Validate command line arguments."""
    if len(sys.argv) < 2:
        raise ValidationError(
            "Incorrect number of arguments",
            suggestion="Usage: python merge_esv_tables.py <esv_table1.tsv> <esv_table2.tsv> ...",
            details=f"Expected at least 1 argument, got {len(sys.argv) - 1}"
        )
    
    table_paths = [Path(f) for f in sys.argv[1:]]
    
    for table_path in table_paths:
        if not table_path.exists():
            raise FileProcessingError(
                f"ESV table file not found: {table_path}",
                filepath=str(table_path),
                suggestion="Check that ESV table file exists and is readable"
            )
    
    return table_paths


def main():
    """Main function to merge ESV tables."""
    try:
        table_paths = validate_arguments()
    except ValidationError as exc:
        print(f"Validation Error: {exc.message}", file=sys.stderr)
        sys.exit(1)
    except FileProcessingError as exc:
        print(f"File Error: {exc.message}", file=sys.stderr)
        sys.exit(1)
    
    # Read all ESV tables
    tables = []
    for table_path in table_paths:
        try:
            df = pd.read_csv(table_path, sep='\t', index_col=0)
            tables.append(df)
        except pd.errors.EmptyDataError as exc:
            raise FileProcessingError(
                f"ESV table is empty: {table_path}",
                filepath=str(table_path),
                suggestion="Check that ESV table contains data"
            ) from exc
        except pd.errors.ParserError as exc:
            raise FileProcessingError(
                f"Failed to parse ESV table: {table_path}",
                filepath=str(table_path),
                suggestion="Ensure file is valid TSV format with proper headers"
            ) from exc
        except Exception as exc:
            raise FileProcessingError(
                f"Failed to read ESV table: {table_path}",
                filepath=str(table_path),
                suggestion="Check file permissions and format"
            ) from exc
    
    if not tables:
        raise FileProcessingError(
            "No valid ESV tables to merge",
            suggestion="Check that input files are valid ESV tables"
        )
    
    # Merge tables along columns, fill NA with 0, convert to int
    try:
        merged = pd.concat(tables, axis=1).fillna(0).astype(int)
    except Exception as exc:
        raise FileProcessingError(
            "Failed to merge ESV tables",
            suggestion="Ensure tables have compatible formats and indices",
            details=str(exc)
        ) from exc
    
    # Output merged table
    merged.to_csv(sys.stdout, sep='\t')


if __name__ == "__main__":
    try:
        main()
    except FileProcessingError as exc:
        print(f"File Processing Error: {exc.message}", file=sys.stderr)
        sys.exit(1)
    except Exception as exc:
        print(f"Unexpected Error: {exc}", file=sys.stderr)
        sys.exit(1)
