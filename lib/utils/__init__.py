"""
Utility functions for MetaWorks pipeline
"""

from .file_utils import (
    ensure_directory,
    safe_remove,
    get_file_extension,
    is_gzipped
)

__all__ = [
    "ensure_directory",
    "safe_remove", 
    "get_file_extension",
    "is_gzipped"
]
