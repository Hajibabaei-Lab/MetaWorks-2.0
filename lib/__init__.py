"""
MetaWorks Python Library

This package contains reusable Python modules for the MetaWorks pipeline.
Organized by functionality to support modular Snakemake workflows.
"""

__version__ = "2.0.0"
__author__ = "Hajibabaei Lab"

# Make subpackages easily importable
from . import preprocessing
from . import taxonomy
from . import pseudogene
from . import utils

__all__ = [
    "preprocessing",
    "taxonomy", 
    "pseudogene",
    "utils"
]
