"""
File utility functions for MetaWorks pipeline
"""

import os
import gzip
from pathlib import Path
from typing import Union, List


def ensure_directory(directory: Union[str, Path]) -> Path:
    """
    Create directory if it doesn't exist.
    
    Args:
        directory: Path to directory
        
    Returns:
        Path object of the directory
    """
    dir_path = Path(directory)
    dir_path.mkdir(parents=True, exist_ok=True)
    return dir_path


def safe_remove(filepath: Union[str, Path]) -> bool:
    """
    Safely remove a file if it exists.
    
    Args:
        filepath: Path to file
        
    Returns:
        True if file was removed, False if it didn't exist
    """
    file_path = Path(filepath)
    if file_path.exists():
        file_path.unlink()
        return True
    return False


def get_file_extension(filepath: Union[str, Path]) -> str:
    """
    Get file extension, handling .gz files specially.
    
    Args:
        filepath: Path to file
        
    Returns:
        File extension (e.g., '.fastq.gz', '.fasta')
    """
    file_path = Path(filepath)
    
    # Handle .gz files
    if file_path.suffix == '.gz':
        # Get the extension before .gz
        stem = file_path.stem
        if '.' in stem:
            return '.' + stem.split('.')[-1] + '.gz'
        return '.gz'
    
    return file_path.suffix


def is_gzipped(filepath: Union[str, Path]) -> bool:
    """
    Check if file is gzipped.
    
    Args:
        filepath: Path to file
        
    Returns:
        True if file is gzipped
    """
    file_path = Path(filepath)
    
    # Check extension
    if file_path.suffix == '.gz':
        return True
    
    # Check magic number (first 2 bytes)
    if file_path.exists():
        try:
            with open(file_path, 'rb') as f:
                return f.read(2) == b'\x1f\x8b'
        except:
            pass
    
    return False


def find_files_by_pattern(directory: Union[str, Path], 
                          pattern: str, 
                          recursive: bool = False) -> List[Path]:
    """
    Find files matching a glob pattern.
    
    Args:
        directory: Directory to search
        pattern: Glob pattern (e.g., '*.fastq.gz')
        recursive: Search recursively if True
        
    Returns:
        List of matching file paths
    """
    dir_path = Path(directory)
    
    if recursive:
        return sorted(dir_path.rglob(pattern))
    else:
        return sorted(dir_path.glob(pattern))


def get_sample_name(filepath: Union[str, Path], 
                    remove_extensions: List[str] = None) -> str:
    """
    Extract sample name from filepath.
    
    Args:
        filepath: Path to file
        remove_extensions: Extensions to remove (default: ['.fastq', '.fq', '.gz'])
        
    Returns:
        Sample name
    """
    if remove_extensions is None:
        remove_extensions = ['.fastq', '.fq', '.fasta', '.fa', '.gz']
    
    file_path = Path(filepath)
    name = file_path.name
    
    # Remove extensions
    for ext in remove_extensions:
        name = name.replace(ext, '')
    
    # Remove common suffixes (R1, R2, etc.)
    name = name.rstrip('_R1').rstrip('_R2')
    name = name.rstrip('_1').rstrip('_2')
    
    return name
