"""Configuration loading and merging utilities."""

import os
from pathlib import Path
from typing import Any, Dict, Optional

import yaml


class ConfigError(Exception):
    """Base exception for configuration errors."""

    pass


class FileNotFoundError(ConfigError):
    """Raised when a configuration file is not found."""

    pass


class ParseError(ConfigError):
    """Raised when a configuration file cannot be parsed."""

    pass


def load_yaml(file_path: str) -> Dict[str, Any]:
    """
    Load a YAML configuration file.

    Args:
        file_path: Path to the YAML file

    Returns:
        Dictionary containing the parsed YAML content

    Raises:
        FileNotFoundError: If the file doesn't exist
        ParseError: If the file cannot be parsed
    """
    path = Path(file_path)

    if not path.exists():
        raise FileNotFoundError(f"Configuration file not found: {file_path}")

    try:
        with open(path, "r") as f:
            config = yaml.safe_load(f)

        if config is None:
            config = {}

        return config
    except yaml.YAMLError as e:
        raise ParseError(f"Error parsing YAML file {file_path}: {str(e)}")
    except Exception as e:
        raise ParseError(f"Error reading file {file_path}: {str(e)}")


def merge_configs(base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
    """
    Deep merge two configuration dictionaries.

    Override values take precedence over base values.

    Args:
        base: Base configuration dictionary
        override: Override configuration dictionary

    Returns:
        Merged configuration dictionary
    """
    if not isinstance(base, dict) or not isinstance(override, dict):
        return override

    result = base.copy()

    for key, value in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            # Recursively merge nested dictionaries
            result[key] = merge_configs(result[key], value)
        else:
            # Override with new value
            result[key] = value

    return result


def validate_parameter(
    value: Any,
    param_type: str,
    min_value: Optional[float] = None,
    max_value: Optional[float] = None,
    allowed_values: Optional[list] = None,
) -> tuple[bool, Optional[str]]:
    """
    Validate a parameter against rules.

    Args:
        value: Value to validate
        param_type: Expected type ('integer', 'float', 'string', 'boolean')
        min_value: Minimum value (for numeric types)
        max_value: Maximum value (for numeric types)
        allowed_values: List of allowed values (for any type)

    Returns:
        Tuple of (is_valid, error_message)
    """
    # Type validation
    if param_type == "integer":
        if not isinstance(value, int):
            return False, f"Expected integer, got {type(value).__name__}"
    elif param_type == "float":
        if not isinstance(value, (int, float)):
            return False, f"Expected number, got {type(value).__name__}"
    elif param_type == "string":
        if not isinstance(value, str):
            return False, f"Expected string, got {type(value).__name__}"
    elif param_type == "boolean":
        if not isinstance(value, bool):
            return False, f"Expected boolean, got {type(value).__name__}"

    # Range validation (for numeric types)
    if min_value is not None and value < min_value:
        return False, f"Value {value} is below minimum {min_value}"

    if max_value is not None and value > max_value:
        return False, f"Value {value} is above maximum {max_value}"

    # Allowed values validation
    if allowed_values is not None and value not in allowed_values:
        return False, f"Value {value} not in allowed values: {allowed_values}"

    return True, None


def find_config_file(filename: str, search_paths: list[str]) -> Optional[str]:
    """
    Find a configuration file in a list of search paths.

    Args:
        filename: Name of the configuration file
        search_paths: List of directories to search

    Returns:
        Full path to the file if found, None otherwise
    """
    for search_path in search_paths:
        path = Path(search_path) / filename
        if path.exists():
            return str(path.absolute())

    return None


def resolve_path(path: str, base_dir: str) -> str:
    """
    Resolve a path relative to a base directory.

    If the path is absolute, it's returned as-is.
    If the path is relative, it's resolved relative to base_dir.

    Args:
        path: Path to resolve
        base_dir: Base directory for relative paths

    Returns:
        Resolved absolute path
    """
    if os.path.isabs(path):
        return path

    return str(Path(base_dir) / path)
