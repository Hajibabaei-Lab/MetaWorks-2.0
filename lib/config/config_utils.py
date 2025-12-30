"""
Configuration utility module for MetaWorks
Provides centralized functions for safely accessing configuration parameters
"""

def get_param(config, module_name, param_name, default=None):
    """
    Safely get a parameter from module configuration.
    
    Handles cases where module config might be a boolean (True/False)
    instead of a dictionary with parameters.
    
    Args:
        config: Full configuration dictionary
        module_name: Name of the module (e.g., "trimming", "denoising")
        param_name: Name of the parameter to retrieve
        default: Default value if parameter not found
    
    Returns:
        Parameter value or default
        
    Example:
        min_length = get_param(config, "trimming", "min_length", 150)
    """
    modules = config.get("modules", {})
    if isinstance(modules, dict):
        module_config = modules.get(module_name, {})
        if isinstance(module_config, dict):
            return module_config.get(param_name, default)
    
    # Fallback to top-level config if modules structure not found
    fallback = config.get(module_name, {})
    if isinstance(fallback, dict):
        return fallback.get(param_name, default)
    
    return default


def get_module_config(config, module_name):
    """
    Safely get entire module configuration dictionary.
    
    Returns empty dict if module config is a boolean or not found.
    
    Args:
        config: Full configuration dictionary
        module_name: Name of the module
    
    Returns:
        Dictionary containing module config, or empty dict
        
    Example:
        TRIMMING_CONFIG = get_module_config(config, "trimming")
    """
    modules = config.get("modules", {})
    if isinstance(modules, dict):
        module_config = modules.get(module_name, {})
        if isinstance(module_config, dict):
            return module_config
    
    # Fallback to top-level config
    fallback = config.get(module_name, {})
    return fallback if isinstance(fallback, dict) else {}


def get_output_dir(config):
    """
    Get the pipeline output directory from configuration.
    
    Args:
        config: Full configuration dictionary
    
    Returns:
        Output directory path as string
        
    Example:
        output_dir = get_output_dir(config)
    """
    pipeline_config = config.get("pipeline", {})
    return pipeline_config.get("output_dir", "")
