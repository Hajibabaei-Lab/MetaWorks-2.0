"""
Configuration utility module for MetaWorks.

Canonical config layout:
  - Module toggles (booleans) live under config["modules"][name].
  - Module parameters live under config[name] (top-level section).

These helpers mirror the functions in modules/common.smk (for Snakemake)
and can be imported by any Python code (API, scripts, tests).
"""


def get_param(config, module_name, param_name, default=None):
    """
    Get a single parameter from a module's top-level config section.

    Args:
        config: Full configuration dictionary.
        module_name: Module name (e.g. "trimming", "denoising").
        param_name: Parameter key inside that section.
        default: Fallback if the key is missing.

    Example:
        min_length = get_param(config, "trimming", "min_length", 150)
    """
    section = config.get(module_name, {})
    if isinstance(section, dict):
        return section.get(param_name, default)
    return default


def get_module_config(config, module_name):
    """
    Get the full parameter dict for a module from its top-level config section.

    Returns an empty dict when the section is missing or not a dict.

    Example:
        TRIMMING_CONFIG = get_module_config(config, "trimming")
    """
    section = config.get(module_name, {})
    return section if isinstance(section, dict) else {}


def is_module_enabled(config, module_name, default=True):
    """
    Check whether a module is enabled via config["modules"][name].

    Example:
        if is_module_enabled(config, "stats"):
            ...
    """
    modules = config.get("modules", {})
    if not isinstance(modules, dict):
        return default
    value = modules.get(module_name, default)
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.strip().lower() in ("1", "true", "yes", "on")
    return bool(value)


def get_output_dir(config):
    """
    Get the pipeline output directory from configuration.

    Example:
        output_dir = get_output_dir(config)
    """
    return config.get("pipeline", {}).get("output_dir", "")


def get_classification_engine(config):
    """
    Canonical classification engine lookup.

    Priority: modules.classification_engine > classification.engine > "rdp".
    """
    modules = config.get("modules", {})
    if isinstance(modules, dict):
        engine = modules.get("classification_engine")
        if engine is not None:
            return str(engine).strip().lower()
    classification = config.get("classification", {})
    if isinstance(classification, dict):
        engine = classification.get("engine")
        if engine is not None:
            return str(engine).strip().lower()
    return "rdp"
