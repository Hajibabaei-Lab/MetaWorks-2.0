# modules/common.smk
# Shared helper functions for all MetaWorks Snakemake modules.
# Include this file FIRST in Snakefile_ESV so every downstream .smk can use these.


def get_module_config(config, module_name):
    """Get module parameters from the top-level config section.

    Module parameters live at config[module_name] (e.g. config["trimming"]),
    NOT under config["modules"][module_name] (which is a boolean toggle).
    """
    section = config.get(module_name, {})
    return section if isinstance(section, dict) else {}


def is_module_enabled(config, module_name, default=True):
    """Check whether a module is enabled via config["modules"][name]."""
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
    """Return the pipeline output directory."""
    return config["pipeline"]["output_dir"]


def get_classification_engine(config):
    """Canonical engine lookup: modules.classification_engine then classification.engine."""
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


def clustering_enabled(config):
    """Check whether OTU clustering is active."""
    return is_module_enabled(config, "clustering", default=False)


def get_sequences_input(config):
    """Return the FASTA file for classification: centroids.fasta if clustering, else denoised."""
    out = config["pipeline"]["output_dir"]
    if clustering_enabled(config):
        return out + "/centroids.fasta"
    return out + "/cat.denoised.nonchimeras"


def get_abundance_table(config):
    """Return the abundance table: OTU.table if clustering, else ESV.table.tmp."""
    out = config["pipeline"]["output_dir"]
    if clustering_enabled(config):
        return out + "/OTU.table"
    return out + "/ESV.table.tmp"
