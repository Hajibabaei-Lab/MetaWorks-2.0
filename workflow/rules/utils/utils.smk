# rules/utils.smk

import os
import sys
import glob
import pandas as pd
import re

_scripts_dir = os.path.join(str(REPO_ROOT), "workflow", "scripts")
if _scripts_dir not in sys.path:
    sys.path.insert(0, _scripts_dir)

from marker_defs import MARKER_TO_CONDITION, HEADERS as MARKER_HEADERS, get_condition, get_header  # noqa: E402

# Extracts sample name from FASTQ file name
def extract_sample_name_from_file(filename):
    base = os.path.basename(filename)
    return re.sub(r"_L001_R[12]_001\.fastq\.gz$", "", base)

# Determines the path to the FASTQ file for a given sample and read
def get_fastq_path(sample, read):
    if config["input"]["sample_source"] == "csv":
        return SAMPLE_PATHS[sample].format(read=read)
    elif config["input"]["sample_source"] == "folder":
        return os.path.join(config["input"]["fastq_dir"], f"{sample}_L001_R{read}_001.fastq.gz")
    else:
        raise ValueError("Invalid sample_source config value.")

# Load sample information based on input method
if config["input"]["sample_source"] == "csv":
    samples_df = pd.read_csv(config["input"]["samples_csv"], sep=",")
    SAMPLES_UNIQUE = samples_df["sample"].tolist()
    SAMPLE_PATHS = dict(zip(samples_df["sample"], samples_df["path"]))

elif config["input"]["sample_source"] == "folder":
    FASTQ_FOLDER = config["input"]["fastq_dir"]
    samples = sorted(glob.glob(os.path.join(FASTQ_FOLDER, "*_R1_001.fastq.gz")))
    SAMPLES_UNIQUE = [extract_sample_name_from_file(f) for f in samples]
else:
    raise ValueError("Unknown sample_source: must be 'csv' or 'folder'")

# Output directory for convenience
OUTPUT_DIR = config["pipeline"]["output_dir"]

def rdp_options(config):
    """Generate RDP classifier options from classification.rdp config."""
    classification = config.get("classification", {})
    rdp_cfg = classification.get("rdp", {}) if isinstance(classification.get("rdp", {}), dict) else {}

    use_custom = rdp_cfg.get("use_custom_classifier", True)
    classifier_path = rdp_cfg.get("classifier_path")
    builtin_classifier = rdp_cfg.get("builtin_classifier", "fungallsu")

    if use_custom:
        if not classifier_path:
            raise ValueError(
                "classification.rdp.classifier_path is required when using a custom RDP classifier"
            )
        return f"-t {classifier_path}"
    else:
        if builtin_classifier in ("fungallsu", "fungalits_unite", "fungalits_warcup"):
            return f"-g {builtin_classifier}"
        else:
            raise ValueError(f"Unknown builtin_classifier: {builtin_classifier}")


HEADERS = MARKER_HEADERS

def condition_key(config):
    """Determine condition key based on marker type."""
    marker = config.get("classification", {}).get("marker", "COI")
    try:
        return get_condition(marker)
    except ValueError:
        print(f"Unknown marker '{marker}' for results configuration.")
        return "default"


def header_value(config):
    """Resolve the header string for results based on marker."""
    key = condition_key(config)
    cfg_headers = config.get("HEADER", {}) if isinstance(config, dict) else {}
    return cfg_headers.get(key, HEADERS.get(key, ""))
