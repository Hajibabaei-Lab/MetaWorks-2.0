# Preprocessing Module - Read Pairing
# Pairs forward and reverse Illumina reads using SeqPrep

import os
import glob
import pandas as pd
import re
from pathlib import Path

# Module configuration
MODULE_NAME = "preprocessing"
MODULE_CONFIG = config.get("modules", {}).get(MODULE_NAME, {})

# Use module config if available, otherwise fall back to main config
SEQPREP_CONFIG = MODULE_CONFIG.get("seqprep", config.get("SEQPREP", {}))

# Input validation
def validate_inputs():
    """Validate required configuration parameters"""
    required = ["fastq_dir", "dir", "sample_source"]
    for req in required:
        if req not in config:
            raise ValueError(f"Missing required config parameter: {req}")
    
    # Validate sample_source value
    if config["sample_source"] not in ["folder", "csv"]:
        raise ValueError("sample_source must be 'folder' or 'csv'")
    
    # Check FASTQ directory exists
    if not os.path.exists(config["fastq_dir"]):
        raise FileNotFoundError(f"FASTQ directory not found: {config['fastq_dir']}")
    
    # If using CSV, check it exists
    if config["sample_source"] == "csv":
        if "samples_csv" not in config:
            raise ValueError("samples_csv must be specified when sample_source='csv'")
        if not os.path.exists(config["samples_csv"]):
            raise FileNotFoundError(f"Samples CSV not found: {config['samples_csv']}")

validate_inputs()

# Sample detection functions
def extract_sample_name_from_file(filename):
    """Extract sample name from FASTQ filename"""
    base = os.path.basename(filename)
    return re.sub(r"_L001_R[12]_001\.fastq\.gz$", "", base)

def get_fastq_path(sample, read):
    """Get path to FASTQ file for a given sample and read number"""
    if config["sample_source"] == "csv":
        return SAMPLE_PATHS[sample].format(read=read)
    elif config["sample_source"] == "folder":
        return os.path.join(config["fastq_dir"], f"{sample}_L001_R{read}_001.fastq.gz")
    else:
        raise ValueError("Invalid sample_source config value")

# Load sample information
if config["sample_source"] == "csv":
    samples_df = pd.read_csv(config["samples_csv"], sep=",")
    SAMPLES_UNIQUE = samples_df["sample"].tolist()
    SAMPLE_PATHS = dict(zip(samples_df["sample"], samples_df["path"]))
elif config["sample_source"] == "folder":
    FASTQ_FOLDER = config["fastq_dir"]
    samples = glob.glob(os.path.join(FASTQ_FOLDER, "*_R1_001.fastq.gz"))
    SAMPLES_UNIQUE = [extract_sample_name_from_file(f) for f in samples]
    
    if not SAMPLES_UNIQUE:
        raise ValueError(f"No FASTQ files found in {FASTQ_FOLDER} matching pattern *_R1_001.fastq.gz")

# Export for use by other modules
config["SAMPLES_UNIQUE"] = SAMPLES_UNIQUE

# Rules
rule pair_reads:
    """Pair forward and reverse reads using SeqPrep"""
    input:
        f = lambda wildcards: get_fastq_path(wildcards.sample, 1),
        r = lambda wildcards: get_fastq_path(wildcards.sample, 2)
    output:
        unpaired_r1 = temp("{sample}_R1.out"),
        unpaired_r2 = temp("{sample}_R2.out"),
        paired = config["dir"] + "/paired/{sample}.fastq.gz"
    params:
        q = SEQPREP_CONFIG.get("q", 13),
        m = SEQPREP_CONFIG.get("m", 0.02),
        n = SEQPREP_CONFIG.get("n", 0.90),
        o = SEQPREP_CONFIG.get("o", 25)
    threads: 1
    resources:
        mem_mb = 2000,
        time_min = 30
    log:
        config["dir"] + "/logs/pairing/{sample}.log"
    benchmark:
        config["dir"] + "/benchmarks/pairing/{sample}.txt"
    shell:
        """
        SeqPrep \
            -f {input.f} \
            -r {input.r} \
            -1 {output.unpaired_r1} \
            -2 {output.unpaired_r2} \
            -q {params.q} \
            -m {params.m} \
            -n {params.n} \
            -s {output.paired} \
            -o {params.o} \
            2>&1 | tee {log}
        """

# Checkpoint to signal completion of pairing
checkpoint pairing_complete:
    """Signal that all samples have been paired"""
    input:
        expand(config["dir"] + "/paired/{sample}.fastq.gz", sample=SAMPLES_UNIQUE)
    output:
        touch(config["dir"] + "/checkpoints/pairing_complete.done")
    message:
        "Read pairing complete for all {0} samples".format(len(SAMPLES_UNIQUE))
