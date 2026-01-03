# Trimming Module - Adapter Trimming
# Trims linked adapters from paired reads using Cutadapt

import os
from pathlib import Path

# Module configuration
MODULE_NAME = "trimming"
MODULE_CONFIG = config.get("modules", {}).get(MODULE_NAME, {})

# Use module config if available, otherwise fall back to main config
CUTADAPT_CONFIG = MODULE_CONFIG.get("cutadapt", config.get("CUTADAPT", {}))

# Input validation
def validate_inputs():
    """Validate required configuration parameters"""
    required = ["dir"]
    for req in required:
        if req not in config:
            raise ValueError(f"Missing required config parameter: {req}")
    
    # Check adapter file exists
    adapter_file = CUTADAPT_CONFIG.get("fasta")
    if not adapter_file:
        raise ValueError("Adapter file (CUTADAPT.fasta) must be specified")
    
    if not os.path.exists(adapter_file):
        raise FileNotFoundError(f"Adapter file not found: {adapter_file}")
    
    # Validate samples list
    if "SAMPLES_UNIQUE" not in config:
        raise ValueError("SAMPLES_UNIQUE not found in config. Did preprocessing module run?")

validate_inputs()

# Get samples from config (set by preprocessing module)
SAMPLES_UNIQUE = config["SAMPLES_UNIQUE"]

# Rules
rule trim_linked_adapters:
    """Trim linked adapters from paired reads using Cutadapt"""
    input:
        adapters = CUTADAPT_CONFIG.get("fasta"),
        paired = config["pipeline"]["output_dir"] + "/paired/{sample}.fastq.gz"
    output:
        config["pipeline"]["output_dir"] + "/trimmed/{sample}.fasta"
    params:
        m = CUTADAPT_CONFIG.get("m", 150),
        q = CUTADAPT_CONFIG.get("q", "20,20"),
        e = CUTADAPT_CONFIG.get("e", 0.1),
        O = CUTADAPT_CONFIG.get("O", 3),
        mn = CUTADAPT_CONFIG.get("mn", 3),
        rc = CUTADAPT_CONFIG.get("rc", "No")
    threads: 1
    resources:
        mem_mb = 4000,
        time_min = 60
    log:
        config["pipeline"]["output_dir"] + "/logs/trimming/{sample}.log"
    benchmark:
        config["pipeline"]["output_dir"] + "/benchmarks/trimming/{sample}.txt"
    shell:
        """
        cutadapt \
            -a file:{input.adapters} \
            -m {params.m} \
            -q {params.q} \
            -e {params.e} \
            -O {params.O} \
            --max-n={params.mn} \
            --prefix {{name}} \
            --discard-untrimmed \
            --output {output} \
            {input.paired} \
            2>&1 | tee {log}
        """

rule gzip_trimmed_fasta:
    """Compress trimmed FASTA files"""
    input:
        config["pipeline"]["output_dir"] + "/trimmed/{sample}.fasta"
    output:
        config["pipeline"]["output_dir"] + "/trimmed/{sample}.fasta.gz"
    threads: 1
    resources:
        mem_mb = 1000,
        time_min = 10
    log:
        config["pipeline"]["output_dir"] + "/logs/trimming/{sample}.gzip.log"
    shell:
        "gzip -c {input} > {output} 2> {log}"

# Checkpoint to signal completion of trimming
checkpoint trimming_complete:
    """Signal that all samples have been trimmed"""
    input:
        expand(config["pipeline"]["output_dir"] + "/trimmed/{sample}.fasta.gz", sample=SAMPLES_UNIQUE)
    output:
        touch(config["pipeline"]["output_dir"] + "/checkpoints/trimming_complete.done")
    message:
        "Adapter trimming complete for all {0} samples".format(len(SAMPLES_UNIQUE))
