# Classification Module - Taxonomic Assignment
# Performs taxonomic assignment using RDP Classifier with parallel processing

import os
from pathlib import Path

# Module configuration
MODULE_NAME = "classification"
MODULE_CONFIG = config.get("modules", {}).get(MODULE_NAME, {})

# Use module config if available, otherwise fall back to main config
RDP_CONFIG = MODULE_CONFIG.get("rdp", config.get("RDP", {}))

def rdp_options(config):
    """Generate RDP classifier options based on configuration"""
    if config["RDP"]["custom"] == "yes":
        return f"-t {config['RDP']['t']}"
    elif config["RDP"]["custom"] == "no":
        if config["marker"] == "16S":
            return f"-c {config['RDP']['c']} -f {config['RDP']['f']}"
        elif config["marker"] == "28S_fungi":
            return f"-g {config['RDP']['g']}"
        else:
            print("ERROR: Indicate whether you are analyzing a prokaryote 16S or fungal ITS or 28S (LSU) marker in the config")
            return ""
    else:
        print("ERROR: Indicate whether you are working with a custom reference sequence database in the config")
        return ""

# Input validation
def validate_inputs():
    """Validate required configuration parameters"""
    required = ["dir"]
    for req in required:
        if req not in config:
            raise ValueError(f"Missing required config parameter: {req}")
    
    # Validate RDP custom option
    if RDP_CONFIG.get("custom") not in ["yes", "no"]:
        raise ValueError(f"RDP custom must be 'yes' or 'no', got: {RDP_CONFIG.get('custom')}")
    
    # If using custom classifier, validate properties file exists
    if RDP_CONFIG.get("custom") == "yes":
        props_file = RDP_CONFIG.get("t", "path/to/rRNAClassifier.properties")
        if not os.path.exists(props_file):
            raise FileNotFoundError(f"RDP classifier properties file not found: {props_file}")

validate_inputs()

# Rule for taxonomic assignment using RDP classifier
rule taxonomic_assignment:
    """Perform taxonomic assignment using RDP Classifier with parallel processing"""
    input:
        config["dir"] + "/cat.denoised.nonchimeras"
    output:
        config["dir"] + "/rdp.out.tmp"
    threads: 4
    resources:
        mem_mb = 20000,
        time_min = 240
    log:
        config["dir"] + "/logs/classification.log"
    params:
        memory = RDP_CONFIG.get("memory", "20g"),
        options = lambda wildcards: rdp_options(config)
    shell:
        """
        python3 python_scripts/parallel_rdp.py \
            --input {input} \
            --output {output} \
            --threads {threads} \
            --memory '{params.memory}' \
            --options '{params.options}' \
            2>&1 | tee {log}
        """

# Checkpoint to signal completion of classification
checkpoint classification_complete:
    """Signal that taxonomic classification is complete"""
    input:
        config["dir"] + "/rdp.out.tmp"
    output:
        touch(config["dir"] + "/checkpoints/classification_complete.done")
    message:
        "Taxonomic classification complete for {0} sequences".format(
            count_sequences(config["dir"] + "/cat.denoised.nonchimeras")
        )

def count_sequences(fasta_file):
    """Count the number of sequences in a FASTA file"""
    count = 0
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                count += 1
    return count
