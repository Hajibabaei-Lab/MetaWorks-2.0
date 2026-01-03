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
    classification_config = config.get("modules", {}).get("classification", config.get("classification", {}))
    
    if classification_config.get("use_custom_classifier", False):
        classifier_path = classification_config.get("classifier_path", "classifiers/COI.properties")
        return f"-t {classifier_path}"
    else:
        marker = config.get("marker", classification_config.get("marker", "COI"))
        builtin_classifier = classification_config.get("builtin_classifier", "fungallsu")
        
        if marker == "16S":
            return f"-c {builtin_classifier}"
        elif marker in ["ITS", "28S_fungi", "fungal_ITS"]:
            return f"-g {builtin_classifier}"
        else:
            # For other markers, use default classifier
            return f"-t {builtin_classifier}"

# Input validation
def validate_inputs():
    """Validate required configuration parameters"""
    required = ["dir"]
    for req in required:
        if req not in config:
            raise ValueError(f"Missing required config parameter: {req}")
    
    # If using custom classifier, validate properties file exists
    if RDP_CONFIG.get("use_custom_classifier", False):
        props_file = RDP_CONFIG.get("classifier_path", "classifiers/COI.properties")
        if not os.path.exists(props_file):
            raise FileNotFoundError(f"RDP classifier properties file not found: {props_file}")

validate_inputs()

# Rule for taxonomic assignment using RDP classifier
rule taxonomic_assignment:
    """Perform taxonomic assignment using RDP Classifier with parallel processing"""
    input:
        config["pipeline"]["output_dir"] + "/cat.denoised.nonchimeras"
    output:
        config["pipeline"]["output_dir"] + "/rdp.out.tmp"
    threads: 4
    resources:
        mem_mb = 20000,
        time_min = 240
    log:
        config["pipeline"]["output_dir"] + "/logs/classification.log"
    params:
        memory = lambda wc: f"{config.get('modules', {}).get('classification', config.get('classification', {})).get('memory_gb', 20)}g",
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
        config["pipeline"]["output_dir"] + "/rdp.out.tmp"
    output:
        touch(config["pipeline"]["output_dir"] + "/checkpoints/classification_complete.done")
    message:
        "Taxonomic classification complete for {0} sequences".format(
            count_sequences(config["pipeline"]["output_dir"] + "/cat.denoised.nonchimeras")
        )

def count_sequences(fasta_file):
    """Count the number of sequences in a FASTA file"""
    count = 0
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                count += 1
    return count
