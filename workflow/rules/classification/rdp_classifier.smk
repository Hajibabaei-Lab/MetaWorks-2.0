# RDP classification backend
# Uses shared helpers from rules/common.smk.
# CLASSIFICATION_CONFIG is set by the parent classifier.smk.

rule taxonomic_assignment:
    input:
        get_sequences_input(config)
    output:
        config["pipeline"]["output_dir"] + "/rdp.out.tmp"
    threads: 4
    log:
        config["pipeline"]["output_dir"] + "/logs/classification.log"
    params:
        memory = lambda wc: f"{(CLASSIFICATION_CONFIG.get('rdp', {}) if isinstance(CLASSIFICATION_CONFIG.get('rdp', {}), dict) else {}).get('memory_gb', CLASSIFICATION_CONFIG.get('memory_gb', 20))}G",
        options = lambda wildcards: rdp_options(config)
    shell:
        """
        set -euo pipefail
        python3 workflow/scripts/parallel_rdp.py \
            --input {input} \
            --output {output} \
            --threads {threads} \
            --memory '{params.memory}' \
            --options '{params.options}' \
            2>&1 | tee {log}
        """
