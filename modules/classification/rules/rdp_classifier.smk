# RDP classification backend
# Uses shared helpers from modules/common.smk (get_module_config).
# CLASSIFICATION_CONFIG is set by the parent classifier.smk.

rule taxonomic_assignment:
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
        memory = lambda wc: f"{(CLASSIFICATION_CONFIG.get('rdp', {}) if isinstance(CLASSIFICATION_CONFIG.get('rdp', {}), dict) else {}).get('memory_gb', CLASSIFICATION_CONFIG.get('memory_gb', 20))}G",
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
