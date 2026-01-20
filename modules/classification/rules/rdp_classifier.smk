# Helper function to safely get module config
def get_module_config(config, module_name):
    """Safely get module configuration, handling boolean values"""
    modules = config.get("modules", {})
    if isinstance(modules, dict):
        module_config = modules.get(module_name, {})
        if isinstance(module_config, dict):
            return module_config
    # Fallback to top-level config
    fallback = config.get(module_name, {})
    return fallback if isinstance(fallback, dict) else {}

CLASSIFICATION_CONFIG = get_module_config(config, "classification")

rule taxonomic_assignment:
    input:
        config["pipeline"]["output_dir"] + "/cat.denoised.nonchimeras"
    output:
        config["pipeline"]["output_dir"] + "/rdp.out.tmp"
    threads: 4
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
            --options '{params.options}'
        """
