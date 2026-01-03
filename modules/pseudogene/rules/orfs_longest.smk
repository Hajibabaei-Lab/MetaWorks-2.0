# rules/orfs_longest.smk

# Get output directory from module config
OUTPUT_DIR = config.get("dir", "")

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

PSEUDOGENE_CONFIG = get_module_config(config, "pseudogene_filtering")

rule get_orfs_nt_longest:
    input: OUTPUT_DIR + "/chimera.denoised.nonchimeras.taxon"
    output: OUTPUT_DIR + "/orfs.fasta.nt.longest"
    params:
        g = lambda wc: PSEUDOGENE_CONFIG.get("genetic_code", 5),
        s = lambda wc: PSEUDOGENE_CONFIG.get("orf_start_codon", 2),
        ml = lambda wc: PSEUDOGENE_CONFIG.get("min_orf_length", 30),
        n = lambda wc: str(PSEUDOGENE_CONFIG.get("ignore_nested_orfs", True)).lower(),
        strand = lambda wc: PSEUDOGENE_CONFIG.get("strand", "plus")
    shell:
        """
        ORFfinder -in {input} \
            -g {params.g} \
            -s {params.s} \
            -ml {params.ml} \
            -n {params.n} \
            -strand {params.strand} \
            -outfmt 1 > {output}
        """

rule get_longest_orfs:
    input: OUTPUT_DIR + "/orfs.fasta.nt.longest"
    output: OUTPUT_DIR + "/longest.orfs.fasta"
    shell:
        "python3 python_scripts/parse_orfs3.py {input} {output}"
