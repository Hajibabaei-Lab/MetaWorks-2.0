# VSEARCH SINTAX classification backend
#
# Produces MetaWorks' RDP-like `rdp.out.tmp` so downstream modules (pseudogene, results)
# continue to work unchanged.


def get_module_config(config, module_name):
    """Safely get module configuration, handling boolean values."""
    modules = config.get("modules", {})
    if isinstance(modules, dict):
        module_config = modules.get(module_name, {})
        if isinstance(module_config, dict):
            return module_config
    fallback = config.get(module_name, {})
    return fallback if isinstance(fallback, dict) else {}


CLASSIFICATION_CONFIG = get_module_config(config, "classification")
SINTAX_CONFIG = CLASSIFICATION_CONFIG.get("sintax", {}) if isinstance(CLASSIFICATION_CONFIG.get("sintax", {}), dict) else {}

db_fasta = SINTAX_CONFIG.get("db_fasta")
if not db_fasta:
    raise ValueError(
        "modules.classification_engine is 'sintax' but classification.sintax.db_fasta is not set"
    )

cutoff = SINTAX_CONFIG.get("cutoff")
if cutoff is None:
    cutoff = CLASSIFICATION_CONFIG.get("min_confidence", 0.8)

threads = int(SINTAX_CONFIG.get("threads", 4))
strand = SINTAX_CONFIG.get("strand", "both")


rule taxonomic_assignment:
    input:
        config["pipeline"]["output_dir"] + "/cat.denoised.nonchimeras"
    output:
        raw = temp(config["pipeline"]["output_dir"] + "/sintax.out.tmp"),
        rdp = config["pipeline"]["output_dir"] + "/rdp.out.tmp"
    threads: threads
    params:
        db = db_fasta,
        cutoff = cutoff,
        strand = strand,
        marker = lambda wc: CLASSIFICATION_CONFIG.get("marker", "COI"),
    shell:
        r"""
        vsearch \
            --sintax {input} \
            --db {params.db} \
            --tabbedout {output.raw} \
            --sintax_cutoff {params.cutoff} \
            --strand {params.strand} \
            --threads {threads}

        python3 python_scripts/sintax_to_rdp_out.py \
            --input {output.raw} \
            --marker {params.marker} \
            --output {output.rdp}
        """
