# rules/adapter_trimming.smk

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

TRIMMING_CONFIG = get_module_config(config, "trimming")
PREPROCESSING_CONFIG = get_module_config(config, "preprocessing")

rule pair_reads:
    input:
        f = lambda wildcards: get_fastq_path(wildcards.sample, 1),
        r = lambda wildcards: get_fastq_path(wildcards.sample, 2)
    output:
        X1 = temp("{sample}_R1.out"),
        X2 = temp("{sample}_R2.out"),
        s = config["pipeline"]["output_dir"] + "/paired/{sample}.fastq.gz"
    params:
        q = lambda wc: PREPROCESSING_CONFIG.get("quality_score", 13),
        m = lambda wc: PREPROCESSING_CONFIG.get("min_overlap", 25),
        n = lambda wc: PREPROCESSING_CONFIG.get("max_mismatch", 0.02),
        o = lambda wc: PREPROCESSING_CONFIG.get("min_match", 0.90)
    shell:
        """
        SeqPrep \
            -f {input.f} \
            -r {input.r} \
            -1 {output.X1} \
            -2 {output.X2} \
            -q {params.q} \
            -m {params.m} \
            -n {params.n} \
            -s {output.s} \
            -o {params.o}
        """

rule trim_linked_adapters:
    input:
        adapters = TRIMMING_CONFIG.get("adapters"),
        paired = config["pipeline"]["output_dir"] + "/paired/{sample}.fastq.gz"
    output:
        config["pipeline"]["output_dir"] + "/trimmed/{sample}.fasta"
    params:
        m = lambda wc: TRIMMING_CONFIG.get("min_length", 150),
        q = lambda wc: TRIMMING_CONFIG.get("quality_cutoff", "20,20"),
        e = lambda wc: TRIMMING_CONFIG.get("error_rate", 0.1),
        O = lambda wc: TRIMMING_CONFIG.get("min_adapter_overlap", 3),
        mn = lambda wc: TRIMMING_CONFIG.get("max_n_bases", 3),
        rc = lambda wc: TRIMMING_CONFIG.get("enable_rc", False)
    shell:
        """
        cutadapt \
            -a file:{input.adapters} \
            -m {params.m} \
            -q {params.q} \
            -e {params.e} \
            -O {params.O} \
            --max-n={params.mn} \
            {("--reverse-complement" if {params.rc} else "")} \
            --prefix {{name}} \
            --discard-untrimmed \
            --output {output} \
            {input.paired} \
        """

rule gzip_trimmed_fasta:
    input:
        config["pipeline"]["output_dir"] + "/trimmed/{sample}.fasta"
    output:
        config["pipeline"]["output_dir"] + "/trimmed/{sample}.fasta.gz"
    shell:
        "gzip -c {input} > {output}"
