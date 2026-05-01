# rules/adapter_trimming.smk
# Uses shared helpers from modules/common.smk (get_module_config).

TRIMMING_CONFIG = get_module_config(config, "trimming")
PREPROCESSING_CONFIG = get_module_config(config, "preprocessing")

rule pair_reads:
    input:
        f = lambda wildcards: get_fastq_path(wildcards.sample, 1),
        r = lambda wildcards: get_fastq_path(wildcards.sample, 2)
    output:
        X1 = temp(config["pipeline"]["output_dir"] + "/{sample}_R1.out"),
        X2 = temp(config["pipeline"]["output_dir"] + "/{sample}_R2.out"),
        s = config["pipeline"]["output_dir"] + "/paired/{sample}.fastq.gz"
    params:
        q = lambda wc: PREPROCESSING_CONFIG.get("quality_score", 13),
        o = lambda wc: PREPROCESSING_CONFIG.get("min_overlap", 25),
        m = lambda wc: PREPROCESSING_CONFIG.get("max_mismatch", 0.02),
        n = lambda wc: PREPROCESSING_CONFIG.get("min_match", 0.90)
    threads: 1
    log:
        config["pipeline"]["output_dir"] + "/logs/pairing/{sample}.log"
    benchmark:
        config["pipeline"]["output_dir"] + "/benchmarks/pairing/{sample}.txt"
    shell:
        """
        set -euo pipefail
        SeqPrep \
            -f {input.f} \
            -r {input.r} \
            -1 {output.X1} \
            -2 {output.X2} \
            -q {params.q} \
            -o {params.o} \
            -m {params.m} \
            -n {params.n} \
            -s {output.s} \
            2> >(tee "{log}" >&2)
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
        rc = lambda wc: "--revcomp" if TRIMMING_CONFIG.get("enable_rc", False) else ""
    threads: 1
    log:
        config["pipeline"]["output_dir"] + "/logs/trimming/{sample}.log"
    benchmark:
        config["pipeline"]["output_dir"] + "/benchmarks/trimming/{sample}.txt"
    shell:
        """
        set -euo pipefail
        cutadapt \
            -a file:{input.adapters} \
            -m {params.m} \
            -q {params.q} \
            -e {params.e} \
            -O {params.O} \
            --max-n={params.mn} \
            {params.rc} \
            --prefix {{name}} \
            --discard-untrimmed \
            --fasta \
            --output {output} \
            {input.paired} \
            2> >(tee "{log}" >&2)
        """

rule gzip_trimmed_fasta:
    input:
        config["pipeline"]["output_dir"] + "/trimmed/{sample}.fasta"
    output:
        config["pipeline"]["output_dir"] + "/trimmed/{sample}.fasta.gz"
    threads: 1
    log:
        config["pipeline"]["output_dir"] + "/logs/trimming/{sample}.gzip.log"
    shell:
        "gzip -c {input} > {output} 2> {log}"
