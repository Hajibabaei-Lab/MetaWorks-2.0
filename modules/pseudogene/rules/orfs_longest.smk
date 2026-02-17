# rules/orfs_longest.smk
# Loaded as a Snakemake module from pseudogene.smk.
# Receives config with "dir" and "pseudogene_filtering" keys.

OUTPUT_DIR = config.get("dir", "")
PSEUDOGENE_CONFIG = config.get("pseudogene_filtering", {})

rule get_orfs_nt_longest:
    input: OUTPUT_DIR + "/chimera.denoised.nonchimeras.taxon"
    output: OUTPUT_DIR + "/orfs.fasta.nt.longest"
    params:
        g = lambda wc: PSEUDOGENE_CONFIG.get("genetic_code", 5),
        s = lambda wc: PSEUDOGENE_CONFIG.get("orf_start_codon", 2),
        ml = lambda wc: PSEUDOGENE_CONFIG.get("min_orf_length", 30),
        n = lambda wc: str(PSEUDOGENE_CONFIG.get("ignore_nested_orfs", True)).lower(),
        strand = lambda wc: PSEUDOGENE_CONFIG.get("strand", "plus")
    threads: 1
    resources:
        mem_mb = 2000,
        time_min = 30
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
    threads: 1
    resources:
        mem_mb = 1000,
        time_min = 10
    shell:
        "python3 python_scripts/parse_orfs3.py {input} {output}"
