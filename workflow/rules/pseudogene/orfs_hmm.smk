# rules/orfs_hmm.smk
# Loaded as a Snakemake module from pseudogene.smk.
# Receives config with "dir" and "pseudogene_filtering" keys.

OUTPUT_DIR = config.get("dir", "")
PSEUDOGENE_CONFIG = config.get("pseudogene_filtering", {})

rule get_orfs_nt_hmm:
    input: OUTPUT_DIR + "/chimera.denoised.nonchimeras.taxon"
    output: OUTPUT_DIR + "/orfs.fasta.nt.hmm"
    params:
        g=lambda wc: PSEUDOGENE_CONFIG.get("genetic_code", 5),
        s=lambda wc: PSEUDOGENE_CONFIG.get("orf_start_codon", 2),
        ml=lambda wc: PSEUDOGENE_CONFIG.get("min_orf_length", 30),
        n=lambda wc: str(PSEUDOGENE_CONFIG.get("ignore_nested_orfs", True)).lower(),
        strand=lambda wc: PSEUDOGENE_CONFIG.get("strand", "plus")
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

rule get_orfs_aa_hmm:
    input: OUTPUT_DIR + "/chimera.denoised.nonchimeras.taxon"
    output: OUTPUT_DIR + "/orfs.fasta.aa.hmm"
    params:
        g=lambda wc: PSEUDOGENE_CONFIG.get("genetic_code", 5),
        s=lambda wc: PSEUDOGENE_CONFIG.get("orf_start_codon", 2),
        ml=lambda wc: PSEUDOGENE_CONFIG.get("min_orf_length", 30),
        n=lambda wc: str(PSEUDOGENE_CONFIG.get("ignore_nested_orfs", True)).lower(),
        strand=lambda wc: PSEUDOGENE_CONFIG.get("strand", "plus")
    shell:
        """
        ORFfinder -in {input} \
            -g {params.g} \
            -s {params.s} \
            -ml {params.ml} \
            -n {params.n} \
            -strand {params.strand} \
            -outfmt 0 > {output}
        """

rule consolidate_orfs_hmm:
    input:
        nt = OUTPUT_DIR + "/orfs.fasta.nt.hmm",
        aa = OUTPUT_DIR + "/orfs.fasta.aa.hmm"
    output:
        nt2 = OUTPUT_DIR + "/orfs.fasta.nt.filtered.hmm",
        aa2 = OUTPUT_DIR + "/orfs.fasta.aa.filtered.hmm"
    shell:
        "python3 workflow/scripts/parse_orfs4.py {input.nt} {input.aa} {output.nt2} {output.aa2}"
