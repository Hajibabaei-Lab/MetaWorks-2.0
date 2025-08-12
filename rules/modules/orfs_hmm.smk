# rules/orfs_hmm.smk

rule get_orfs_nt_hmm:
    input: config["dir"] + "/chimera.denoised.nonchimeras.taxon"
    output: config["dir"] + "/orfs.fasta.nt.hmm"
    params:
        g=lambda wc: config["ORFFINDER"]["g"],
        s=lambda wc: config["ORFFINDER"]["s"],
        ml=lambda wc: config["ORFFINDER"]["ml"],
        n=lambda wc: config["ORFFINDER"]["n"],
        strand=lambda wc: config["ORFFINDER"]["strand"]
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
    input: config["dir"] + "/chimera.denoised.nonchimeras.taxon"
    output: config["dir"] + "/orfs.fasta.aa.hmm"
    params:
        g=lambda wc: config["ORFFINDER"]["g"],
        s=lambda wc: config["ORFFINDER"]["s"],
        ml=lambda wc: config["ORFFINDER"]["ml"],
        n=lambda wc: config["ORFFINDER"]["n"],
        strand=lambda wc: config["ORFFINDER"]["strand"]
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
        nt = config["dir"] + "/orfs.fasta.nt.hmm",
        aa = config["dir"] + "/orfs.fasta.aa.hmm"
    output:
        nt2 = config["dir"] + "/orfs.fasta.nt.filtered.hmm",
        aa2 = config["dir"] + "/orfs.fasta.aa.filtered.hmm"
    shell:
        "python3 python_scripts/parse_orfs4.py {input.nt} {input.aa} {output.nt2} {output.aa2}"
