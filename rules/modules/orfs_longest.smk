# rules/orfs_longest.smk

rule get_orfs_nt_longest:
    input: config["dir"] + "/chimera.denoised.nonchimeras.taxon"
    output: config["dir"] + "/orfs.fasta.nt.longest"
    params:
        g = lambda wc: config["ORFFINDER"]["g"],
        s = lambda wc: config["ORFFINDER"]["s"],
        ml = lambda wc: config["ORFFINDER"]["ml"],
        n = lambda wc: config["ORFFINDER"]["n"],
        strand = lambda wc: config["ORFFINDER"]["strand"]
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
    input: config["dir"] + "/orfs.fasta.nt.longest"
    output: config["dir"] + "/longest.orfs.fasta"
    shell:
        "python3 python_scripts/parse_orfs3.py {input} {output}"
