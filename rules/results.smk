# rules/results.smk

if config["pseudogene_filtering"] == "yes":

    if config["removal_type"] == 1:
        rule generate_results_longest_orf:
            input:
                table = config["dir"] + "/ESV.table.tmp",
                rdp = config["dir"] + "/taxonomy_ORF.tsv"
            output: config["dir"] + "/results.csv"
            shell:
                "python3 python_scripts/add_abundance_to_rdp_out.py {input.table} {input.rdp} {header} > {output}"

    elif config["removal_type"] == 2 and config["marker"] == "COI":
        rule generate_results_hmm:
            input:
                table = config["dir"] + "/ESV.table.tmp",
                rdp = config["dir"] + "/taxonomy_ORF.tsv"
            output: config["dir"] + "/results.csv"
            params:
                header = condition_key(config)
            shell:
                "python3 python_scripts/add_abundance_to_rdp_out.py {input.table} {input.rdp} '{params.header}' > {output}"

    else:
        rule summarize_results_pseudogene:
            output: config["dir"] + "/results.csv"
            shell: 'echo "ESV.table, filtered ORFs, taxonomy.csv" > {output}'

else:
    rule add_ESV_sequences_to_taxonomy:
        input:
            rdp = config["dir"] + "/rdp.out.tmp",
            esvs = config["dir"] + "/cat.denoised.nonchimeras"
        output:
            temp(config["dir"] + "/taxonomy_seq.tsv")
        shell:
            "python3 python_scripts/add_seqs_to_tax3.plx {input.esvs} {input.rdp} >> {output}"

    if config["report_type"] == 1:
        rule generate_results_basic:
            input:
                table = config["dir"] + "/ESV.table.tmp",
                rdp = config["dir"] + "/taxonomy_seq.tsv"
            output: config["dir"] + "/results.csv"
            params:
                header = condition_key(config)
            shell:
                "python3 python_scripts/add_abundance_to_rdp_out.py {input.table} {input.rdp} '{params.header}' > {output}"

    else:
        rule summarize_results_basic:
            output: config["dir"] + "/results.csv"
            shell: 'echo "ESV.table, cat.denoised.nonchimeras, taxonomy.csv" > {output}'
