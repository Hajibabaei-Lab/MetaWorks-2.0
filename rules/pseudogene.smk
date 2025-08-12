# rules/pseudogene.smk

if config["pseudogene_filtering"] == 'yes':

    ############################
    # Taxon Subsetting
    ############################

    if config["grep_type"] == 1:
        rule subset_taxonomy_by_taxon1:
            input: config["dir"] + "/rdp.out.tmp"
            output: config["dir"] + "/taxon.zotus"
            shell:
                "grep {config[taxon1]} {input} | awk '{{print $1}}' > {output}"
    else:
        rule subset_taxonomy_by_taxon1_and_taxon2:
            input: config["dir"] + "/rdp.out.tmp"
            output: config["dir"] + "/taxon.zotus"
            shell:
                "grep {config[taxon1]} {input} | grep {config[taxon2]} | awk '{{print $1}}' > {output}"

    rule subset_ESV_sequences_by_taxon:
        input:
            tax = config["dir"] + "/taxon.zotus",
            fas = config["dir"] + "/cat.denoised.nonchimeras"
        output: config["dir"] + "/chimera.denoised.nonchimeras.taxon"
        shell:
            "python3 python_scripts/get_taxon_only.py {input.tax} {input.fas} > {output}"

    ############################
    # Strategy 1: Longest ORFs
    ############################

    if config["removal_type"] == 1:
        use rule get_orfs_nt_longest from orfs_longest as get_orfs_nt
        use rule get_longest_orfs from orfs_longest

        rule filter_rdp:
            input:
                orf = config["dir"] + "/longest.orfs.fasta",
                rdp = config["dir"] + "/rdp.out.tmp"
            output: config["dir"] + "/taxonomy.csv"
            shell:
                "python3 python_scripts/filter_rdp_taxonomy.py {input.orf} {input.rdp} {config[marker]} > {output}"

        rule filter_ESV_table:
            input:
                table = config["dir"] + "/ESV.table.tmp",
                orf = config["dir"] + "/longest.orfs.fasta"
            output: config["dir"] + "/ESV.table"
            shell:
                "python3 python_scripts/filter_ESV_table.py {input.table} {input.orf} > {output}"

        rule add_good_orfs_to_taxonomy:
            input:
                orf = config["dir"] + "/longest.orfs.fasta",
                rdp = config["dir"] + "/rdp.out.tmp"
            output: temp(config["dir"] + "/taxonomy_ORF.tsv")
            shell:
                "python3 python_scripts/add_seqs_to_tax4.py {input.orf} {input.rdp} >> {output}"

    ############################
    # Strategy 2: HMM-based ORFs
    ############################

    elif config["removal_type"] == 2 and config["marker"] == "COI":
        use rule get_orfs_nt_hmm from orfs_hmm
        use rule get_orfs_aa_hmm from orfs_hmm
        use rule consolidate_orfs_hmm from orfs_hmm as consolidate_orfs

        rule hmmscan:
            input:
                orf = config["dir"] + "/orfs.fasta.aa.filtered.hmm",
                hmm = config["hmm"]
            output: config["dir"] + "/hmm.txt"
            shell:
                "hmmscan --tblout {output} {input.hmm} {input.orf}"

        rule add_good_orf_sequences_to_taxonomy:
            input:
                hmmer = config["dir"] + "/hmm.txt",
                rdp = config["dir"] + "/rdp.out.tmp",
                orfs = config["dir"] + "/orfs.fasta.nt.filtered.hmm"
            output: config["dir"] + "/taxonomy_ORF.tsv"
            shell:
                "python3 python_scripts/filter_rdp.py {input.hmmer} {input.orfs} {input.rdp} {config[marker]} >> {output}"

        rule filter_rdp:
            input:
                orf = config["dir"] + "/orfs.fasta.nt.filtered.hmm",
                rdp = config["dir"] + "/rdp.out.tmp"
            output: config["dir"] + "/taxonomy.csv"
            shell:
                "python3 python_scripts/filter_rdp_taxonomy.py {input.orf} {input.rdp} {config[marker]} > {output}"

        rule filter_ESV_table:
            input:
                table = config["dir"] + "/ESV.table.tmp",
                orf = config["dir"] + "/orfs.fasta.nt.filtered.hmm"
            output: config["dir"] + "/ESV.table"
            shell:
                "python3 python_scripts/filter_ESV_table.py {input.table} {input.orf} > {output}"

    else:
        print("ERROR: Invalid or missing 'removal_type' or unsupported 'marker' in config.")
