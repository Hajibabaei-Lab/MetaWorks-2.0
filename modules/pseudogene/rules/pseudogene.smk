# rules/pseudogene.smk

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

# Get classification config for marker
CLASSIFICATION_CONFIG = get_module_config(config, "classification")

# Module declarations for ORF processing
module orfs_hmm:
    snakefile: "orfs_hmm.smk"
    config: {"dir": config["pipeline"]["output_dir"]}

module orfs_longest:
    snakefile: "orfs_longest.smk"
    config: {"dir": config["pipeline"]["output_dir"]}

if PSEUDOGENE_CONFIG.get("method", "none") in ["hmm", "orf"]:

    ############################
    # Taxon Subsetting
    ############################

    if PSEUDOGENE_CONFIG.get("grep_type", 1) == 1:
        rule subset_taxonomy_by_taxon1:
            input: config["pipeline"]["output_dir"] + "/rdp.out.tmp"
            output: config["pipeline"]["output_dir"] + "/taxon.zotus"
            params:
                taxon1 = PSEUDOGENE_CONFIG.get("taxon1", "")
            shell:
                "grep {params.taxon1} {input} | awk '{{print $1}}' > {output}"
    else:
        rule subset_taxonomy_by_taxon1_and_taxon2:
            input: config["pipeline"]["output_dir"] + "/rdp.out.tmp"
            output: config["pipeline"]["output_dir"] + "/taxon.zotus"
            params:
                taxon1 = PSEUDOGENE_CONFIG.get("taxon1", ""),
                taxon2 = PSEUDOGENE_CONFIG.get("taxon2", "")
            shell:
                "grep {params.taxon1} {input} | grep {params.taxon2} | awk '{{print $1}}' > {output}"

    rule subset_ESV_sequences_by_taxon:
        input:
            tax = config["pipeline"]["output_dir"] + "/taxon.zotus",
            fas = config["pipeline"]["output_dir"] + "/cat.denoised.nonchimeras"
        output: config["pipeline"]["output_dir"] + "/chimera.denoised.nonchimeras.taxon"
        shell:
            "python3 python_scripts/get_taxon_only.py {input.tax} {input.fas} > {output}"

    ############################
    # Strategy 1: Longest ORFs
    ############################

    if PSEUDOGENE_CONFIG.get("method", "hmm") == "orf":
        use rule get_orfs_nt_longest from orfs_longest as get_orfs_nt
        use rule get_longest_orfs from orfs_longest

        rule filter_rdp:
            input:
                orf = config["pipeline"]["output_dir"] + "/longest.orfs.fasta",
                rdp = config["pipeline"]["output_dir"] + "/rdp.out.tmp"
            output: config["pipeline"]["output_dir"] + "/taxonomy.csv"
            params:
                marker = CLASSIFICATION_CONFIG.get("marker", "COI")
            shell:
                "python3 python_scripts/filter_rdp_taxonomy.py {input.orf} {input.rdp} {params.marker} > {output}"

        rule filter_ESV_table:
            input:
                table = config["pipeline"]["output_dir"] + "/ESV.table.tmp",
                orf = config["pipeline"]["output_dir"] + "/longest.orfs.fasta"
            output: config["pipeline"]["output_dir"] + "/ESV.table"
            shell:
                "python3 python_scripts/filter_ESV_table.py {input.table} {input.orf} > {output}"

        rule add_good_orfs_to_taxonomy:
            input:
                orf = config["pipeline"]["output_dir"] + "/longest.orfs.fasta",
                rdp = config["pipeline"]["output_dir"] + "/rdp.out.tmp"
            output: temp(config["pipeline"]["output_dir"] + "/taxonomy_ORF.tsv")
            shell:
                "python3 python_scripts/add_seqs_to_tax4.py {input.orf} {input.rdp} >> {output}"

    ############################
    # Strategy 2: HMM-based ORFs
    ############################

    elif PSEUDOGENE_CONFIG.get("method", "hmm") == "hmm":
        use rule get_orfs_nt_hmm from orfs_hmm
        use rule get_orfs_aa_hmm from orfs_hmm
        use rule consolidate_orfs_hmm from orfs_hmm as consolidate_orfs

        rule hmmscan:
            input:
                orf = config["pipeline"]["output_dir"] + "/orfs.fasta.aa.filtered.hmm",
                hmm = PSEUDOGENE_CONFIG.get("hmm_profile", "config/hmm/bold.hmm")
            output: config["pipeline"]["output_dir"] + "/hmm.txt"
            shell:
                "hmmscan --tblout {output} {input.hmm} {input.orf}"

        rule add_good_orf_sequences_to_taxonomy:
            input:
                hmmer = config["pipeline"]["output_dir"] + "/hmm.txt",
                rdp = config["pipeline"]["output_dir"] + "/rdp.out.tmp",
                orfs = config["pipeline"]["output_dir"] + "/orfs.fasta.nt.filtered.hmm"
            output: config["pipeline"]["output_dir"] + "/taxonomy_ORF.tsv"
            params:
                marker = CLASSIFICATION_CONFIG.get("marker", "COI")
            shell:
                "python3 python_scripts/filter_rdp.py {input.hmmer} {input.orfs} {input.rdp} {params.marker} >> {output}"

        rule filter_rdp:
            input:
                orf = config["pipeline"]["output_dir"] + "/orfs.fasta.nt.filtered.hmm",
                rdp = config["pipeline"]["output_dir"] + "/rdp.out.tmp"
            output: config["pipeline"]["output_dir"] + "/taxonomy.csv"
            params:
                marker = CLASSIFICATION_CONFIG.get("marker", "COI")
            shell:
                "python3 python_scripts/filter_rdp_taxonomy.py {input.orf} {input.rdp} {params.marker} > {output}"

        rule filter_ESV_table:
            input:
                table = config["pipeline"]["output_dir"] + "/ESV.table.tmp",
                orf = config["pipeline"]["output_dir"] + "/orfs.fasta.nt.filtered.hmm"
            output: config["pipeline"]["output_dir"] + "/ESV.table"
            shell:
                "python3 python_scripts/filter_ESV_table.py {input.table} {input.orf} > {output}"

    else:
        print("ERROR: Invalid or missing 'method' in pseudogene_filtering config. Expected 'hmm' or 'orf'")
