# rules/results.smk
# Uses shared helpers from modules/common.smk (get_module_config, is_module_enabled).

PSEUDOGENE_CONFIG = get_module_config(config, "pseudogene_filtering")
CLASSIFICATION_CONFIG = get_module_config(config, "classification")
UTILS_CONFIG = get_module_config(config, "utils")

def header_value(config):
    """Get header value from config"""
    if UTILS_CONFIG.get("include_header", True):
        return config.get("header_value", "OTU")
    return ""

# Check if pseudogene filtering is enabled
pseudogene_enabled = (
    is_module_enabled(config, "pseudogene_filtering", default=True)
    and PSEUDOGENE_CONFIG.get("method", "none") in ["hmm", "orf"]
)

if pseudogene_enabled:

    # Strategy 1: Longest ORFs
    if PSEUDOGENE_CONFIG.get("method", "hmm") == "orf":
        rule generate_results_longest_orf:
            input:
                table = config["pipeline"]["output_dir"] + "/ESV.table.tmp",
                rdp = config["pipeline"]["output_dir"] + "/taxonomy_ORF.tsv"
            output: config["pipeline"]["output_dir"] + "/results.csv"
            params:
                header = header_value(config)
            shell:
                "python3 workflow/scripts/add_abundance_to_rdp_out.py {input.table} {input.rdp} '{params.header}' > {output}"

    # Strategy 2: HMM-based ORFs
    elif PSEUDOGENE_CONFIG.get("method", "hmm") == "hmm" and CLASSIFICATION_CONFIG.get("marker", "COI") == "COI":
        rule generate_results_hmm:
            input:
                table = config["pipeline"]["output_dir"] + "/ESV.table.tmp",
                rdp = config["pipeline"]["output_dir"] + "/taxonomy_ORF.tsv"
            output: config["pipeline"]["output_dir"] + "/results.csv"
            params:
                header = header_value(config)
            shell:
                "python3 workflow/scripts/add_abundance_to_rdp_out.py {input.table} {input.rdp} '{params.header}' > {output}"

    else:
        rule summarize_results_pseudogene:
            output: config["pipeline"]["output_dir"] + "/results.csv"
            shell: 'echo "ESV.table, filtered ORFs, taxonomy.csv" > {output}'

else:
    rule add_ESV_sequences_to_taxonomy:
        input:
            rdp = config["pipeline"]["output_dir"] + "/rdp.out.tmp",
            esvs = config["pipeline"]["output_dir"] + "/cat.denoised.nonchimeras"
        output:
            temp(config["pipeline"]["output_dir"] + "/taxonomy_seq.tsv")
        shell:
            "python3 workflow/scripts/add_seqs_to_tax3.py {input.esvs} {input.rdp} >> {output}"

    if UTILS_CONFIG.get("report_type", 1) == 1:
        rule generate_results_basic:
            input:
                table = config["pipeline"]["output_dir"] + "/ESV.table.tmp",
                rdp = config["pipeline"]["output_dir"] + "/taxonomy_seq.tsv"
            output: config["pipeline"]["output_dir"] + "/results.csv"
            params:
                header = header_value(config)
            shell:
                "python3 workflow/scripts/add_abundance_to_rdp_out.py {input.table} {input.rdp} '{params.header}' > {output}"

    else:
        rule summarize_results_basic:
            output: config["pipeline"]["output_dir"] + "/results.csv"
            shell: 'echo "ESV.table, cat.denoised.nonchimeras, taxonomy.csv" > {output}'
