# rules/results.smk
# Uses shared helpers from rules/common.smk and lib.config.module_registry.

PSEUDOGENE_CONFIG = get_module_config(config, "pseudogene_filtering")
CLASSIFICATION_CONFIG = get_module_config(config, "classification")
ITSX_CONFIG = get_module_config(config, "itsx_extraction")

pseudogene_enabled = (
    is_module_enabled(config, "pseudogene_filtering")
    and PSEUDOGENE_CONFIG.get("method", "none") in ["hmm", "orf"]
)

itsx_active = is_module_enabled(config, "itsx_extraction")

ITSX_PART = ITSX_CONFIG.get("its_part", "ITS2")

if pseudogene_enabled:

    # Strategy 1: Longest ORFs
    if PSEUDOGENE_CONFIG.get("method", "hmm") == "orf":
        rule generate_results_longest_orf:
            input:
                table = config["pipeline"]["output_dir"] + "/ESV.table",
                rdp = config["pipeline"]["output_dir"] + "/taxonomy_ORF.tsv"
            output: config["pipeline"]["output_dir"] + "/results.csv"
            params:
                header = header_value(config)
            log: config["pipeline"]["output_dir"] + "/logs/results/generate_results.log"
            shell:
                "set -euo pipefail; python3 workflow/scripts/add_abundance_to_rdp_out.py {input.table} {input.rdp} '{params.header}' > {output}"

    # Strategy 2: HMM-based ORFs
    elif PSEUDOGENE_CONFIG.get("method", "hmm") == "hmm":
        rule generate_results_hmm:
            input:
                table = config["pipeline"]["output_dir"] + "/ESV.table",
                rdp = config["pipeline"]["output_dir"] + "/taxonomy_ORF.tsv"
            output: config["pipeline"]["output_dir"] + "/results.csv"
            params:
                header = header_value(config)
            log: config["pipeline"]["output_dir"] + "/logs/results/generate_results.log"
            shell:
                "set -euo pipefail; python3 workflow/scripts/add_abundance_to_rdp_out.py {input.table} {input.rdp} '{params.header}' > {output}"

    else:
        rule summarize_results_pseudogene:
            output: config["pipeline"]["output_dir"] + "/results.csv"
            log: config["pipeline"]["output_dir"] + "/logs/results/summarize.log"
            shell: 'set -euo pipefail; echo "ESV.table, filtered ORFs, taxonomy.csv" > {output}'

elif itsx_active:

    itsx_stripped = config["pipeline"]["output_dir"] + f"/ITSx_out.{ITSX_PART}.fasta.stripped"

    rule filter_itsx_taxonomy:
        input:
            its = itsx_stripped,
            rdp = config["pipeline"]["output_dir"] + "/rdp.out.tmp"
        output:
            config["pipeline"]["output_dir"] + "/taxonomy.csv"
        params:
            marker = CLASSIFICATION_CONFIG.get("marker", "COI")
        log:
            config["pipeline"]["output_dir"] + "/logs/results/filter_taxonomy.log"
        shell:
            "set -euo pipefail; python3 workflow/scripts/filter_rdp_taxonomy.py {input.its} {input.rdp} {params.marker} > {output}"

    rule filter_itsx_esv_table:
        input:
            table = get_abundance_table(config),
            its = itsx_stripped
        output:
            config["pipeline"]["output_dir"] + "/ESV.table"
        log:
            config["pipeline"]["output_dir"] + "/logs/results/filter_esv_table.log"
        shell:
            "set -euo pipefail; python3 workflow/scripts/filter_ESV_table.py {input.table} {input.its} > {output}"

    rule add_itsx_sequences_to_taxonomy:
        input:
            rdp = config["pipeline"]["output_dir"] + "/rdp.out.tmp",
            esvs = itsx_stripped
        output:
            temp(config["pipeline"]["output_dir"] + "/taxonomy_seq.tsv")
        log:
            config["pipeline"]["output_dir"] + "/logs/results/add_seqs_taxonomy.log"
        shell:
            "set -euo pipefail; python3 workflow/scripts/add_seqs_to_tax3.py {input.esvs} {input.rdp} > {output}"

    if config.get("output", {}).get("report_type", 1) == 1:
        rule generate_results_itsx:
            input:
                table = config["pipeline"]["output_dir"] + "/ESV.table",
                rdp = config["pipeline"]["output_dir"] + "/taxonomy_seq.tsv"
            output:
                config["pipeline"]["output_dir"] + "/results.csv"
            params:
                header = header_value(config)
            log:
                config["pipeline"]["output_dir"] + "/logs/results/generate_results.log"
            shell:
                "set -euo pipefail; python3 workflow/scripts/add_abundance_to_rdp_out.py {input.table} {input.rdp} '{params.header}' > {output}"

    else:
        rule summarize_results_itsx:
            output:
                config["pipeline"]["output_dir"] + "/results.csv"
            log:
                config["pipeline"]["output_dir"] + "/logs/results/summarize.log"
            shell: 'set -euo pipefail; echo "ITSx filtered ESV.table, taxonomy.csv" > {output}'

else:
    rule add_ESV_sequences_to_taxonomy:
        input:
            rdp = config["pipeline"]["output_dir"] + "/rdp.out.tmp",
            esvs = get_sequences_input(config)
        output:
            temp(config["pipeline"]["output_dir"] + "/taxonomy_seq.tsv")
        log: config["pipeline"]["output_dir"] + "/logs/results/add_seqs_taxonomy.log"
        shell:
            "set -euo pipefail; python3 workflow/scripts/add_seqs_to_tax3.py {input.esvs} {input.rdp} > {output}"

    if config.get("output", {}).get("report_type", 1) == 1:
        rule generate_results_basic:
            input:
                table = get_abundance_table(config),
                rdp = config["pipeline"]["output_dir"] + "/taxonomy_seq.tsv"
            output: config["pipeline"]["output_dir"] + "/results.csv"
            params:
                header = header_value(config)
            log: config["pipeline"]["output_dir"] + "/logs/results/generate_results.log"
            shell:
                "set -euo pipefail; python3 workflow/scripts/add_abundance_to_rdp_out.py {input.table} {input.rdp} '{params.header}' > {output}"

    else:
        rule summarize_results_basic:
            output: config["pipeline"]["output_dir"] + "/results.csv"
            log: config["pipeline"]["output_dir"] + "/logs/results/summarize.log"
            shell: 'set -euo pipefail; echo "ESV.table, cat.denoised.nonchimeras, taxonomy.csv" > {output}'
