include: "global_common.smk"

if is_module_enabled(config, "global_esv"):

    GLOBAL_ESV_CONFIG = config.get("global_esv", {})
    TRIAL_DIRS = sorted(GLOBAL_ESV_CONFIG.get("trial_dirs", []))
    GLOBAL_ESV_TARGETS = expand("{trial}/global_results.csv", trial=TRIAL_DIRS)

    rule grab_seqs_from_results:
        input:
            "{trial}/results.csv"
        output:
            "{trial}/global_results.fasta"
        log:
            "{trial}/logs/grab_seqs.log"
        shell:
            """
            set -euo pipefail
            python3 workflow/scripts/grab_seqs_from_results.py {input} > {output} 2> {log}
            """

    rule map_to_global_esvs:
        input:
            fasta = "{trial}/global_results.fasta",
            db = GLOBAL_OUTPUT_DIR + "/global_derep.fasta"
        output:
            "{trial}/global.uc"
        threads: GLOBAL_ESV_CONFIG.get("threads", 4)
        log:
            "{trial}/logs/map_to_global_esvs.log"
        shell:
            """
            set -euo pipefail
            vsearch --threads {threads} --usearch_global {input.fasta} --db {input.db} --uc {output} --id 1.0 --log {log}
            """

    rule map_global_esvs_to_results:
        input:
            uc = "{trial}/global.uc",
            res = "{trial}/results.csv"
        output:
            "{trial}/global_results.csv"
        log:
            "{trial}/logs/map_global_esvs.log"
        shell:
            """
            set -euo pipefail
            python3 workflow/scripts/map_global_esvs_to_results.py {input.uc} {input.res} > {output} 2> {log}
            """
