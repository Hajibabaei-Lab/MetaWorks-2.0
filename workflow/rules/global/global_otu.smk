include: "global_common.smk"

if is_module_enabled(config, "global_otu"):

    GLOBAL_OTU_CONFIG = config.get("global_otu", {})
    TRIAL_DIRS = sorted(GLOBAL_OTU_CONFIG.get("trial_dirs", []))
    GLOBAL_OTU_TARGETS = expand("{trial}/global_results.csv", trial=TRIAL_DIRS)

    rule cluster_global_otus:
        input:
            GLOBAL_OUTPUT_DIR + "/global_derep.fasta"
        output:
            GLOBAL_OUTPUT_DIR + "/global_otus.fasta"
        threads: GLOBAL_OTU_CONFIG.get("threads", 4)
        params:
            cluster_id=GLOBAL_OTU_CONFIG.get("cluster_id", 0.97)
        log:
            GLOBAL_OUTPUT_DIR + "/logs/cluster_global_otus.log"
        shell:
            """
            set -euo pipefail
            vsearch --threads {threads} --cluster_smallmem {input} --usersort --id {params.cluster_id} --centroids {output} --log {log}; \
            sed -i 's/^>Zotu/>Otu/g' {output}
            """

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

    rule map_to_global_otus:
        input:
            fasta="{trial}/global_results.fasta",
            db=GLOBAL_OUTPUT_DIR + "/global_otus.fasta"
        output:
            "{trial}/global.uc"
        threads: GLOBAL_OTU_CONFIG.get("threads", 4)
        params:
            map_id=GLOBAL_OTU_CONFIG.get("cluster_id", 0.97)
        log:
            "{trial}/logs/map_to_global_otus.log"
        shell:
            """
            set -euo pipefail
            vsearch --threads {threads} --usearch_global {input.fasta} --db {input.db} --uc {output} --id {params.map_id} --log {log}
            """

    rule map_global_otus_to_results:
        input:
            uc="{trial}/global.uc",
            res="{trial}/results.csv"
        output:
            "{trial}/global_results.csv"
        log:
            "{trial}/logs/map_global_otus.log"
        shell:
            """
            set -euo pipefail
            python3 workflow/scripts/map_global_otus_to_results.py {input.uc} {input.res} > {output} 2> {log}
            """
