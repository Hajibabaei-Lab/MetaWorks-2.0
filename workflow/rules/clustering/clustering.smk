# rules/clustering/clustering.smk
# OTU clustering at 97% similarity using VSEARCH cluster_smallmem.
# Takes denoised non-chimeras and produces OTU centroids + OTU table.
# Depends on denoising module.

CLUSTERING_CONFIG = get_module_config(config, "clustering")

rule otu_clustering:
    input:
        config["pipeline"]["output_dir"] + "/cat.denoised.nonchimeras"
    output:
        config["pipeline"]["output_dir"] + "/centroids.fasta"
    threads: CLUSTERING_CONFIG.get("threads", 4)
    log:
        config["pipeline"]["output_dir"] + "/logs/clustering/otu_clustering.log"
    params:
        cluster_id = CLUSTERING_CONFIG.get("cluster_id", 0.97)
    shell:
        """
        set -euo pipefail
        vsearch --cluster_smallmem {input} --usersort \
            --id {params.cluster_id} \
            --centroids {output}.tmp \
            --threads {threads} \
            --log {log}
        sed -e '/^>/s/^>Zotu/>Otu/g' {output}.tmp > {output}
        rm -f {output}.tmp
        """

rule create_OTU_table:
    input:
        vsearch_global = config["pipeline"]["output_dir"] + "/cat.fasta.gz",
        db = config["pipeline"]["output_dir"] + "/centroids.fasta"
    output:
        config["pipeline"]["output_dir"] + "/OTU.table"
    threads: CLUSTERING_CONFIG.get("threads", 4)
    log:
        config["pipeline"]["output_dir"] + "/logs/clustering/create_otu_table.log"
    params:
        cluster_id = CLUSTERING_CONFIG.get("cluster_id", 0.97)
    shell:
        """
        set -euo pipefail
        vsearch --threads {threads} \
            --usearch_global {input.vsearch_global} \
            --db {input.db} \
            --strand plus \
            --id {params.cluster_id} \
            --otutabout {output} \
            --log {log}
        """
