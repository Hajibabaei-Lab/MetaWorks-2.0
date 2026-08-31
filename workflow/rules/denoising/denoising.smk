# rules/denoising.smk
# Uses shared helpers from rules/common.smk.

DENOISING_CONFIG = get_module_config(config, "denoising")

rule edit_fasta_header1:
    input: config["pipeline"]["output_dir"] + "/trimmed/{sample}.fasta.gz"
    output: temp(config["pipeline"]["output_dir"] + "/{sample}.prepared.fasta.gz")
    shell:
        "set -euo pipefail; python3 workflow/scripts/prepare_pooled_reads.py {input} --sample-name {wildcards.sample} > {output}"

# Global pooling path
if DENOISING_CONFIG.get("pool_samples", True):

    rule concatenate_for_global_analysis:
        input:
            expand(config["pipeline"]["output_dir"] + "/{sample}.prepared.fasta.gz", sample=SAMPLES_UNIQUE)
        output:
            config["pipeline"]["output_dir"] + "/cat.fasta.gz"
        shell:
            "set -euo pipefail; cat {input} > {output}"

    rule dereplicate:
        input:
            config["pipeline"]["output_dir"] + "/cat.fasta.gz"
        output:
            config["pipeline"]["output_dir"] + "/cat.uniques"
        log:
            config["pipeline"]["output_dir"] + "/dereplication.log"
        shell:
            "set -euo pipefail; vsearch --fastx_uniques {input} --fastaout {output} --sizein --sizeout --log {log}"

    rule denoise:
        input:
            config["pipeline"]["output_dir"] + "/cat.uniques"
        output:
            config["pipeline"]["output_dir"] + "/cat.denoised"
        log:
            config["pipeline"]["output_dir"] + "/denoising.log"
        params:
            minsize = DENOISING_CONFIG.get("min_cluster_size", 8)
        shell:
            "set -euo pipefail; vsearch --cluster_unoise {input} --sizein --sizeout --minsize {params.minsize} --centroids {output} --log {log}"

    rule chimera_removal:
        input:
            config["pipeline"]["output_dir"] + "/cat.denoised"
        output:
            config["pipeline"]["output_dir"] + "/cat.denoised.nonchimeras"
        log:
            config["pipeline"]["output_dir"] + "/chimeraRemoval.log"
        shell:
            "set -euo pipefail; vsearch --uchime3_denovo {input} --sizein --xsize --nonchimeras {output} --relabel 'Zotu' --log {log}"

    rule create_ESV_table:
        input:
            vsearch_global = config["pipeline"]["output_dir"] + "/cat.fasta.gz",
            db = config["pipeline"]["output_dir"] + "/cat.denoised.nonchimeras"
        output:
            temp(config["pipeline"]["output_dir"] + "/ESV.table.tmp")
        threads: DENOISING_CONFIG.get("threads", 4)
        log:
            config["pipeline"]["output_dir"] + "/table.log"
        shell:
            "set -euo pipefail; vsearch --threads {threads} --search_exact {input.vsearch_global} --db {input.db} --otutabout {output} --log {log}"


# Per-sample processing path
else:

    rule dereplicate_single:
        input:
            config["pipeline"]["output_dir"] + "/{sample}.prepared.fasta.gz"
        output:
            temp(config["pipeline"]["output_dir"] + "/{sample}.uniques.tmp")
        shell:
            "set -euo pipefail; vsearch --fastx_uniques {input} --fastaout {output} --sizein --sizeout"

    rule denoise:
        input:
            config["pipeline"]["output_dir"] + "/{sample}.uniques.tmp"
        output:
            temp(config["pipeline"]["output_dir"] + "/{sample}.denoised")
        log:
            config["pipeline"]["output_dir"] + "/logs/{sample}.denoising.log"
        params:
            minsize = DENOISING_CONFIG.get("min_cluster_size", 8)
        shell:
            "set -euo pipefail; vsearch --cluster_unoise {input} --sizein --sizeout --minsize {params.minsize} --centroids {output} --log {log}"

    rule concatenate_for_global_analysis:
        input:
            expand(config["pipeline"]["output_dir"] + "/{sample}.denoised", sample=SAMPLES_UNIQUE)
        output:
            config["pipeline"]["output_dir"] + "/cat.denoised.tmp"
        shell:
            "set -euo pipefail; cat {input} > {output}"

    rule dereplicate:
        input:
            config["pipeline"]["output_dir"] + "/cat.denoised.tmp"
        output:
            config["pipeline"]["output_dir"] + "/cat.uniques"
        log:
            config["pipeline"]["output_dir"] + "/dereplication.log"
        shell:
            "set -euo pipefail; vsearch --fastx_uniques {input} --fastaout {output} --sizein --sizeout --log {log}"

    rule compress:
        input:
            config["pipeline"]["output_dir"] + "/cat.uniques"
        output:
            config["pipeline"]["output_dir"] + "/cat.uniques.gz"
        shell:
            "set -euo pipefail; gzip -c {input} > {output}"

    rule chimera_removal:
        input:
            config["pipeline"]["output_dir"] + "/cat.uniques.gz"
        output:
            config["pipeline"]["output_dir"] + "/cat.denoised.nonchimeras"
        log:
            config["pipeline"]["output_dir"] + "/chimeraRemoval.log"
        shell:
            "set -euo pipefail; vsearch --uchime3_denovo {input} --sizein --xsize --nonchimeras {output} --relabel 'Zotu' --log {log}"

    rule create_ESV_table:
        input:
            vsearch_global = config["pipeline"]["output_dir"] + "/{sample}.prepared.fasta.gz",
            db = config["pipeline"]["output_dir"] + "/cat.denoised.nonchimeras"
        output:
            temp(config["pipeline"]["output_dir"] + "/{sample}.esv.tmp")
        threads: DENOISING_CONFIG.get("threads", 4)
        log:
            config["pipeline"]["output_dir"] + "/{sample}.table.log"
        shell:
            "set -euo pipefail; vsearch --threads {threads} --search_exact {input.vsearch_global} --db {input.db} --otutabout {output} --log {log}"

    rule merge_sample_ESV_tables:
        input:
            esv_tables = expand(config["pipeline"]["output_dir"] + "/{sample}.esv.tmp", sample=SAMPLES_UNIQUE)
        output:
            config["pipeline"]["output_dir"] + "/ESV.table.tmp"
        shell:
            "set -euo pipefail; python3 workflow/scripts/merge_esv_tables.py {input.esv_tables} > {output}"
