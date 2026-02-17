# rules/denoising.smk
# Uses shared helpers from modules/common.smk (get_module_config).

DENOISING_CONFIG = get_module_config(config, "denoising")

rule edit_fasta_header1:
    input: config["pipeline"]["output_dir"] + "/trimmed/{sample}.fasta.gz"
    output: temp(config["pipeline"]["output_dir"] + "/{sample}.fasta.tmp")
    shell: "python3 python_scripts/rename_fasta_gzip.py {input} > {output}"

# Global pooling path
if DENOISING_CONFIG.get("pool_samples", True):

    rule concatenate_for_global_analysis:
        input:
            expand(config["pipeline"]["output_dir"] + "/{sample}.fasta.tmp", sample=SAMPLES_UNIQUE)
        output:
            temp(config["pipeline"]["output_dir"] + "/cat.fasta.tmp")
        threads: 1
        resources:
            mem_mb = 2000,
            time_min = 10
        shell:
            "cat {input} > {output}"

    rule edit_fasta_header2:
        input:
            config["pipeline"]["output_dir"] + "/cat.fasta.tmp"
        output:
            temp(config["pipeline"]["output_dir"] + "/cat.fasta")
        shell:
            "sed -e 's/-/_/g' {input} > {output}"

    rule compress:
        input:
            config["pipeline"]["output_dir"] + "/cat.fasta"
        output:
            config["pipeline"]["output_dir"] + "/cat.fasta.gz"
        threads: 1
        resources:
            mem_mb = 1000,
            time_min = 5
        shell:
            "gzip -c {input} > {output}"

    rule dereplicate:
        input:
            config["pipeline"]["output_dir"] + "/cat.fasta.gz"
        output:
            config["pipeline"]["output_dir"] + "/cat.uniques"
        threads: 1
        resources:
            mem_mb = 4000,
            time_min = 30
        log:
            config["pipeline"]["output_dir"] + "/dereplication.log"
        shell:
            "vsearch --fastx_uniques {input} --fastaout {output} --sizein --sizeout --log {log}"

    rule denoise:
        input:
            config["pipeline"]["output_dir"] + "/cat.uniques"
        output:
            config["pipeline"]["output_dir"] + "/cat.denoised"
        threads: 1
        resources:
            mem_mb = 4000,
            time_min = 60
        log:
            config["pipeline"]["output_dir"] + "/denoising.log"
        params:
            minsize = DENOISING_CONFIG.get("min_cluster_size", 8)
        shell:
            "vsearch --cluster_unoise {input} --sizein --sizeout --minsize {params.minsize} --centroids {output} --log {log}"

    rule chimera_removal:
        input:
            config["pipeline"]["output_dir"] + "/cat.denoised"
        output:
            config["pipeline"]["output_dir"] + "/cat.denoised.nonchimeras"
        threads: 1
        resources:
            mem_mb = 4000,
            time_min = 45
        log:
            config["pipeline"]["output_dir"] + "/chimeraRemoval.log"
        shell:
            "vsearch --uchime3_denovo {input} --sizein --xsize --nonchimeras {output} --relabel 'Zotu' --log {log}"

    rule create_ESV_table:
        input:
            vsearch_global = config["pipeline"]["output_dir"] + "/cat.fasta.gz",
            db = config["pipeline"]["output_dir"] + "/cat.denoised.nonchimeras"
        output:
            temp(config["pipeline"]["output_dir"] + "/ESV.table.tmp")
        threads: DENOISING_CONFIG.get("threads", 4)
        resources:
            mem_mb = 6000,
            time_min = 60
        log:
            config["pipeline"]["output_dir"] + "/table.log"
        shell:
            "vsearch --threads {threads} --search_exact {input.vsearch_global} --db {input.db} --otutabout {output} --log {log}"


# Per-sample processing path
else:

    rule edit_fasta_header2:
        input:
            config["pipeline"]["output_dir"] + "/{sample}.fasta.tmp"
        output:
            temp(config["pipeline"]["output_dir"] + "/{sample}.tagged")
        shell:
            "sed -e 's/-/_/g' {input} >> {output}"

    rule dereplicate_single:
        input:
            config["pipeline"]["output_dir"] + "/{sample}.tagged"
        output:
            temp(config["pipeline"]["output_dir"] + "/{sample}.uniques.tmp")
        threads: 4
        resources:
            mem_mb = 2000,
            time_min = 20
        shell:
            "vsearch --fastx_uniques {input} --fastaout {output} --sizein --sizeout"

    rule denoise:
        input:
            config["pipeline"]["output_dir"] + "/{sample}.uniques.tmp"
        output:
            temp(config["pipeline"]["output_dir"] + "/{sample}.denoised")
        threads: 4
        resources:
            mem_mb = 2000,
            time_min = 30
        log:
            config["pipeline"]["output_dir"] + "/logs/{sample}.denoising.log"
        params:
            minsize = DENOISING_CONFIG.get("min_cluster_size", 8)
        shell:
            "vsearch --cluster_unoise {input} --sizein --sizeout --minsize {params.minsize} --centroids {output} --log {log}"

    rule concatenate_for_global_analysis:
        input:
            expand(config["pipeline"]["output_dir"] + "/{sample}.denoised", sample=SAMPLES_UNIQUE)
        output:
            config["pipeline"]["output_dir"] + "/cat.denoised.tmp"
        threads: 4
        resources:
            mem_mb = 2000,
            time_min = 10
        shell:
            "cat {input} > {output}"

    rule dereplicate:
        input:
            config["pipeline"]["output_dir"] + "/cat.denoised.tmp"
        output:
            config["pipeline"]["output_dir"] + "/cat.uniques"
        threads: 4
        resources:
            mem_mb = 4000,
            time_min = 30
        log:
            config["pipeline"]["output_dir"] + "/dereplication.log"
        shell:
            "vsearch --fastx_uniques {input} --fastaout {output} --sizein --sizeout --log {log}"

    rule compress:
        input:
            config["pipeline"]["output_dir"] + "/cat.uniques"
        output:
            config["pipeline"]["output_dir"] + "/cat.uniques.gz"
        threads: 1
        resources:
            mem_mb = 1000,
            time_min = 5
        shell:
            "gzip -c {input} > {output}"

    rule chimera_removal:
        input:
            config["pipeline"]["output_dir"] + "/cat.uniques.gz"
        output:
            config["pipeline"]["output_dir"] + "/cat.denoised.nonchimeras"
        threads: 4
        resources:
            mem_mb = 4000,
            time_min = 45
        log:
            config["pipeline"]["output_dir"] + "/chimeraRemoval.log"
        shell:
            "vsearch --uchime3_denovo {input} --sizein --xsize --nonchimeras {output} --relabel 'Zotu' --log {log}"

    rule create_ESV_table:
        input:
            vsearch_global = config["pipeline"]["output_dir"] + "/{sample}.tagged",
            db = config["pipeline"]["output_dir"] + "/cat.denoised.nonchimeras"
        output:
            temp(config["pipeline"]["output_dir"] + "/{sample}.esv.tmp")
        threads: DENOISING_CONFIG.get("threads", 4)
        resources:
            mem_mb = 2000,
            time_min = 30
        log:
            config["pipeline"]["output_dir"] + "/{sample}.table.log"
        shell:
            "vsearch --threads {threads} --search_exact {input.vsearch_global} --db {input.db} --otutabout {output} --log {log}"

    rule merge_sample_ESV_tables:
        input:
            esv_tables = expand(config["pipeline"]["output_dir"] + "/{sample}.esv.tmp", sample=SAMPLES_UNIQUE)
        output:
            config["pipeline"]["output_dir"] + "/ESV.table.tmp"
        threads: 1
        resources:
            mem_mb = 2000,
            time_min = 15
        shell:
            "python3 python_scripts/merge_esv_tables.py {input.esv_tables} > {output}"
