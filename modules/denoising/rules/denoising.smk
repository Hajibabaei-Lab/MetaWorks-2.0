# rules/denoising.smk

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

DENOISING_CONFIG = get_module_config(config, "denoising")

# Global pooling path
if DENOISING_CONFIG.get("pool_samples", True):

    rule concatenate_for_global_analysis:
        input:
            expand(config["pipeline"]["output_dir"] + "/{sample}.fasta.tmp", sample=SAMPLES_UNIQUE)
        output:
            temp(config["pipeline"]["output_dir"] + "/cat.fasta.tmp")
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
        shell:
            "gzip -c {input} > {output}"

    rule dereplicate:
        input:
            config["pipeline"]["output_dir"] + "/cat.fasta.gz"
        output:
            config["pipeline"]["output_dir"] + "/cat.uniques"
        log:
            config["pipeline"]["output_dir"] + "/dereplication.log"
        shell:
            "vsearch --fastx_uniques {input} --fastaout {output} --sizein --sizeout --log {log}"

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
            "vsearch --cluster_unoise {input} --sizein --sizeout --minsize {params.minsize} --centroids {output} --log {log}"

    rule chimera_removal:
        input:
            config["pipeline"]["output_dir"] + "/cat.denoised"
        output:
            config["pipeline"]["output_dir"] + "/cat.denoised.nonchimeras"
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
        threads:4
        shell:
            "vsearch --fastx_uniques {input} --fastaout {output} --sizein --sizeout"

    rule denoise:
        input:
            config["pipeline"]["output_dir"] + "/{sample}.uniques.tmp"
        output:
            temp(config["pipeline"]["output_dir"] + "/{sample}.denoised")
        log:
            config["pipeline"]["output_dir"] + "/logs/{sample}.denoising.log"
        threads:4
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
        shell:
            "cat {input} > {output}"

    rule dereplicate:
        input:
            config["pipeline"]["output_dir"] + "/cat.denoised.tmp"
        output:
            config["pipeline"]["output_dir"] + "/cat.uniques"
        log:
            config["pipeline"]["output_dir"] + "/dereplication.log"
        threads:4
        shell:
            "vsearch --fastx_uniques {input} --fastaout {output} --sizein --sizeout --log {log}"

    rule compress:
        input:
            config["pipeline"]["output_dir"] + "/cat.uniques"
        output:
            config["pipeline"]["output_dir"] + "/cat.uniques.gz"
        shell:
            "gzip -c {input} > {output}"

    rule chimera_removal:
        input:
            config["pipeline"]["output_dir"] + "/cat.uniques.gz"
        output:
            config["pipeline"]["output_dir"] + "/cat.denoised.nonchimeras"
        log:
            config["pipeline"]["output_dir"] + "/chimeraRemoval.log"
        threads:4
        shell:
            "vsearch --uchime3_denovo {input} --sizein --xsize --nonchimeras {output} --relabel 'Zotu' --log {log}"

    rule create_ESV_table:
        input:
            vsearch_global = config["pipeline"]["output_dir"] + "/{sample}.tagged",
            db = config["pipeline"]["output_dir"] + "/cat.denoised.nonchimeras"
        output:
            temp(config["pipeline"]["output_dir"] + "/{sample}.esv.tmp")
        threads: DENOISING_CONFIG.get("threads", 4)
        log:
            config["pipeline"]["output_dir"] + "/{sample}.table.log"
        shell:
            "vsearch --threads {threads} --search_exact {input.vsearch_global} --db {input.db} --otutabout {output} --log {log}"

    rule merge_sample_ESV_tables:
        input:
            esv_tables = expand(config["pipeline"]["output_dir"] + "/{sample}.esv.tmp", sample=SAMPLES_UNIQUE)
        output:
            config["pipeline"]["output_dir"] + "/ESV.table.tmp"
        shell:
            "python3 python_scripts/merge_esv_tables.py {input.esv_tables} > {output}"
