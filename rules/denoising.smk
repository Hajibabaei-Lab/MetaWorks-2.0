# rules/denoising.smk

# Global pooling path
if config["pooling"] == 'Yes':

    rule concatenate_for_global_analysis:
        input:
            expand(config["dir"] + "/{sample}.fasta.tmp", sample=SAMPLES_UNIQUE)
        output:
            temp(config["dir"] + "/cat.fasta.tmp")
        shell:
            "cat {input} > {output}"

    rule edit_fasta_header2:
        input:
            config["dir"] + "/cat.fasta.tmp"
        output:
            temp(config["dir"] + "/cat.fasta")
        shell:
            "sed -e 's/-/_/g' {input} > {output}"

    rule compress:
        input:
            config["dir"] + "/cat.fasta"
        output:
            config["dir"] + "/cat.fasta.gz"
        shell:
            "gzip -c {input} > {output}"

    rule dereplicate:
        input:
            config["dir"] + "/cat.fasta.gz"
        output:
            config["dir"] + "/cat.uniques"
        log:
            config["dir"] + "/dereplication.log"
        shell:
            "vsearch --fastx_uniques {input} --fastaout {output} --sizein --sizeout --log {log}"

    rule denoise:
        input:
            config["dir"] + "/cat.uniques"
        output:
            config["dir"] + "/cat.denoised"
        log:
            config["dir"] + "/denoising.log"
        params:
            minsize = config["VSEARCH_DENOISE"]["minsize"]
        shell:
            "vsearch --cluster_unoise {input} --sizein --sizeout --minsize {params.minsize} --centroids {output} --log {log}"

    rule chimera_removal:
        input:
            config["dir"] + "/cat.denoised"
        output:
            config["dir"] + "/cat.denoised.nonchimeras"
        log:
            config["dir"] + "/chimeraRemoval.log"
        shell:
            "vsearch --uchime3_denovo {input} --sizein --xsize --nonchimeras {output} --relabel 'Zotu' --log {log}"

    rule create_ESV_table:
        input:
            vsearch_global = config["dir"] + "/cat.fasta.gz",
            db = config["dir"] + "/cat.denoised.nonchimeras"
        output:
            temp(config["dir"] + "/ESV.table.tmp")
        threads: config["VSEARCH_TABLE"]["t"]
        log:
            config["dir"] + "/table.log"
        shell:
            "vsearch --threads {threads} --search_exact {input.vsearch_global} --db {input.db} --otutabout {output} --log {log}"


# Per-sample processing path
else:

    rule edit_fasta_header2:
        input:
            config["dir"] + "/{sample}.fasta.tmp"
        output:
            temp(config["dir"] + "/{sample}.tagged")
        shell:
            "sed -e 's/-/_/g' {input} >> {output}"

    rule dereplicate_single:
        input:
            config["dir"] + "/{sample}.tagged"
        output:
            temp(config["dir"] + "/{sample}.uniques.tmp")
        threads:4
        shell:
            "vsearch --fastx_uniques {input} --fastaout {output} --sizein --sizeout"

    rule denoise:
        input:
            config["dir"] + "/{sample}.uniques.tmp"
        output:
            temp(config["dir"] + "/{sample}.denoised")
        log:
            config["dir"] + "/logs/{sample}.denoising.log"
        threads:4
        shell:
            "vsearch --cluster_unoise {input} --sizein --sizeout --minsize {config[VSEARCH_DENOISE][minsize]} --centroids {output} --log {log}"

    rule concatenate_for_global_analysis:
        input:
            expand(config["dir"] + "/{sample}.denoised", sample=SAMPLES_UNIQUE)
        output:
            config["dir"] + "/cat.denoised.tmp"
        threads: 4
        shell:
            "cat {input} > {output}"

    rule dereplicate:
        input:
            config["dir"] + "/cat.denoised.tmp"
        output:
            config["dir"] + "/cat.uniques"
        log:
            config["dir"] + "/dereplication.log"
        threads:4
        shell:
            "vsearch --fastx_uniques {input} --fastaout {output} --sizein --sizeout --log {log}"

    rule compress:
        input:
            config["dir"] + "/cat.uniques"
        output:
            config["dir"] + "/cat.uniques.gz"
        shell:
            "gzip -c {input} > {output}"

    rule chimera_removal:
        input:
            config["dir"] + "/cat.uniques.gz"
        output:
            config["dir"] + "/cat.denoised.nonchimeras"
        log:
            config["dir"] + "/chimeraRemoval.log"
        threads:4
        shell:
            "vsearch --uchime3_denovo {input} --sizein --xsize --nonchimeras {output} --relabel 'Zotu' --log {log}"

    rule create_ESV_table:
        input:
            vsearch_global = config["dir"] + "/{sample}.tagged",
            db = config["dir"] + "/cat.denoised.nonchimeras"
        output:
            temp(config["dir"] + "/{sample}.esv.tmp")
        threads: config["VSEARCH_TABLE"]["t"]
        log:
            config["dir"] + "/{sample}.table.log"
        shell:
            "vsearch --threads {threads} --search_exact {input.vsearch_global} --db {input.db} --otutabout {output} --log {log}"

    rule merge_sample_ESV_tables:
        input:
            esv_tables = expand(config["dir"] + "/{sample}.esv.tmp", sample=SAMPLES_UNIQUE)
        output:
            config["dir"] + "/ESV.table.tmp"
        shell:
            "python3 python_scripts/merge_esv_tables.py {input.esv_tables} > {output}"
