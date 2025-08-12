# rules/preprocessing.smk

rule preprocess_reads:
    input:
        read1=lambda wildcards: get_fastq_path(wildcards.sample, 1),
        read2=lambda wildcards: get_fastq_path(wildcards.sample, 2)
    output:
        temp(config["dir"] + "/{sample}.processed")
    shell:
        "echo Processing {input.read1} and {input.read2} > {output}"
