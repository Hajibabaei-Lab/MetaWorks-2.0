# rules/adapter_trimming.smk

rule pair_reads:
    input:
        f = lambda wildcards: get_fastq_path(wildcards.sample, 1),
        r = lambda wildcards: get_fastq_path(wildcards.sample, 2)
    output:
        X1 = temp("{sample}_R1.out"),
        X2 = temp("{sample}_R2.out"),
        s = config["dir"] + "/paired/{sample}.fastq.gz"
    params:
        q = lambda wc: config["SEQPREP"]["q"],
        m = lambda wc: config["SEQPREP"]["m"],
        n = lambda wc: config["SEQPREP"]["n"],
        o = lambda wc: config["SEQPREP"]["o"]
    shell:
        """
        SeqPrep \
            -f {input.f} \
            -r {input.r} \
            -1 {output.X1} \
            -2 {output.X2} \
            -q {params.q} \
            -m {params.m} \
            -n {params.n} \
            -s {output.s} \
            -o {params.o}
        """

rule trim_linked_adapters:
    input:
        adapters = config["CUTADAPT"]["fasta"],
        paired = config["dir"] + "/paired/{sample}.fastq.gz"
    output:
        config["dir"] + "/trimmed/{sample}.fasta"
    params:
        m = lambda wc: config["CUTADAPT"]["m"],
        q = lambda wc: config["CUTADAPT"]["q"],
        e = lambda wc: config["CUTADAPT"]["e"],
        O = lambda wc: config["CUTADAPT"]["O"],
        mn = lambda wc: config["CUTADAPT"]["mn"]
    shell:
        """
        cutadapt \
            -a file:{input.adapters} \
            -m {params.m} \
            -q {params.q} \
            -e {params.e} \
            -O {params.O} \
            --max-n={params.mn} \
            --prefix {{name}} \
            --discard-untrimmed \
            --output {output} \
            {input.paired} \
        """

rule gzip_trimmed_fasta:
    input:
        config["dir"] + "/trimmed/{sample}.fasta"
    output:
        config["dir"] + "/trimmed/{sample}.fasta.gz"
    shell:
        "gzip -c {input} > {output}"