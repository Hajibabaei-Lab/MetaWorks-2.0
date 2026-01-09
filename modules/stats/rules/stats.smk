# rules/stats.smk

# -------------------------------------
# Raw FASTQ R1
# -------------------------------------
rule raw1_stats:
    input: lambda wildcards: get_fastq_path(wildcards.sample, 1)
    output: temp(config["pipeline"]["output_dir"] + "/stats/{sample}.R1stats")
    shell: "python3 python_scripts/fastq_gz_stats.py {input} > {output}"

rule cat_raw1_stats:
    input: expand(config["pipeline"]["output_dir"] + "/stats/{sample}.R1stats", sample=SAMPLES_UNIQUE)
    output: temp(config["pipeline"]["output_dir"] + "/stats/R1.stats.tmp")
    shell: "cat {input} > {output}"

rule cat_raw1_stats_add_headers:
    input: config["pipeline"]["output_dir"] + "/stats/R1.stats.tmp"
    output: config["pipeline"]["output_dir"] + "/stats/R1.stats"
    shell: "printf 'FileName\tTotSeqs\tMinLength\tMaxLength\tMeanLength\tMedianLength\tModeLength\n' > {output} && cat {input} >> {output}"


# -------------------------------------
# Raw FASTQ R2
# -------------------------------------
rule raw2_stats:
    input: lambda wildcards: get_fastq_path(wildcards.sample, 2)
    output: temp(config["pipeline"]["output_dir"] + "/stats/{sample}.R2stats")
    shell: "python3 python_scripts/fastq_gz_stats.py {input} > {output}"

rule cat_raw2_stats:
    input: expand(config["pipeline"]["output_dir"] + "/stats/{sample}.R2stats", sample=SAMPLES_UNIQUE)
    output: temp(config["pipeline"]["output_dir"] + "/stats/R2.stats.tmp")
    shell: "cat {input} > {output}"

rule cat_raw2_stats_add_headers:
    input: config["pipeline"]["output_dir"] + "/stats/R2.stats.tmp"
    output: config["pipeline"]["output_dir"] + "/stats/R2.stats"
    shell: "printf 'FileName\tTotSeqs\tMinLength\tMaxLength\tMeanLength\tMedianLength\tModeLength\n' > {output} && cat {input} >> {output}"


# -------------------------------------
# Paired FASTQ
# -------------------------------------
rule paired_stats:
    input: config["pipeline"]["output_dir"] + "/paired/{sample}.fastq.gz"
    output: temp(config["pipeline"]["output_dir"] + "/stats/{sample}.Pstats")
    shell: "python3 python_scripts/fastq_gz_stats.py {input} > {output}"

rule cat_paired_stats:
    input: expand(config["pipeline"]["output_dir"] + "/stats/{sample}.Pstats", sample=SAMPLES_UNIQUE)
    output: temp(config["pipeline"]["output_dir"] + "/stats/paired.stats.tmp")
    shell: "cat {input} > {output}"

rule cat_paired_stats_add_headers:
    input: config["pipeline"]["output_dir"] + "/stats/paired.stats.tmp"
    output: config["pipeline"]["output_dir"] + "/stats/paired.stats"
    shell: "printf 'FileName\tTotSeqs\tMinLength\tMaxLength\tMeanLength\tMedianLength\tModeLength\n' > {output} && cat {input} >> {output}"


# -------------------------------------
# Trimmed FASTA
# -------------------------------------


rule trimmed_stats:
    input: config["pipeline"]["output_dir"] + "/trimmed/{sample}.fasta.gz"
    output: temp(config["pipeline"]["output_dir"] + "/stats/{sample}.trimmedstats")
    shell: "python3 python_scripts/fasta_gz_stats.py {input} >> {output}"

rule cat_trimmed_stats:
    input: expand(config["pipeline"]["output_dir"] + "/stats/{sample}.trimmedstats", sample=SAMPLES_UNIQUE)
    output: temp(config["pipeline"]["output_dir"] + "/stats/trimmed.stats.tmp")
    shell: "cat {input} > {output}"

rule cat_trimmed_stats_add_headers:
    input: config["pipeline"]["output_dir"] + "/stats/trimmed.stats.tmp"
    output: config["pipeline"]["output_dir"] + "/stats/trimmed.stats"
    shell: "printf 'FileName\tTotSeqs\tMinLength\tMaxLength\tMeanLength\tMedianLength\tModeLength\n' > {output} && cat {input} >> {output}"

rule edit_fasta_header1:
    input: config["pipeline"]["output_dir"] + "/trimmed/{sample}.fasta.gz"
    output: temp(config["pipeline"]["output_dir"] + "/{sample}.fasta.tmp")
    shell: "python3 python_scripts/rename_fasta_gzip.py {input} > {output}"
