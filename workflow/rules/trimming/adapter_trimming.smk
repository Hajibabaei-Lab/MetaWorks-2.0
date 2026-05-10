TRIMMING_CONFIG = get_module_config(config, "trimming")
PREPROCESSING_CONFIG = get_module_config(config, "preprocessing")
READ_MODE = TRIMMING_CONFIG.get("read_mode", "paired")
ADAPTER_SOURCE = TRIMMING_CONFIG.get("adapter_source", "file")

if READ_MODE == "paired" and ADAPTER_SOURCE == "file":

    rule pair_reads:
        input:
            f = lambda wildcards: get_fastq_path(wildcards.sample, 1),
            r = lambda wildcards: get_fastq_path(wildcards.sample, 2)
        output:
            X1 = temp(config["pipeline"]["output_dir"] + "/{sample}_R1.out"),
            X2 = temp(config["pipeline"]["output_dir"] + "/{sample}_R2.out"),
            s = config["pipeline"]["output_dir"] + "/paired/{sample}.fastq.gz"
        params:
            q = lambda wc: PREPROCESSING_CONFIG.get("quality_score", 13),
            o = lambda wc: PREPROCESSING_CONFIG.get("min_overlap", 25),
            m = lambda wc: PREPROCESSING_CONFIG.get("max_mismatch", 0.02),
            n = lambda wc: PREPROCESSING_CONFIG.get("min_match", 0.90)
        threads: 1
        log:
            config["pipeline"]["output_dir"] + "/logs/pairing/{sample}.log"
        benchmark:
            config["pipeline"]["output_dir"] + "/benchmarks/pairing/{sample}.txt"
        shell:
            """
            set -euo pipefail
            SeqPrep \
                -f {input.f} \
                -r {input.r} \
                -1 {output.X1} \
                -2 {output.X2} \
                -q {params.q} \
                -o {params.o} \
                -m {params.m} \
                -n {params.n} \
                -s {output.s} \
                2> >(tee "{log}" >&2)
            """

    rule trim_linked_adapters:
        input:
            adapters = TRIMMING_CONFIG.get("adapters"),
            paired = config["pipeline"]["output_dir"] + "/paired/{sample}.fastq.gz"
        output:
            config["pipeline"]["output_dir"] + "/trimmed/{sample}.fasta"
        params:
            m = lambda wc: TRIMMING_CONFIG.get("min_length", 150),
            q = lambda wc: TRIMMING_CONFIG.get("quality_cutoff", "20,20"),
            e = lambda wc: TRIMMING_CONFIG.get("error_rate", 0.1),
            O = lambda wc: TRIMMING_CONFIG.get("min_adapter_overlap", 3),
            mn = lambda wc: TRIMMING_CONFIG.get("max_n_bases", 3),
            rc = lambda wc: "--revcomp" if TRIMMING_CONFIG.get("enable_rc", False) else ""
        threads: 1
        log:
            config["pipeline"]["output_dir"] + "/logs/trimming/{sample}.log"
        benchmark:
            config["pipeline"]["output_dir"] + "/benchmarks/trimming/{sample}.txt"
        shell:
            """
            set -euo pipefail
            cutadapt \
                -a file:{input.adapters} \
                -m {params.m} \
                -q {params.q} \
                -e {params.e} \
                -O {params.O} \
                --max-n={params.mn} \
                {params.rc} \
                --prefix {{name}} \
                --discard-untrimmed \
                --fasta \
                --output {output} \
                {input.paired} \
                2> >(tee "{log}" >&2)
            """

elif READ_MODE == "paired" and ADAPTER_SOURCE == "csv":

    rule pair_reads:
        input:
            f = lambda wildcards: get_fastq_path(wildcards.sample, 1),
            r = lambda wildcards: get_fastq_path(wildcards.sample, 2)
        output:
            X1 = temp(config["pipeline"]["output_dir"] + "/{sample}_R1.out"),
            X2 = temp(config["pipeline"]["output_dir"] + "/{sample}_R2.out"),
            s = config["pipeline"]["output_dir"] + "/paired/{sample}.fastq.gz"
        params:
            q = lambda wc: PREPROCESSING_CONFIG.get("quality_score", 13),
            o = lambda wc: PREPROCESSING_CONFIG.get("min_overlap", 25),
            m = lambda wc: PREPROCESSING_CONFIG.get("max_mismatch", 0.02),
            n = lambda wc: PREPROCESSING_CONFIG.get("min_match", 0.90)
        threads: 1
        log:
            config["pipeline"]["output_dir"] + "/logs/pairing/{sample}.log"
        benchmark:
            config["pipeline"]["output_dir"] + "/benchmarks/pairing/{sample}.txt"
        shell:
            """
            set -euo pipefail
            SeqPrep \
                -f {input.f} \
                -r {input.r} \
                -1 {output.X1} \
                -2 {output.X2} \
                -q {params.q} \
                -o {params.o} \
                -m {params.m} \
                -n {params.n} \
                -s {output.s} \
                2> >(tee "{log}" >&2)
            """

    rule format_base_adapters:
        input:
            TRIMMING_CONFIG.get("adapter_csv", "")
        output:
            config["pipeline"]["output_dir"] + "/adapters.fasta"
        log:
            config["pipeline"]["output_dir"] + "/logs/format_adapters.log"
        shell:
            """
            set -euo pipefail
            python3 workflow/scripts/formatAdapters_anchored.py {input} > {output} 2> {log}
            """

    rule generate_per_sample_adapters:
        input:
            config["pipeline"]["output_dir"] + "/adapters.fasta"
        output:
            expand(config["pipeline"]["output_dir"] + "/{sample}_adapters.fasta", sample=SAMPLES_UNIQUE)
        log:
            config["pipeline"]["output_dir"] + "/logs/generate_per_sample_adapters.log"
        params:
            samples_file = config["pipeline"]["output_dir"] + "/unique_samples.txt",
            output_dir = config["pipeline"]["output_dir"],
            samples = " ".join(SAMPLES_UNIQUE)
        shell:
            """
            set -euo pipefail
            printf '%%s\\n' {params.samples} > {params.samples_file}
            python3 workflow/scripts/formatAdapters_anchored_filename.py {input} {params.samples_file} --output-dir {params.output_dir} 2> {log}
            """

    rule trim_per_sample_adapters:
        input:
            adapters = config["pipeline"]["output_dir"] + "/{sample}_adapters.fasta",
            paired = config["pipeline"]["output_dir"] + "/paired/{sample}.fastq.gz"
        output:
            config["pipeline"]["output_dir"] + "/trimmed/{sample}.fasta"
        params:
            m = lambda wc: TRIMMING_CONFIG.get("min_length", 150),
            q = lambda wc: TRIMMING_CONFIG.get("quality_cutoff", "20,20"),
            e = lambda wc: TRIMMING_CONFIG.get("error_rate", 0.1),
            O = lambda wc: TRIMMING_CONFIG.get("min_adapter_overlap", 3),
            mn = lambda wc: TRIMMING_CONFIG.get("max_n_bases", 3),
            rc = lambda wc: "--revcomp" if TRIMMING_CONFIG.get("enable_rc", False) else ""
        threads: 1
        log:
            config["pipeline"]["output_dir"] + "/logs/trimming/{sample}.log"
        benchmark:
            config["pipeline"]["output_dir"] + "/benchmarks/trimming/{sample}.txt"
        shell:
            """
            set -euo pipefail
            cutadapt \
                -g file:{input.adapters} \
                -m {params.m} \
                -q {params.q} \
                -e {params.e} \
                -O {params.O} \
                --max-n={params.mn} \
                {params.rc} \
                --prefix {{name}} \
                --discard-untrimmed \
                --fasta \
                --output {output} \
                {input.paired} \
                2> >(tee "{log}" >&2)
            """

elif READ_MODE == "single":

    PROCESS_AS = TRIMMING_CONFIG.get("process_as", "F")
    PRIMER = TRIMMING_CONFIG.get("primer", "")

    if PROCESS_AS == "R":

        rule reverse_complement:
            input:
                lambda wildcards: get_fastq_path(wildcards.sample, 2)
            output:
                config["pipeline"]["output_dir"] + "/rc/{sample}.fastq.gz"
            threads: 1
            log:
                config["pipeline"]["output_dir"] + "/logs/rc/{sample}.log"
            shell:
                """
                set -euo pipefail
                python3 workflow/scripts/rc.py {input} {output} 2> {log}
                """

        rule trim_single_primer:
            input:
                config["pipeline"]["output_dir"] + "/rc/{sample}.fastq.gz"
            output:
                config["pipeline"]["output_dir"] + "/trimmed/{sample}.fasta"
            params:
                m = lambda wc: TRIMMING_CONFIG.get("min_length", 150),
                q = lambda wc: TRIMMING_CONFIG.get("quality_cutoff", "20,20"),
                e = lambda wc: TRIMMING_CONFIG.get("error_rate", 0.1),
                mn = lambda wc: TRIMMING_CONFIG.get("max_n_bases", 3),
                primer = PRIMER,
                flag = "-a"
            threads: 1
            log:
                config["pipeline"]["output_dir"] + "/logs/trimming/{sample}.log"
            benchmark:
                config["pipeline"]["output_dir"] + "/benchmarks/trimming/{sample}.txt"
            shell:
                """
                set -euo pipefail
                cutadapt \
                    {params.flag} {params.primer} \
                    -m {params.m} \
                    -q {params.q} \
                    -e {params.e} \
                    --max-n={params.mn} \
                    --discard-untrimmed \
                    --fasta \
                    --output {output} \
                    {input} \
                    2> >(tee "{log}" >&2)
                """

    else:

        rule trim_single_primer:
            input:
                lambda wildcards: get_fastq_path(wildcards.sample, 1)
            output:
                config["pipeline"]["output_dir"] + "/trimmed/{sample}.fasta"
            params:
                m = lambda wc: TRIMMING_CONFIG.get("min_length", 150),
                q = lambda wc: TRIMMING_CONFIG.get("quality_cutoff", "20,20"),
                e = lambda wc: TRIMMING_CONFIG.get("error_rate", 0.1),
                mn = lambda wc: TRIMMING_CONFIG.get("max_n_bases", 3),
                primer = PRIMER,
                flag = "-g"
            threads: 1
            log:
                config["pipeline"]["output_dir"] + "/logs/trimming/{sample}.log"
            benchmark:
                config["pipeline"]["output_dir"] + "/benchmarks/trimming/{sample}.txt"
            shell:
                """
                set -euo pipefail
                cutadapt \
                    {params.flag} {params.primer} \
                    -m {params.m} \
                    -q {params.q} \
                    -e {params.e} \
                    --max-n={params.mn} \
                    --discard-untrimmed \
                    --fasta \
                    --output {output} \
                    {input} \
                    2> >(tee "{log}" >&2)
                """

rule gzip_trimmed_fasta:
    input:
        config["pipeline"]["output_dir"] + "/trimmed/{sample}.fasta"
    output:
        config["pipeline"]["output_dir"] + "/trimmed/{sample}.fasta.gz"
    threads: 1
    log:
        config["pipeline"]["output_dir"] + "/logs/trimming/{sample}.gzip.log"
    shell:
        "gzip -c {input} > {output} 2> {log}"
