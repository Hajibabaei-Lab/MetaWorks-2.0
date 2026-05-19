def get_global_trial_dirs(config):
    if is_module_enabled(config, "global_otu"):
        return config.get("global_otu", {}).get("trial_dirs", [])
    return config.get("global_esv", {}).get("trial_dirs", [])


def get_global_threads(config):
    if is_module_enabled(config, "global_otu"):
        return config.get("global_otu", {}).get("threads", 4)
    return config.get("global_esv", {}).get("threads", 4)


GLOBAL_OUTPUT_DIR = config["pipeline"]["output_dir"] + "/global"
GLOBAL_TRIAL_DIRS = sorted(get_global_trial_dirs(config))
GLOBAL_THREADS = get_global_threads(config)

rule cat_global_denoised:
    input:
        expand("{trial}/cat.denoised.nonchimeras", trial=GLOBAL_TRIAL_DIRS)
    output:
        GLOBAL_OUTPUT_DIR + "/global_denoised.fasta"
    log:
        GLOBAL_OUTPUT_DIR + "/logs/cat_global_denoised.log"
    shell:
        """
        set -euo pipefail
        cat {input} > {output} 2> {log}
        """

rule dereplicate_global:
    input:
        GLOBAL_OUTPUT_DIR + "/global_denoised.fasta"
    output:
        GLOBAL_OUTPUT_DIR + "/global_derep.fasta"
    threads: GLOBAL_THREADS
    log:
        GLOBAL_OUTPUT_DIR + "/logs/dereplicate_global.log"
    shell:
        """
        set -euo pipefail
        vsearch --threads {threads} --fastx_uniques {input} --fastaout {output} --relabel_sha1 --log {log}
        """
