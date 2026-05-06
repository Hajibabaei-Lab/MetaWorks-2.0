# rules/itsx/itsx_extraction.smk
# ITSx extraction: extracts ITS1 or ITS2 regions from denoised sequences.
# Active only when modules.itsx_extraction is enabled (ITS_fungi/ITS_plants markers).

ITSX_CONFIG = get_module_config(config, "itsx_extraction")

ITSX_PART = ITSX_CONFIG.get("its_part", "ITS2")

if is_module_enabled(config, "itsx_extraction"):

    rule itsx_extract:
        input:
            config["pipeline"]["output_dir"] + "/cat.denoised.nonchimeras"
        output:
            config["pipeline"]["output_dir"] + f"/ITSx_out.{ITSX_PART}.fasta"
        params:
            prefix = config["pipeline"]["output_dir"] + "/ITSx_out",
        threads: ITSX_CONFIG.get("threads", 4)
        log:
            config["pipeline"]["output_dir"] + "/logs/itsx_extract.log"
        shell:
            """
            set -euo pipefail
            ITSx -i {input} -o {params.prefix} --cpu {threads} 2> {log}
            """

    rule strip_fasta_headers:
        input:
            config["pipeline"]["output_dir"] + f"/ITSx_out.{ITSX_PART}.fasta"
        output:
            config["pipeline"]["output_dir"] + f"/ITSx_out.{ITSX_PART}.fasta.stripped"
        log:
            config["pipeline"]["output_dir"] + "/logs/strip_fasta_headers.log"
        shell:
            """
            set -euo pipefail
            sed 's/^\\(>Zotu[0-9]\\+\\).\\+/\\1/' {input} > {output}
            """
