# rules/pseudogene.smk
# Uses shared helpers from modules/common.smk (get_module_config, is_module_enabled).

PSEUDOGENE_CONFIG = get_module_config(config, "pseudogene_filtering")
CLASSIFICATION_CONFIG = get_module_config(config, "classification")
PSEUDOGENE_SUBMODULE_CONFIG = {
    "dir": config["pipeline"]["output_dir"],
    "pseudogene_filtering": config.get("pseudogene_filtering", {}),
}

module orfs_hmm:
    snakefile: "orfs_hmm.smk"
    config: PSEUDOGENE_SUBMODULE_CONFIG

module orfs_longest:
    snakefile: "orfs_longest.smk"
    config: PSEUDOGENE_SUBMODULE_CONFIG

pseudogene_enabled = (
    is_module_enabled(config, "pseudogene_filtering", default=True)
    and PSEUDOGENE_CONFIG.get("method", "none") in ["hmm", "orf"]
)

if pseudogene_enabled:

    ############################
    # Taxon Subsetting
    ############################

    if PSEUDOGENE_CONFIG.get("grep_type", 1) == 1:
        rule subset_taxonomy_by_taxon1:
            input: config["pipeline"]["output_dir"] + "/rdp.out.tmp"
            output: config["pipeline"]["output_dir"] + "/taxon.zotus"
            params:
                taxon1 = PSEUDOGENE_CONFIG.get("taxon1", "")
            log: config["pipeline"]["output_dir"] + "/logs/pseudogene/taxon_subset.log"
            shell:
                "set -euo pipefail; set +o pipefail; grep {params.taxon1} {input} | awk '{{print $1}}' > \"{output}\" || true"
    else:
        rule subset_taxonomy_by_taxon1_and_taxon2:
            input: config["pipeline"]["output_dir"] + "/rdp.out.tmp"
            output: config["pipeline"]["output_dir"] + "/taxon.zotus"
            params:
                taxon1 = PSEUDOGENE_CONFIG.get("taxon1", ""),
                taxon2 = PSEUDOGENE_CONFIG.get("taxon2", "")
            log: config["pipeline"]["output_dir"] + "/logs/pseudogene/taxon_subset.log"
            shell:
                "set -euo pipefail; set +o pipefail; grep {params.taxon1} {input} | grep {params.taxon2} | awk '{{print $1}}' > \"{output}\" || true"

    rule subset_ESV_sequences_by_taxon:
        input:
            tax = config["pipeline"]["output_dir"] + "/taxon.zotus",
            fas = config["pipeline"]["output_dir"] + "/cat.denoised.nonchimeras"
        output: config["pipeline"]["output_dir"] + "/chimera.denoised.nonchimeras.taxon"
        log: config["pipeline"]["output_dir"] + "/logs/pseudogene/esv_subset.log"
        shell:
            "set -euo pipefail; python3 workflow/scripts/get_taxon_only.py {input.tax} {input.fas} > {output}"

    ############################
    # Strategy 1: Longest ORFs
    ############################

    if PSEUDOGENE_CONFIG.get("method", "hmm") == "orf":
        use rule get_orfs_nt_longest from orfs_longest as get_orfs_nt
        use rule get_longest_orfs from orfs_longest

        rule filter_rdp:
            input:
                orf = config["pipeline"]["output_dir"] + "/longest.orfs.fasta",
                rdp = config["pipeline"]["output_dir"] + "/rdp.out.tmp"
            output: config["pipeline"]["output_dir"] + "/taxonomy.csv"
            params:
                marker = CLASSIFICATION_CONFIG.get("marker", "COI")
            log: config["pipeline"]["output_dir"] + "/logs/pseudogene/filter_rdp.log"
            shell:
                "set -euo pipefail; python3 workflow/scripts/filter_rdp_taxonomy.py {input.orf} {input.rdp} {params.marker} > {output}"

        rule filter_ESV_table:
            input:
                table = config["pipeline"]["output_dir"] + "/ESV.table.tmp",
                orf = config["pipeline"]["output_dir"] + "/longest.orfs.fasta"
            output: config["pipeline"]["output_dir"] + "/ESV.table"
            log: config["pipeline"]["output_dir"] + "/logs/pseudogene/filter_esv_table.log"
            shell:
                "set -euo pipefail; python3 workflow/scripts/filter_ESV_table.py {input.table} {input.orf} > {output}"

        rule add_good_orfs_to_taxonomy:
            input:
                orf = config["pipeline"]["output_dir"] + "/longest.orfs.fasta",
                rdp = config["pipeline"]["output_dir"] + "/rdp.out.tmp"
            output: temp(config["pipeline"]["output_dir"] + "/taxonomy_ORF.tsv")
            log: config["pipeline"]["output_dir"] + "/logs/pseudogene/add_orfs_taxonomy.log"
            shell:
                "set -euo pipefail; python3 workflow/scripts/add_seqs_to_tax4.py {input.orf} {input.rdp} > {output}"

    ############################
    # Strategy 2: HMM-based ORFs
    ############################

    elif PSEUDOGENE_CONFIG.get("method", "hmm") == "hmm":
        use rule get_orfs_nt_hmm from orfs_hmm
        use rule get_orfs_aa_hmm from orfs_hmm
        use rule consolidate_orfs_hmm from orfs_hmm as consolidate_orfs

        rule hmmscan:
            input:
                orf = config["pipeline"]["output_dir"] + "/orfs.fasta.aa.filtered.hmm",
                hmm = PSEUDOGENE_CONFIG.get("hmm_profile", "config/hmm/bold.hmm")
            output: config["pipeline"]["output_dir"] + "/hmm.txt"
            log: config["pipeline"]["output_dir"] + "/logs/pseudogene/hmmscan.log"
            shell:
                "set -euo pipefail; hmmscan --tblout {output} {input.hmm} {input.orf} 2> {log}"

        rule add_good_orf_sequences_to_taxonomy:
            input:
                hmmer = config["pipeline"]["output_dir"] + "/hmm.txt",
                rdp = config["pipeline"]["output_dir"] + "/rdp.out.tmp",
                orfs = config["pipeline"]["output_dir"] + "/orfs.fasta.nt.filtered.hmm"
            output: config["pipeline"]["output_dir"] + "/taxonomy_ORF.tsv"
            log: config["pipeline"]["output_dir"] + "/logs/pseudogene/add_orfs_taxonomy.log"
            shell:
                "set -euo pipefail; python3 workflow/scripts/filter_rdp.py {input.hmmer} {input.orfs} {input.rdp} > {output}"

        rule filter_rdp:
            input:
                orf = config["pipeline"]["output_dir"] + "/orfs.fasta.nt.filtered.hmm",
                rdp = config["pipeline"]["output_dir"] + "/rdp.out.tmp"
            output: config["pipeline"]["output_dir"] + "/taxonomy.csv"
            params:
                marker = CLASSIFICATION_CONFIG.get("marker", "COI")
            log: config["pipeline"]["output_dir"] + "/logs/pseudogene/filter_rdp.log"
            shell:
                "set -euo pipefail; python3 workflow/scripts/filter_rdp_taxonomy.py {input.orf} {input.rdp} {params.marker} > {output}"

        rule filter_ESV_table:
            input:
                table = config["pipeline"]["output_dir"] + "/ESV.table.tmp",
                orf = config["pipeline"]["output_dir"] + "/orfs.fasta.nt.filtered.hmm"
            output: config["pipeline"]["output_dir"] + "/ESV.table"
            log: config["pipeline"]["output_dir"] + "/logs/pseudogene/filter_esv_table.log"
            shell:
                "set -euo pipefail; python3 workflow/scripts/filter_ESV_table.py {input.table} {input.orf} > {output}"

    else:
        raise ValueError(
            f"Invalid pseudogene_filtering.method: {PSEUDOGENE_CONFIG.get('method')!r}. "
            "Expected 'hmm' or 'orf'."
        )
else:
    rule copy_unfiltered_taxonomy:
        input:
            rdp = config["pipeline"]["output_dir"] + "/rdp.out.tmp",
            fas = config["pipeline"]["output_dir"] + "/cat.denoised.nonchimeras"
        output: config["pipeline"]["output_dir"] + "/taxonomy.csv"
        params:
            marker = CLASSIFICATION_CONFIG.get("marker", "COI")
        log: config["pipeline"]["output_dir"] + "/logs/pseudogene/copy_taxonomy.log"
        shell:
            "set -euo pipefail; python3 workflow/scripts/filter_rdp_taxonomy.py {input.fas} {input.rdp} {params.marker} > {output}"

    rule copy_unfiltered_esv_table:
        input: config["pipeline"]["output_dir"] + "/ESV.table.tmp"
        output: config["pipeline"]["output_dir"] + "/ESV.table"
        log: config["pipeline"]["output_dir"] + "/logs/pseudogene/copy_esv_table.log"
        shell:
            "set -euo pipefail; cp {input} {output}"
