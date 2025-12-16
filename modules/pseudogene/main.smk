# Pseudogene Module - Pseudogene Filtering
# Filters putative pseudogenes using ORF analysis (longest ORF or HMM-based methods)

import os
from pathlib import Path

# Module configuration
MODULE_NAME = "pseudogene"
MODULE_CONFIG = config.get("modules", {}).get(MODULE_NAME, {})

# Use module config if available, otherwise fall back to main config
ORFFINDER_CONFIG = MODULE_CONFIG.get("orffinder", config.get("ORFFINDER", {}))
HMM_CONFIG = MODULE_CONFIG.get("hmm", config.get("hmm", {}))
GREP_CONFIG = MODULE_CONFIG.get("grep", {})

# Input validation
def validate_inputs():
    """Validate required configuration parameters"""
    required = ["dir"]
    for req in required:
        if req not in config:
            raise ValueError(f"Missing required config parameter: {req}")
    
    # Validate pseudogene filtering option
    pseudogene_filtering = config.get("pseudogene_filtering", "no")
    if pseudogene_filtering not in ["yes", "no"]:
        raise ValueError(f"pseudogene_filtering must be 'yes' or 'no', got: {pseudogene_filtering}")
    
    # If pseudogene filtering is enabled, validate removal type
    if pseudogene_filtering == "yes":
        removal_type = config.get("removal_type", 1)
        if removal_type not in [1, 2]:
            raise ValueError(f"removal_type must be 1 or 2, got: {removal_type}")
        
        # If using HMM method, validate HMM file exists
        if removal_type == 2:
            hmm_file = HMM_CONFIG.get("path", config.get("hmm", "hmm/bold.hmm"))
            if not os.path.exists(hmm_file):
                raise FileNotFoundError(f"HMM profile file not found: {hmm_file}")

validate_inputs()

# Get pseudogene filtering setting
PSEUDOGENE_FILTERING = config.get("pseudogene_filtering", "no")

# Conditional execution based on pseudogene filtering setting
if PSEUDOGENE_FILTERING == 'yes':
    # Taxon subsetting based on grep type
    if GREP_CONFIG.get("type", config.get("grep_type", 1)) == 1:
        rule subset_taxonomy_by_taxon1:
            """Subset taxonomy by single taxonomic filter"""
            input:
                config["dir"] + "/rdp.out.tmp"
            output:
                config["dir"] + "/taxon.zotus"
            threads: 1
            resources:
                mem_mb = 1000,
                time_min = 5
            shell:
                "grep {config[taxon1]} {input} | awk '{{print $1}}' > {output}"
    else:
        rule subset_taxonomy_by_taxon1_and_taxon2:
            """Subset taxonomy by two taxonomic filters"""
            input:
                config["dir"] + "/rdp.out.tmp"
            output:
                config["dir"] + "/taxon.zotus"
            threads: 1
            resources:
                mem_mb = 1000,
                time_min = 5
            shell:
                "grep {config[taxon1]} {input} | grep {config[taxon2]} | awk '{{print $1}}' > {output}"

    rule subset_ESV_sequences_by_taxon:
        """Subset ESV sequences by taxonomic group"""
        input:
            tax = config["dir"] + "/taxon.zotus",
            fas = config["dir"] + "/cat.denoised.nonchimeras"
        output:
            config["dir"] + "/chimera.denoised.nonchimeras.taxon"
        threads: 1
        resources:
            mem_mb = 2000,
            time_min = 10
        shell:
            "python3 python_scripts/get_taxon_only.py {input.tax} {input.fas} > {output}"

    # Pseudogene removal strategies based on removal type
    if config.get("removal_type", 1) == 1:
        # Strategy 1: Longest ORFs
        rule get_orfs_nt_longest:
            """Get nucleotide ORFs using ORFfinder (longest strategy)"""
            input:
                config["dir"] + "/chimera.denoised.nonchimeras.taxon"
            output:
                config["dir"] + "/orfs.fasta.nt.longest"
            params:
                g = ORFFINDER_CONFIG.get("g", 5),
                s = ORFFINDER_CONFIG.get("s", 2),
                ml = ORFFINDER_CONFIG.get("ml", 30),
                n = ORFFINDER_CONFIG.get("n", 'true'),
                strand = ORFFINDER_CONFIG.get("strand", 'plus')
            threads: 1
            resources:
                mem_mb = 2000,
                time_min = 30
            shell:
                """
                ORFfinder -in {input} \
                    -g {params.g} \
                    -s {params.s} \
                    -ml {params.ml} \
                    -n {params.n} \
                    -strand {params.strand} \
                    -outfmt 1 > {output}
                """

        rule get_longest_orfs:
            """Extract longest ORFs from nucleotide ORFs"""
            input:
                config["dir"] + "/orfs.fasta.nt.longest"
            output:
                config["dir"] + "/longest.orfs.fasta"
            threads: 1
            resources:
                mem_mb = 1000,
                time_min = 10
            shell:
                "python3 python_scripts/parse_orfs3.py {input} {output}"

        rule filter_rdp_longest:
            """Filter RDP results using longest ORFs"""
            input:
                orf = config["dir"] + "/longest.orfs.fasta",
                rdp = config["dir"] + "/rdp.out.tmp"
            output:
                config["dir"] + "/taxonomy.csv"
            threads: 1
            resources:
                mem_mb = 2000,
                time_min = 15
            shell:
                "python3 python_scripts/filter_rdp_taxonomy.py {input.orf} {input.rdp} {config[marker]} > {output}"

        rule filter_ESV_table_longest:
            """Filter ESV table using longest ORFs"""
            input:
                table = config["dir"] + "/ESV.table.tmp",
                orf = config["dir"] + "/longest.orfs.fasta"
            output:
                config["dir"] + "/ESV.table"
            threads: 1
            resources:
                mem_mb = 2000,
                time_min = 10
            shell:
                "python3 python_scripts/filter_ESV_table.py {input.table} {input.orf} > {output}"

        rule add_good_orfs_to_taxonomy_longest:
            """Add good ORFs to taxonomy output"""
            input:
                orf = config["dir"] + "/longest.orfs.fasta",
                rdp = config["dir"] + "/rdp.out.tmp"
            output:
                temp(config["dir"] + "/taxonomy_ORF.tsv")
            threads: 1
            resources:
                mem_mb = 1000,
                time_min = 5
            shell:
                "python3 python_scripts/add_seqs_to_tax4.py {input.orf} {input.rdp} >> {output}"

    elif config.get("removal_type", 1) == 2 and config.get("marker", "COI") == "COI":
        # Strategy 2: HMM-based ORFs (only for COI marker)
        rule get_orfs_nt_hmm:
            """Get nucleotide ORFs using ORFfinder (HMM strategy)"""
            input:
                config["dir"] + "/chimera.denoised.nonchimeras.taxon"
            output:
                config["dir"] + "/orfs.fasta.nt.hmm"
            params:
                g = ORFFINDER_CONFIG.get("g", 5),
                s = ORFFINDER_CONFIG.get("s", 2),
                ml = ORFFINDER_CONFIG.get("ml", 30),
                n = ORFFINDER_CONFIG.get("n", 'true'),
                strand = ORFFINDER_CONFIG.get("strand", 'plus')
            threads: 1
            resources:
                mem_mb = 2000,
                time_min = 30
            shell:
                """
                ORFfinder -in {input} \
                    -g {params.g} \
                    -s {params.s} \
                    -ml {params.ml} \
                    -n {params.n} \
                    -strand {params.strand} \
                    -outfmt 1 > {output}
                """

        rule get_orfs_aa_hmm:
            """Get amino acid ORFs using ORFfinder (HMM strategy)"""
            input:
                config["dir"] + "/chimera.denoised.nonchimeras.taxon"
            output:
                config["dir"] + "/orfs.fasta.aa.hmm"
            params:
                g = ORFFINDER_CONFIG.get("g", 5),
                s = ORFFINDER_CONFIG.get("s", 2),
                ml = ORFFINDER_CONFIG.get("ml", 30),
                n = ORFFINDER_CONFIG.get("n", 'true'),
                strand = ORFFINDER_CONFIG.get("strand", 'plus')
            threads: 1
            resources:
                mem_mb = 2000,
                time_min = 30
            shell:
                """
                ORFfinder -in {input} \
                    -g {params.g} \
                    -s {params.s} \
                    -ml {params.ml} \
                    -n {params.n} \
                    -strand {params.strand} \
                    -outfmt 0 > {output}
                """

        rule consolidate_orfs_hmm:
            """Consolidate nucleotide and amino acid ORFs (HMM strategy)"""
            input:
                nt = config["dir"] + "/orfs.fasta.nt.hmm",
                aa = config["dir"] + "/orfs.fasta.aa.hmm"
            output:
                nt2 = config["dir"] + "/orfs.fasta.nt.filtered.hmm",
                aa2 = config["dir"] + "/orfs.fasta.aa.filtered.hmm"
            threads: 1
            resources:
                mem_mb = 1000,
                time_min = 10
            shell:
                "python3 python_scripts/parse_orfs4.py {input.nt} {input.aa} {output.nt2} {output.aa2}"

        rule hmmscan:
            """Scan ORFs against HMM profile"""
            input:
                orf = config["dir"] + "/orfs.fasta.aa.filtered.hmm",
                hmm = HMM_CONFIG.get("path", config.get("hmm", "hmm/bold.hmm"))
            output:
                config["dir"] + "/hmm.txt"
            threads: 1
            resources:
                mem_mb = 2000,
                time_min = 45
            shell:
                "hmmscan --tblout {output} {input.hmm} {input.orf}"

        rule add_good_orf_sequences_to_taxonomy_hmm:
            """Add good ORF sequences to taxonomy output (HMM strategy)"""
            input:
                hmmer = config["dir"] + "/hmm.txt",
                rdp = config["dir"] + "/rdp.out.tmp",
                orfs = config["dir"] + "/orfs.fasta.nt.filtered.hmm"
            output:
                config["dir"] + "/taxonomy_ORF.tsv"
            threads: 1
            resources:
                mem_mb = 1000,
                time_min = 10
            shell:
                "python3 python_scripts/filter_rdp.py {input.hmmer} {input.orfs} {input.rdp} {config[marker]} >> {output}"

        rule filter_rdp_hmm:
            """Filter RDP results using HMM-based ORFs"""
            input:
                orf = config["dir"] + "/orfs.fasta.nt.filtered.hmm",
                rdp = config["dir"] + "/rdp.out.tmp"
            output:
                config["dir"] + "/taxonomy.csv"
            threads: 1
            resources:
                mem_mb = 2000,
                time_min = 15
            shell:
                "python3 python_scripts/filter_rdp_taxonomy.py {input.orf} {input.rdp} {config[marker]} > {output}"

        rule filter_ESV_table_hmm:
            """Filter ESV table using HMM-based ORFs"""
            input:
                table = config["dir"] + "/ESV.table.tmp",
                orf = config["dir"] + "/orfs.fasta.nt.filtered.hmm"
            output:
                config["dir"] + "/ESV.table"
            threads: 1
            resources:
                mem_mb = 2000,
                time_min = 10
            shell:
                "python3 python_scripts/filter_ESV_table.py {input.table} {input.orf} > {output}"

    else:
        print("ERROR: Invalid or missing 'removal_type' or unsupported 'marker' in config.")

# If pseudogene filtering is not enabled, create the final outputs directly
else:
    # Copy unfiltered taxonomy and ESV table
    rule copy_unfiltered_taxonomy:
        """Copy unfiltered taxonomy when pseudogene filtering is disabled"""
        input:
            config["dir"] + "/rdp.out.tmp"
        output:
            config["dir"] + "/taxonomy.csv"
        threads: 1
        resources:
            mem_mb = 1000,
            time_min = 5
        shell:
            "python3 python_scripts/filter_rdp_taxonomy.py {input} {input} {config[marker]} > {output}"
    
    rule copy_unfiltered_esv_table:
        """Copy unfiltered ESV table when pseudogene filtering is disabled"""
        input:
            config["dir"] + "/ESV.table.tmp"
        output:
            config["dir"] + "/ESV.table"
        threads: 1
        resources:
            mem_mb = 1000,
            time_min = 5
        shell:
            "cp {input} {output}"

# Checkpoint to signal completion of pseudogene filtering
checkpoint pseudogene_filtering_complete:
    """Signal that pseudogene filtering is complete"""
    input:
        config["dir"] + "/taxonomy.csv",
        config["dir"] + "/ESV.table"
    output:
        touch(config["dir"] + "/checkpoints/pseudogene_complete.done")
    message:
        "Pseudogene filtering complete (filtering: {0}, method: {1})".format(
            PSEUDOGENE_FILTERING,
            config.get("removal_type", "none")
        )
