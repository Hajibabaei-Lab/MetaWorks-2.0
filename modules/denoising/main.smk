# Denoising Module - ESV Generation
# Generates ESVs using VSEARCH with dereplication, denoising, and chimera removal

import os
from pathlib import Path

# Module configuration
MODULE_NAME = "denoising"
MODULE_CONFIG = config.get("modules", {}).get(MODULE_NAME, {})

# Use module config if available, otherwise fall back to main config
VSEARCH_DENOISE_CONFIG = MODULE_CONFIG.get("vsearch_denoise", config.get("VSEARCH_DENOISE", {}))
VSEARCH_TABLE_CONFIG = MODULE_CONFIG.get("vsearch_table", config.get("VSEARCH_TABLE", {}))

# Input validation
def validate_inputs():
    """Validate required configuration parameters"""
    required = ["dir", "SAMPLES_UNIQUE"]
    for req in required:
        if req not in config:
            raise ValueError(f"Missing required config parameter: {req}")
    
    # Validate pooling value
    pooling = config.get("pooling", "Yes")
    if pooling not in ["Yes", "No"]:
        raise ValueError(f"pooling must be 'Yes' or 'No', got: {pooling}")
    
    # Validate minsize parameter
    minsize = VSEARCH_DENOISE_CONFIG.get("minsize", 8)
    if minsize <= 0:
        raise ValueError(f"minsize must be positive, got: {minsize}")

validate_inputs()

# Get samples from config (set by trimming module)
SAMPLES_UNIQUE = config["SAMPLES_UNIQUE"]

# Get pooling option
POOLING = config.get("pooling", "Yes")

# Rule to rename FASTA headers (common to both paths)
rule rename_fasta_headers:
    """Rename FASTA headers to replace hyphens with underscores"""
    input:
        config["dir"] + "/{sample}.fasta.tmp"
    output:
        temp(config["dir"] + "/{sample}.tagged")
    shell:
        "sed -e 's/-/_/g' {input} > {output}"

# Conditional rules based on pooling option
if POOLING == 'Yes':
    # Global pooling path - concatenate all samples first
    
    rule concatenate_samples:
        """Concatenate all samples for global analysis"""
        input:
            expand(config["dir"] + "/{sample}.fasta.tmp", sample=SAMPLES_UNIQUE)
        output:
            temp(config["dir"] + "/cat.fasta.tmp")
        threads: 1
        resources:
            mem_mb = 2000,
            time_min = 10
        shell:
            "cat {input} > {output}"
    
    rule rename_fasta_headers_pooled:
        """Rename FASTA headers for pooled data"""
        input:
            config["dir"] + "/cat.fasta.tmp"
        output:
            temp(config["dir"] + "/cat.fasta")
        shell:
            "sed -e 's/-/_/g' {input} > {output}"
    
    rule compress_pooled_data:
        """Compress pooled FASTA data"""
        input:
            config["dir"] + "/cat.fasta"
        output:
            config["dir"] + "/cat.fasta.gz"
        threads: 1
        resources:
            mem_mb = 1000,
            time_min = 5
        shell:
            "gzip -c {input} > {output}"
    
    rule dereplicate_pooled:
        """Dereplicate pooled sequences to identify unique sequences"""
        input:
            config["dir"] + "/cat.fasta.gz"
        output:
            config["dir"] + "/cat.uniques"
        threads: 1
        resources:
            mem_mb = 4000,
            time_min = 30
        log:
            config["dir"] + "/dereplication.log"
        shell:
            "vsearch --fastx_uniques {input} --fastaout {output} --sizein --sizeout --log {log}"
    
    rule denoise_pooled:
        """Denoise sequences using UNOISE3 algorithm"""
        input:
            config["dir"] + "/cat.uniques"
        output:
            config["dir"] + "/cat.denoised"
        threads: 1
        resources:
            mem_mb = 4000,
            time_min = 60
        log:
            config["dir"] + "/denoising.log"
        params:
            minsize = VSEARCH_DENOISE_CONFIG.get("minsize", 8)
        shell:
            "vsearch --cluster_unoise {input} --sizein --sizeout --minsize {params.minsize} --centroids {output} --log {log}"
    
    rule chimera_removal_pooled:
        """Remove chimeric sequences using de novo chimera detection"""
        input:
            config["dir"] + "/cat.denoised"
        output:
            config["dir"] + "/cat.denoised.nonchimeras"
        threads: 1
        resources:
            mem_mb = 400,
            time_min = 45
        log:
            config["dir"] + "/chimeraRemoval.log"
        shell:
            "vsearch --uchime3_denovo {input} --sizein --xsize --nonchimeras {output} --relabel 'Zotu' --log {log}"
    
    rule create_esv_table_pooled:
        """Create ESV table using exact sequence matching"""
        input:
            vsearch_global = config["dir"] + "/cat.fasta.gz",
            db = config["dir"] + "/cat.denoised.nonchimeras"
        output:
            temp(config["dir"] + "/ESV.table.tmp")
        threads: VSEARCH_TABLE_CONFIG.get("t", 4)
        resources:
            mem_mb = 6000,
            time_min = 60
        log:
            config["dir"] + "/table.log"
        shell:
            "vsearch --threads {threads} --search_exact {input.vsearch_global} --db {input.db} --otutabout {output} --log {log}"

else:
    # Per-sample processing path - process each sample separately then combine
    
    rule dereplicate_per_sample:
        """Dereplicate sequences per sample"""
        input:
            config["dir"] + "/{sample}.tagged"
        output:
            temp(config["dir"] + "/{sample}.uniques.tmp")
        threads: 4
        resources:
            mem_mb = 2000,
            time_min = 20
        shell:
            "vsearch --fastx_uniques {input} --fastaout {output} --sizein --sizeout"
    
    rule denoise_per_sample:
        """Denoise sequences per sample using UNOISE3 algorithm"""
        input:
            config["dir"] + "/{sample}.uniques.tmp"
        output:
            temp(config["dir"] + "/{sample}.denoised")
        threads: 4
        resources:
            mem_mb = 2000,
            time_min = 30
        log:
            config["dir"] + "/logs/{sample}.denoising.log"
        params:
            minsize = VSEARCH_DENOISE_CONFIG.get("minsize", 8)
        shell:
            "vsearch --cluster_unoise {input} --sizein --sizeout --minsize {params.minsize} --centroids {output} --log {log}"
    
    rule concatenate_denoised_samples:
        """Concatenate all denoised samples"""
        input:
            expand(config["dir"] + "/{sample}.denoised", sample=SAMPLES_UNIQUE)
        output:
            config["dir"] + "/cat.denoised.tmp"
        threads: 4
        resources:
            mem_mb = 2000,
            time_min = 10
        shell:
            "cat {input} > {output}"
    
    rule dereplicate_combined:
        """Dereplicate combined denoised sequences"""
        input:
            config["dir"] + "/cat.denoised.tmp"
        output:
            config["dir"] + "/cat.uniques"
        threads: 4
        resources:
            mem_mb = 4000,
            time_min = 30
        log:
            config["dir"] + "/dereplication.log"
        shell:
            "vsearch --fastx_uniques {input} --fastaout {output} --sizein --sizeout --log {log}"
    
    rule compress_uniques:
        """Compress dereplicated uniques"""
        input:
            config["dir"] + "/cat.uniques"
        output:
            config["dir"] + "/cat.uniques.gz"
        threads: 1
        resources:
            mem_mb = 1000,
            time_min = 5
        shell:
            "gzip -c {input} > {output}"
    
    rule chimera_removal_per_sample:
        """Remove chimeric sequences from combined denoised data"""
        input:
            config["dir"] + "/cat.uniques.gz"
        output:
            config["dir"] + "/cat.denoised.nonchimeras"
        threads: 4
        resources:
            mem_mb = 4000,
            time_min = 45
        log:
            config["dir"] + "/chimeraRemoval.log"
        shell:
            "vsearch --uchime3_denovo {input} --sizein --xsize --nonchimeras {output} --relabel 'Zotu' --log {log}"
    
    rule create_esv_table_per_sample:
        """Create ESV table for each sample"""
        input:
            vsearch_global = config["dir"] + "/{sample}.tagged",
            db = config["dir"] + "/cat.denoised.nonchimeras"
        output:
            temp(config["dir"] + "/{sample}.esv.tmp")
        threads: VSEARCH_TABLE_CONFIG.get("t", 4)
        resources:
            mem_mb = 2000,
            time_min = 30
        log:
            config["dir"] + "/{sample}.table.log"
        shell:
            "vsearch --threads {threads} --search_exact {input.vsearch_global} --db {input.db} --otutabout {output} --log {log}"
    
    rule merge_esv_tables:
        """Merge per-sample ESV tables into a single table"""
        input:
            esv_tables = expand(config["dir"] + "/{sample}.esv.tmp", sample=SAMPLES_UNIQUE)
        output:
            config["dir"] + "/ESV.table.tmp"
        threads: 1
        resources:
            mem_mb = 2000,
            time_min = 15
        shell:
            "python3 python_scripts/merge_esv_tables.py {input.esv_tables} > {output}"

# Checkpoint to signal completion of denoising
checkpoint denoising_complete:
    """Signal that denoising and ESV table creation is complete"""
    input:
        config["dir"] + "/ESV.table.tmp",
        config["dir"] + "/cat.denoised.nonchimeras"
    output:
        touch(config["dir"] + "/checkpoints/denoising_complete.done")
    message:
        "Denoising complete for {0} samples (pooling: {1})".format(
            len(SAMPLES_UNIQUE), POOLING
        )
