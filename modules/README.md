# MetaWorks Modules

This directory contains modular components of the MetaWorks pipeline. Each module is self-contained and follows standardized conventions for inputs, outputs, and configuration.

## Module Structure

Each module should have:
- One or more `.smk` (Snakemake) files containing rules
- A `module.yaml` metadata file describing the module
- Optional `envs/` directory for conda environments
- Documentation of inputs, outputs, and parameters

## Available Modules

### Core Pipeline Modules
- `preprocessing/` - Read pairing and quality filtering
- `trimming/` - Adapter trimming with Cutadapt
- `denoising/` - ESV generation with VSEARCH
- `classification/` - Taxonomic assignment with RDP
- `pseudogene/` - Pseudogene filtering (ORF detection, HMM)
- `postprocessing/` - ESV table generation and merging
- `reporting/` - Results formatting and statistics

### Utility Modules
- `common/` - Shared utilities (I/O, validation, checkpoints)

## Module Standards

See `docs/MODULE_STANDARDS.md` for detailed module development guidelines.
