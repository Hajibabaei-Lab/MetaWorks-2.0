# MetaWorks Module Standards

## Overview

This document describes the actual module architecture of the MetaWorks 2.0 pipeline. Modules are Snakemake rule files (`.smk`) organized under `workflow/rules/`, backed by Python helper scripts in `workflow/scripts/`, and governed by a central registry in `lib/config/module_registry.py`.

## Pipeline Rules

### Directory Structure

Rule files live under `workflow/rules/` in subdirectories organized by pipeline stage:

```
workflow/rules/
├── common.smk                          # Shared helpers (included first by Snakefile)
├── trimming/
│   └── adapter_trimming.smk
├── denoising/
│   └── denoising.smk
├── clustering/
│   └── clustering.smk
├── itsx/
│   └── itsx_extraction.smk
├── classification/
│   ├── classifier.smk
│   ├── rdp_classifier.smk
│   └── sintax_classifier.smk
├── pseudogene/
│   ├── pseudogene.smk
│   ├── orfs_hmm.smk
│   └── orfs_longest.smk
├── stats/
│   └── stats.smk
├── utils/
│   ├── utils.smk
│   └── results.smk
└── global/
    ├── global_common.smk
    ├── global_esv.smk
    └── global_otu.smk
```

There are 17 `.smk` files across these subdirectories. The root `workflow/Snakefile` loads `common.smk` first, then dynamically includes module snakefiles based on the module registry.

### Dynamic Loading

The Snakefile (`workflow/Snakefile`) uses the module registry to determine which `.smk` files to include:

```python
from lib.config.module_registry import (
    resolve_module_include_paths,
    resolve_terminal_targets,
    validate_module_selection,
)

include: "rules/common.smk"

dependency_errors = validate_module_selection(config)

for module_include_path in resolve_module_include_paths(config, repo_root=str(REPO_ROOT)):
    include: module_include_path
```

Modules are loaded in topological dependency order with cycle detection.

## Helper Scripts

### `workflow/scripts/`

There are 24 Python scripts handling data processing tasks. All scripts use `argparse` for CLI argument parsing and are invoked from shell directives in `.smk` rules.

Key scripts by function:

| Script | Purpose |
|--------|---------|
| `marker_defs.py` | Single source of truth for marker definitions |
| `formatAdapters_anchored.py` | Format adapter sequences for Cutadapt |
| `formatAdapters_anchored_filename.py` | Format per-sample adapters from CSV |
| `filter_rdp.py` | Filter RDP classification results |
| `filter_rdp_taxonomy.py` | Filter taxonomy by confidence |
| `parallel_rdp.py` | Parallel RDP classifier — splits FASTA into chunks, dispatches concurrent JVMs, concatenates results |
| `add_abundance_to_rdp_out.py` | Merge abundance data with RDP output |
| `add_seqs_to_tax3.py` | Add sequences to taxonomy output |
| `add_seqs_to_tax4.py` | Add sequences to taxonomy output (alt format) |
| `parse_orfs3.py` | Parse ORF predictions (v3 format) |
| `parse_orfs4.py` | Parse ORF predictions (v4 format) |
| `sintax_to_rdp_out.py` | Convert SINTAX output to RDP format |
| `rdp_tsv_to_csv.py` | Convert RDP TSV to CSV |
| `get_taxon_only.py` | Extract taxon-only column |
| `grab_seqs_from_results.py` | Extract sequences from results file |
| `filter_ESV_table.py` | Filter ESV abundance table |
| `merge_esv_tables.py` | Merge multiple ESV tables |
| `map_global_esvs_to_results.py` | Map global ESVs to per-trial results |
| `map_global_otus_to_results.py` | Map global OTUs to per-trial results |
| `map_global_to_results.py` | Map global identifiers to results |
| `prepare_pooled_reads.py` | Prepare reads for pooled denoising (sample-prefix rename + header sanitize + gzip) |
| `fastq_gz_stats.py` | Compute statistics for FASTQ files |
| `fasta_gz_stats.py` | Compute statistics for FASTA files |
| `rc.py` | Reverse-complement sequences |

### `marker_defs.py` — Single Source of Truth

`workflow/scripts/marker_defs.py` defines marker properties (taxonomy regions, primer sequences, expected lengths) for the 13 supported markers. All pipeline components that need marker information consult this file.

### `parallel_rdp.py` — Parallel RDP Classifier Execution

`workflow/scripts/parallel_rdp.py` is the production entry point for RDP taxonomic assignment. It is invoked by the `taxonomic_assignment` rule in `workflow/rules/classification/rdp_classifier.smk` and provides multi-core parallelism for the Java-based RDP classifier.

**Execution model:**

1. Reads the entire input FASTA using Biopython and splits records evenly into *N* chunk files (where *N* = `--threads`).
2. Dispatches each chunk to a separate `rdp_classifier` JVM subprocess via `concurrent.futures.ThreadPoolExecutor` with a 12-hour per-chunk timeout.
3. Each JVM receives `rdp_classifier -Xmx{memory}G classify {options} -o {result} {chunk}`.
4. Resolves classifier properties paths from multiple locations (absolute, CWD-relative, `runtime/classifiers/`, `config/classifiers/`).
5. Concatenates all chunk results into a single output file and cleans up temporary files.

**CLI interface:**

```
python3 workflow/scripts/parallel_rdp.py \
    --input <fasta> \
    --output <result_file> \
    --threads <N> \
    --memory <Xg> \
    --options '<rdp_options>'
```

**Resource note:** Total JVM memory consumption is `threads × memory`. See the Configuration Guide for tuning guidance.

## Module Registry

### `lib/config/module_registry.py`

The `MODULE_REGISTRY` dictionary maps 11 module names to their configuration entries. Each entry specifies:

- `module` — name, version, description, author, `enabled_by_default`
- `directory` — subdirectory under `workflow/rules/`
- `snakefile` — path to the `.smk` file (relative to repo root)
- `config_section` — top-level config key for module parameters
- `activation` — how the module is activated (see below)
- `depends_on` — list of module names that must also be enabled
- `terminal_outputs` — output file patterns used to build `rule all` targets
- `resources` — thread/memory/time defaults
- `validation` — parameter constraints (type, range, allowed values)

### Activation Modes

| Mode | Behavior | Modules |
|------|----------|---------|
| `always` | Always included regardless of config | `shared_utils` |
| `enabled` | Included when `modules.<name>: true` | `trimming`, `denoising`, `clustering`, `itsx_extraction`, `classification`, `stats`, `global_esv`, `global_otu` |
| `classification_stage` | Included when both denoising and classification are enabled | `utils`, `pseudogene_filtering` |

The `global_esv` and `global_otu` modules have `"stage": "post_classification"` in their registry entries, but their activation mode is `enabled` — they are activated by the `modules.global_esv` / `modules.global_otu` toggles.

### Key Functions

- **`is_module_enabled(config, module_name)`** — Single source of truth for checking module state. Used by the Snakefile, `ConfigManager`, and `schema_builder`.
- **`clustering_enabled(config)`** — Returns whether OTU clustering is active.
- **`should_include_module(config, module_name)`** — Resolves activation mode against config.
- **`validate_module_selection(config)`** — Validates dependency constraints.
- **`resolve_module_include_paths(config, repo_root)`** — Returns ordered list of `.smk` paths to include.
- **`resolve_terminal_targets(config, samples)`** — Builds the final target list for `rule all`.

## Config System

### Defaults — `config/defaults.yaml`

Base defaults for all pipeline parameters. Module toggles live under the `modules:` key:

```yaml
modules:
  trimming: true
  denoising: true
  clustering: false
  itsx_extraction: false
  classification: true
  classification_engine: "rdp"
  pseudogene_filtering: false
  stats: true
  global_esv: false
  global_otu: false
```

Each enabled module has a corresponding top-level section with its parameters (e.g., `trimming:`, `denoising:`, `classification:`).

### Presets — `config/presets/`

Marker-specific override files (14 presets) that layer on top of `defaults.yaml`. Users select a preset to configure primer sequences, marker type, and classifier settings for a specific barcode gene.

### Config Flow

```
CLI:  user_config.yaml → Snakefile (loads directly)
API:  defaults.yaml → preset.yaml → user_overrides → ConfigManager.merge() → ResolvedConfig (frozen)
```

## Common Helpers — `workflow/rules/common.smk`

Included first by the Snakefile. Provides shared functions used across all modules:

| Function | Purpose |
|----------|---------|
| `get_module_config(config, module_name)` | Get parameter dict for a module (e.g., `config["trimming"]`) |
| `get_output_dir(config)` | Resolve pipeline output directory |
| `get_sequences_input(config)` | Resolve input FASTA — routes to centroids if clustering, ITSx output if active, else denoised |
| `get_abundance_table(config)` | Resolve abundance table — `OTU.table` if clustering, else `ESV.table.tmp` |
| `get_classification_engine(config)` | Resolve classifier backend (`rdp` or `sintax`) |

## Snakemake Rule Conventions

### Naming

Rule names use lowercase with underscores:

```python
rule pair_reads:
rule quality_filter:
rule vsearch_denoise:
```

### Shell Blocks

All shell commands follow these patterns:

```python
shell:
    """
    set -euo pipefail
    some_command {input} > {output}
    """
```

- Always start with `set -euo pipefail`
- Use `>` for output (write-only), never `>>` (append)
- Properly quote Snakemake wildcards and config values
- Every rule includes a `log:` directive

### Log Directives

Every rule must include a log directive:

```python
rule example:
    input: ...
    output: ...
    log: "{dir}/logs/module_name/{sample}.log"
    shell: "set -euo pipefail\ncommand {input} > {output} 2> {log}"
```

## Adding a New Module

1. Create `.smk` file(s) in the appropriate `workflow/rules/` subdirectory
2. Add a registry entry in `lib/config/module_registry.py` with:
   - Module metadata (`module` dict)
   - `directory`, `snakefile`, `config_section`
   - `activation` mode (`always`, `enabled`, or `classification_stage`)
   - `depends_on` list
   - `terminal_outputs` patterns
   - `validation` constraints for parameters
3. Add config defaults in `config/defaults.yaml` under the module name
4. Add a toggle under `modules:` in `config/defaults.yaml`
5. Validate with `make test-backend`
