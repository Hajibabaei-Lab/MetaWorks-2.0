# MetaWorks Module Parameters

This document provides a detailed reference for every configurable module parameter in MetaWorks 2.0. For an overview of the configuration system, profile selection, and pipeline data flow, see the [Configuration Guide](CONFIGURATION_GUIDE.md).

---

## Table of Contents

1. [Input Data](#input-data)
2. [Module Toggles](#module-toggles)
3. [Pipeline Settings](#pipeline-settings)
4. [Module Parameters](#module-parameters)
5. [Complete Example](#complete-example)
6. [Common Use Cases](#common-use-cases)

---

## Input Data

Tell MetaWorks where your data is:

```yaml
input:
  sample_source: "folder"
  samples_csv: "samples.csv"
  fastq_dir: "data/reads"
```

| Option | When to Use | Description |
|--------|--------------|-------------|
| `folder` | Most cases | Auto-detects all FASTQ files in directory |
| `csv` | Special cases | Use `samples.csv` to specify exactly which files |

---

## Module Toggles

Enable or disable pipeline modules using the 10 toggles from `config/defaults.yaml`:

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

Module data flow:
1. Trimming → Denoising → Classification → Results
2. Optional paths: Clustering (OTU mode), ITSx extraction, Pseudogene filtering
3. Post-results: Global ESV, Global OTU (cross-trial analyses)

---

## Pipeline Settings

```yaml
pipeline:
  parallel_jobs: 4
  output_dir: "COI_results"
```

- `parallel_jobs`: Number of samples processed in parallel (depends on your resources, 1-32 recommended).
- `output_dir`: Optional override for the output directory. Defaults to `{WORKFLOW}_results` (e.g., `ESV_results` or `OTU_results`) if not specified.

---

## Module Parameters

### Trimming

```yaml
trimming:
  read_mode: "paired"
  adapter_source: "file"
  adapters: "tests/adapters_anchored.fasta"
  adapter_csv: ""
  primer: ""
  process_as: "F"
  min_length: 150
  quality_cutoff: "20,20"
  error_rate: 0.1
  min_adapter_overlap: 3
  max_n_bases: 3
  enable_rc: true
```

For read-pairing parameters (quality score, overlap, mismatch, match), configure under `preprocessing:` — see below.

### Preprocessing (read-pairing)

```yaml
preprocessing:
  quality_score: 13
  min_overlap: 25
  max_mismatch: 0.02
  min_match: 0.90
```

This is **not** a module toggle — it is a parameter section for read-pairing quality filtering, consumed by the trimming module's adapter_trimming rule.

| Parameter | Range | Effect |
|-----------|-------|--------|
| `quality_score` | 0-40 | Phred score cutoff for read merging |
| `min_overlap` | 10-100 bp | Minimum read overlap for pairing |
| `max_mismatch` | 0.0-0.5 | Maximum fraction mismatches in overlap |
| `min_match` | 0.0-1.0 | Minimum fraction matching in overlap |

**When to adjust:**
- Lower quality data → decrease `quality_score`
- Shorter amplicons → decrease `min_overlap`

### Denoising

```yaml
denoising:
  pool_samples: true       # Pool all samples (better for rare ESVs)
  min_cluster_size: 8      # Minimum reads per ESV cluster
  threads: 4               # Number of threads
```

| Strategy | When to Use | Pros | Cons |
|----------|--------------|-------|-------|
| `pool_samples: true` | <100 samples, similar libraries | Better rare ESV detection | Slower |
| `pool_samples: false` | >100 samples, diverse libraries | Faster | May miss rare ESVs |

### ITSx Extraction

```yaml
itsx_extraction:
  its_part: "ITS2"
  threads: 4
```

Enable with `modules.itsx_extraction: true`. When active, denoised sequences pass through ITSx to extract the specified ITS region before classification.

### Clustering

```yaml
clustering:
  cluster_id: 0.97
  threads: 4
```

Enable with `modules.clustering: true`. Produces OTU clusters instead of ESVs.

### Classification

```yaml
modules:
  classification: true
  classification_engine: "rdp"

classification:
  min_confidence: 0.8
  rdp:
    memory_gb: 20
    use_custom_classifier: true
    classifier_path: null
    builtin_classifier: "fungallsu"
  sintax:
    db_fasta: null
    cutoff: null
    threads: 4
```

Note the correct nested paths: `classification.rdp.use_custom_classifier` (not `classification.use_custom_classifier`).

**Engine selection** (one per run, set via `modules.classification_engine`):
- `"rdp"` — RDP Classifier (default)
- `"sintax"` — VSEARCH SINTAX with output converted to RDP-like table

| Option | When to Use |
|--------|-------------|
| `use_custom_classifier: true` | COI, 12S, rbcL, custom markers |
| `use_custom_classifier: false` | 16S, ITS fungi (built-in classifiers) |

#### Parallel RDP Classification

When using the RDP engine, MetaWorks runs taxonomic assignment through `workflow/scripts/parallel_rdp.py`, a wrapper that parallelises the RDP classifier across multiple JVM processes:

1. The input FASTA is split into *N* chunks (where *N* = the thread count from the active Snakemake profile).
2. Each chunk is dispatched to a separate `rdp_classifier` JVM process via a `ThreadPoolExecutor`.
3. After all chunks finish, results are concatenated into a single `rdp.out.tmp` file.

Two parameters control resource usage:

| Parameter    | Where Set                                      | Default | Effect                    |
|--------------|------------------------------------------------|---------|---------------------------|
| `memory_gb`  | `classification.rdp.memory_gb`                 | 20      | `-Xmx` per JVM            |
| `threads`    | Snakemake profile `set-threads`                | 4       | Number of concurrent JVMs |

**Total JVM memory = `threads x memory_gb`.** With defaults this is 4 x 20 GB = 80 GB.

Memory requirements scale with reference dataset size — larger reference databases need more RAM per JVM. If total memory is constrained on your system, reduce `threads` rather than lowering `memory_gb` below what the reference requires, because an undersized JVM heap can cause classification failures or very slow performance.

Example — reducing to 2 threads on a machine with 48 GB available:

```yaml
# In config/presets/coi.yaml or user overrides
classification:
  rdp:
    memory_gb: 20    # keep 20 GB per JVM for the COI reference
```

Then run Snakemake with `--threads 2` or override via a Snakemake profile:

```yaml
# In workflow/profiles/local/profile.yaml
set-threads:
  taxonomic_assignment: 2
```

### Pseudogene Filtering

```yaml
pseudogene_filtering:
  method: "hmm"                  # "hmm" or "orf"
  grep_type: 1                   # 1=simple, 2=compound
  taxon1: "-e Arthropoda"        # First grep pattern
  taxon2: "-v Chordata"          # Second grep pattern (compound only)
  hmm_profile: "config/hmm/bold.hmm"
  genetic_code: 5                # Genetic code (2 or 5 for COI)
  orf_start_codon: 2
  min_orf_length: 30
  ignore_nested_orfs: true
  strand: "plus"
```

Pseudogene filtering always runs after classification. When `modules.pseudogene_filtering: false`, it passes data through unchanged. When enabled, it applies HMM or ORF filtering to remove putative pseudogenes. Available for all protein-coding markers (COI, 12S, rbcL, etc.).

| Method | Description | Best For |
|--------|-------------|-----------|
| `hmm` | HMM profile scoring | COI arthropoda, well-curated HMMs |
| `orf` | ORF length filtering | Any protein-coding gene |

| Genetic Code | Use For |
|-------------|---------|
| `2` | Vertebrate mitochondrial |
| `5` | Invertebrate mitochondrial |

### Output

```yaml
output:
  report_type: 1                 # 1=combined CSV, 2=separate files
  include_intermediate: false     # Keep or delete intermediate files
  compress_output: true           # Gzip output files
  html_report: false              # Generate HTML report (not yet implemented)
```

| Type | When to Use | Output |
|------|--------------|--------|
| `1` (combined) | <100 samples | Single CSV with all results |
| `2` (separate) | 100+ samples | Separate ESV.table, taxonomy.csv, sequences.fasta |

---

## Complete Example

```yaml
# MetaWorks User Configuration

input:
  sample_source: "folder"
  fastq_dir: "data/reads"

modules:
  trimming: true
  denoising: true
  clustering: false
  itsx_extraction: false
  classification: true
  classification_engine: "rdp"
  pseudogene_filtering: true
  stats: true
  global_esv: false
  global_otu: false

preprocessing:
  quality_score: 13
  min_overlap: 25

trimming:
  adapters: "adapters/COI.fasta"
  min_length: 150
  quality_cutoff: "20,20"

denoising:
  pool_samples: true
  min_cluster_size: 8

classification:
  marker: "COI"
  min_confidence: 0.8
  rdp:
    use_custom_classifier: true
    classifier_path: "runtime/classifiers/COI.properties"

pseudogene_filtering:
  method: "hmm"
  taxon1: "-e Arthropoda"
  taxon2: "-v Chordata"
  hmm_profile: "config/hmm/bold.hmm"
  genetic_code: 5

output:
  report_type: 1
  include_intermediate: false
  compress_output: true
```

---

## Common Use Cases

### Case 1: COI Invertebrates

```yaml
input:
  fastq_dir: "data/reads"

modules:
  pseudogene_filtering: true

classification:
  marker: "COI"
  rdp:
    use_custom_classifier: true
    classifier_path: "runtime/classifiers/COI.properties"

pseudogene_filtering:
  method: "hmm"
  taxon1: "-e Arthropoda"
  taxon2: "-v Chordata"
  hmm_profile: "config/hmm/bold.hmm"
  genetic_code: 5
```

Or simply use the `coi` profile and override just the input directory.

### Case 2: COI Vertebrates

```yaml
input:
  fastq_dir: "data/reads"

classification:
  marker: "COI"
  rdp:
    use_custom_classifier: true
    classifier_path: "runtime/classifiers/COI_vertebrates.properties"

pseudogene_filtering:
  method: "hmm"
  taxon1: "-e Chordata"
  genetic_code: 2
```

Or use the `coi_vertebrate` profile.

### Case 3: ITS Fungi with ITSx

```yaml
input:
  fastq_dir: "data/reads"

modules:
  itsx_extraction: true
  classification_engine: "rdp"

itsx_extraction:
  its_part: "ITS2"

classification:
  rdp:
    use_custom_classifier: false
    builtin_classifier: "fungalits_unite"
```

Or use the `its` profile (which may already include ITSx settings).

### Case 4: Large Dataset (100+ Samples)

```yaml
input:
  fastq_dir: "data/large_dataset"

pipeline:
  parallel_jobs: 8

denoising:
  pool_samples: false

output:
  report_type: 2
```

### Case 5: Low Quality Data

```yaml
preprocessing:
  quality_score: 10

trimming:
  quality_cutoff: "15,15"
  min_length: 100

denoising:
  min_cluster_size: 4
```
