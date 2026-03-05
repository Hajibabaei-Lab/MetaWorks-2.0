# MetaWorks Configuration Guide

This guide explains how to use the profile-based configuration system in MetaWorks v2.0.

---

## Table of Contents

1. [Quick Start](#quick-start)
2. [Configuration Overview](#configuration-overview)
3. [Profile System](#profile-system)
4. [Creating Your Config](#creating-your-config)
5. [Configuration Sections](#configuration-sections)
6. [Common Use Cases](#common-use-cases)
7. [Troubleshooting](#troubleshooting)

---

## Quick Start

### Option A: Use a Profile (Recommended)

```python
from lib.config import ConfigManager

# Load with a profile - simplest approach!
config = ConfigManager.load(
    profile="coi",           # Marker-specific profile
    workflow="esv",          # "esv" or "otu"
    repo_root="/path/to/MetaWorks-2.0"
)
```

### Option B: Web UI

1. Select a **Profile** (e.g., "coi" for COI arthropods)
2. Select a **Workflow** (ESV or OTU)
3. Fill in your input directory and run name
4. Click "Submit Run"

### Option C: CLI with Config File

```bash
# Create minimal user config
cat > my_run.yaml << EOF
input:
  fastq_dir: "tests/testing_data"
  sample_source: "folder"
EOF

# Run with profile
snakemake --configfile config/defaults.yaml config/profiles/coi.yaml my_run.yaml
```

---

## Configuration Overview

The configuration system uses four layers, merged in order:

```
┌─────────────────────────────────────┐
│   User Overrides                   │  ← Your minimal config
│   Just what you want to change     │
└─────────────────────────────────────┘
              ↓ merges with
┌─────────────────────────────────────┐
│   Profile Configuration            │  ← Marker-specific settings
│   COI, 16S, ITS, 12S presets       │
└─────────────────────────────────────┘
              ↓ merges with
┌─────────────────────────────────────┐
│   Pipeline Defaults                │  ← Base parameter defaults
│   config/defaults.yaml             │
└─────────────────────────────────────┘
              ↓ merges with
┌─────────────────────────────────────┐
│   System & Module Configs          │  ← Infrastructure & modules
│   System settings, validation      │
└─────────────────────────────────────┘
```

### Layer Responsibilities

| Layer | File(s) | Purpose | Edit? |
|--------|---------|---------|-------|
| **User Overrides** | Your YAML or API request | Your specific choices | ✅ Yes |
| **Profile** | `config/profiles/*.yaml` | Marker-specific presets | Rarely |
| **Defaults** | `config/defaults.yaml` | Pipeline-wide defaults | Rarely |
| **System/Modules** | `config/system_config.yaml`, `modules/*/module_config.yaml` | Infrastructure | No |

---

## Profile System

Profiles are pre-configured settings for common marker genes. They dramatically reduce the amount of configuration needed.

### Available Profiles

| Profile | Description | Marker | Use Case |
|---------|-------------|--------|----------|
| `coi` | COI arthropods | COI | Insect/crustacean metabarcoding |
| `coi_vertebrate` | COI vertebrates | COI | Fish/bird/mammal metabarcoding |
| `16s` | 16S rRNA | 16S | Bacterial/archaeal microbiome |
| `its` | ITS fungi | ITS_fungi | Fungal community analysis |
| `12s` | 12S vertebrates | 12S_vertebrate | eDNA for vertebrates |

### Profile Contents

Each profile contains:

1. **Marker information** - The genetic marker type
2. **Classification settings** - Pre-configured classifier options
3. **Pseudogene filtering** - Appropriate genetic code and HMM settings

Example: `config/profiles/coi.yaml`:

```yaml
profile:
  name: "coi"
  description: "COI metabarcoding for arthropods and invertebrates"
  marker: "COI"

classification:
  marker: "COI"
  rdp:
    use_custom_classifier: true
    classifier_path: "runtime/classifiers/COI.properties"

pseudogene_filtering:
  method: "hmm"
  genetic_code: 5  # Invertebrate mitochondrial
  hmm_profile: "config/hmm/bold.hmm"
  taxon1: "-e Arthropoda"
  taxon2: "-v Chordata"
```

### Using Profiles

#### Web UI

The Web UI now includes a profile selector. Simply:
1. Choose your profile from the dropdown
2. Fill in input directory and run name
3. Submit!

#### API

```python
from lib.config import ConfigManager

# Load with profile
config = ConfigManager.load(
    profile="coi",
    workflow="esv",
    repo_root="/path/to/MetaWorks-2.0"
)

# Add your overrides
config.user_config = {
    "input": {
        "fastq_dir": "/data/my_samples"
    }
}

# Merge and export
config.merge(workflow="esv")
```

#### CLI

```bash
# Layer configs: defaults → profile → your_config
snakemake \
  --configfile config/defaults.yaml \
  --configfile config/profiles/coi.yaml \
  --configfile my_run.yaml
```

### Creating Custom Profiles

Create a new file in `config/profiles/`:

```yaml
# config/profiles/my_custom.yaml
profile:
  name: "my_custom"
  description: "My custom marker configuration"
  marker: "CUSTOM"

classification:
  marker: "CUSTOM"
  rdp:
    use_custom_classifier: true
    classifier_path: "runtime/classifiers/CUSTOM.properties"

pseudogene_filtering:
  method: "orf"
  genetic_code: 1
```

Then use it:

```python
config = ConfigManager.load(profile="my_custom", workflow="esv")
```

---

## Creating Your Config

### Step 1: Pipeline Selection

Tell MetaWorks what to run and where to put results:

```yaml
pipeline:
  name: "esv"              # or "otu"
  output_dir: "COI_results"   # Your output directory
  parallel_jobs: 4            # Samples to process in parallel
```

**Important:**
- `name`: Must be "esv" or "otu"
- `output_dir`: Use short, simple names (no spaces)
- `parallel_jobs`: Depends on your resources (1-32 recommended)

---

### Step 2: Specify Input Data

Tell MetaWorks where your data is:

```yaml
input:
  sample_source: "folder"      # "folder" (auto-detect) or "csv" (manual)
  samples_csv: "samples.csv"    # Required if sample_source="csv"
  fastq_dir: "data/reads"       # Directory with FASTQ files
```

**Sample Source Options:**

| Option | When to Use | Description |
|--------|--------------|-------------|
| `folder` | Most cases | Auto-detects all FASTQ files in directory |
| `csv` | Special cases | Use `samples.csv` to specify exactly which files |

---

### Step 3: Choose Modules

Enable or disable pipeline modules:

```yaml
modules:
  preprocessing: true              # Quality filtering & read merging
  trimming: true                   # Adapter removal
  denoising: true                 # ESV generation
  classification: true              # Taxonomic assignment
  pseudogene_filtering: false        # Pseudogene removal (optional)
  stats: true                     # Statistics & reports
```

**Module Order:**
Modules run in this order:
1. Preprocessing → Trimming → Denoising → Classification
2. Pseudogene filtering (optional) → Stats → Utils

---

### Step 4: Configure Module Parameters

Each module has optional parameters you can override:

#### Preprocessing

```yaml
preprocessing:
  quality_score: 13        # Phred score cutoff (0-40)
  min_overlap: 25          # Minimum read overlap (10-100 bp)
  max_mismatch: 0.02       # Max fraction mismatches (0.0-0.5)
  min_match: 0.90         # Min fraction matching (0.0-1.0)
```

**When to Adjust:**
- Different sequencing platforms (Illumina, IonTorrent, etc.)
- Lower quality data → decrease `quality_score`
- Shorter amplicons → decrease `min_overlap`

---

#### Trimming

```yaml
trimming:
  adapters: "adapters/COI.fasta"   # Adapter sequences file
  min_length: 150                  # Minimum sequence length
  quality_cutoff: "20,20"          # Quality at 5' and 3' ends
  error_rate: 0.1                   # Adapter matching error rate
  enable_rc: true                    # Enable reverse complement
```

**Required:**
- `adapters`: Must point to existing FASTA file

**Common Issues:**
- **No reads after trimming** → Check `adapters` path matches your primers
- **Too few reads** → Decrease `quality_cutoff` or `min_length`

---

#### Denoising

```yaml
denoising:
  pool_samples: true       # Pool all samples (better for rare ESVs)
  min_cluster_size: 8      # Minimum reads per ESV cluster
  threads: 4              # Number of threads
```

**Pooling Strategy:**

| Strategy | When to Use | Pros | Cons |
|----------|--------------|-------|-------|
| `pool_samples: true` | <100 samples, similar libraries | Better rare ESV detection | Slower |
| `pool_samples: false` | >100 samples, diverse libraries | Faster | May miss rare ESVs |

---

#### Classification

```yaml
modules:
  classification: true
  classification_engine: "rdp"    # "rdp" or "sintax" (one per run)

classification:
  marker: "COI"                  # Marker gene type
  min_confidence: 0.8               # Confidence threshold

  # RDP engine settings
  rdp:
    memory_gb: 20                    # Memory for RDP classifier
    use_custom_classifier: true      # Use custom or built-in
    classifier_path: "runtime/classifiers/COI.properties"  # Custom classifier
    builtin_classifier: "fungallsu"  # Built-in choice

  # VSEARCH SINTAX engine settings (optional)
  sintax:
    db_fasta: null                   # SINTAX-formatted DB FASTA (headers include `;tax=...;`)
    cutoff: null                     # Defaults to min_confidence
    threads: 4
```

**Marker Types:**
- `COI` - Cytochrome oxidase I (animals)
- `16S` - 16S rRNA (bacteria, archaea)
- `ITS` - Internal transcribed spacer (fungi)
- `12S` - 12S rRNA (vertebrates)
- `rbcL` - RuBisCO large subunit (plants)

**Classifier Options:**

| Option | When to Use |
|--------|--------------|
| `use_custom_classifier: true` | COI, 12S, rbcL, custom markers |
| `use_custom_classifier: false` | 16S, ITS fungi (built-in classifiers) |

**Multiple tool integration (one per run):**
- Use `modules.classification_engine: "rdp"` (default) for RDP Classifier.
- Use `modules.classification_engine: "sintax"` to classify with VSEARCH SINTAX and convert output into the RDP-like table MetaWorks expects downstream.

Legacy configs may still use `classification.engine`, but `modules.classification_engine` takes precedence when set.

---

#### Pseudogene Filtering (Optional)

```yaml
pseudogene_filtering:
  method: "hmm"                  # "hmm" or "orf"

  # Grep filtering for taxonomy targeting
  grep_type: 1                     # 1=simple, 2=compound
  taxon1: "-e Arthropoda"           # First grep pattern
  taxon2: ""                         # Second grep pattern (compound only)

  hmm_profile: "config/hmm/bold.hmm"  # HMM profile path
  genetic_code: 5                   # Genetic code (2 or 5 for COI)
```

**Methods:**

| Method | Description | Best For |
|--------|-------------|-----------|
| `hmm` | HMM profile scoring | COI arthropoda, well-curated HMMs |
| `orf` | ORF length filtering | Any protein-coding gene |

**Grep Filtering Options:**

| Type | When to Use | Example |
|------|--------------|---------|
| `grep_type: 1` (simple) | Target one taxon | `taxon1: "-e Arthropoda"` includes Arthropoda |
| `grep_type: 2` (compound) | Include/exclude taxa | `taxon1: "-e Metazoa" + `taxon2: "-v Chordata"` (includes Metazoa, excludes Chordata) |

**Genetic Codes:**
- `2` - Vertebrate mitochondrial (vertebrate COI)
- `5` - Invertebrate mitochondrial (invertebrate COI)

---

### Step 5: Output Options

```yaml
output:
  report_type: 1                 # 1=combined CSV, 2=separate files
  include_intermediate: false     # Keep or delete intermediate files
  compress_output: true          # Gzip output files
  html_report: false            # Generate HTML report (not yet implemented)
```

**Report Type:**

| Type | When to Use | Output |
|------|--------------|---------|
| `1` (combined) | <100 samples | Single CSV with all results |
| `2` (separate) | 100+ samples | Separate ESV.table, taxonomy.csv, sequences.fasta |

**HTML Report:** Currently not implemented. This feature is planned for future release.

---

## Configuration Sections

### Complete Example

```yaml
# MetaWorks User Configuration

pipeline:
  name: "esv"
  output_dir: "COI_results"
  parallel_jobs: 4

input:
  sample_source: "folder"
  samples_csv: "samples.csv"
  fastq_dir: "data/reads"

modules:
  preprocessing: true
  trimming: true
  denoising: true
  classification: true
  pseudogene_filtering: true
  stats: true

preprocessing:
  quality_score: 13
  min_overlap: 25

trimming:
  adapters: "adapters/COI.fasta"
  min_length: 150

denoising:
  pool_samples: true
  min_cluster_size: 8

classification:
  marker: "COI"
  memory_gb: 20
  use_custom_classifier: true
  classifier_path: "classifiers/COI.properties"

pseudogene_filtering:
  method: "hmm"
  target_taxon: "Arthropoda"
  hmm_profile: "config/hmm/bold.hmm"
  genetic_code: 5

output:
  report_type: 1
  include_intermediate: false
  compress_output: true
```

**Total: ~50 lines** (vs. ~150 lines in old format!)

---

## Common Use Cases

### Case 1: COI Analysis (Invertebrates)

```yaml
pipeline:
  name: "esv"
  output_dir: "COI_invertebrates"

classification:
  marker: "COI"
  use_custom_classifier: true
  classifier_path: "classifiers/COI.properties"

pseudogene_filtering:
  method: "hmm"
  target_taxon: "Arthropoda"
  genetic_code: 5  # Invertebrate mitochondrial
```

---

### Case 2: COI Analysis (Vertebrates)

```yaml
pipeline:
  name: "esv"
  output_dir: "COI_vertebrates"

classification:
  marker: "COI"
  use_custom_classifier: true
  classifier_path: "classifiers/COI_vertebrates.properties"

pseudogene_filtering:
  method: "hmm"
  target_taxon: "Chordata"
  genetic_code: 2  # Vertebrate mitochondrial
```

---

### Case 3: ITS Fungi

```yaml
pipeline:
  name: "esv"
  output_dir: "ITS_fungi"

classification:
  marker: "ITS_fungi"
  use_custom_classifier: false  # Use built-in
  builtin_classifier: "fungalits_unite"

# No pseudogene filtering needed for ITS
modules:
  pseudogene_filtering: false
```

---

### Case 4: Large Dataset (100+ Samples)

```yaml
pipeline:
  parallel_jobs: 8  # More parallel jobs

denoising:
  pool_samples: false  # Don't pool for speed

output:
  report_type: 2  # Separate files for memory
```

---

### Case 5: Low Quality Data

```yaml
preprocessing:
  quality_score: 10  # Lower threshold

trimming:
  quality_cutoff: "15,15"  # Less stringent
  min_length: 100  # Shorter minimum

denoising:
  min_cluster_size: 4  # More sensitive
```

---

## Troubleshooting

### Validation Errors

**Error:** `preprocessing.quality_score must be between 0 and 40`

**Solution:** Check parameter is within valid range.

---

### File Not Found

**Error:** `Configuration file not found: config/user_config.yaml`

**Solution:** Check file path is correct and file exists.

---

### Module Config Missing

**Error:** `Module config not found: modules/preprocessing/module_config.yaml`

**Solution:** Ensure module directory exists and contains `module_config.yaml`.

---

### Merging Failed

**Error:** `User config missing required section: pipeline`

**Solution:** Ensure user config has required sections: `pipeline`, `input`.

---

### Runtime Errors

**Issue:** Pipeline runs but produces no results

**Checklist:**
1. Adapter file matches your primers
2. FASTQ files are in correct directory
3. Classifier path is correct
4. Quality scores are appropriate for your data

---

## Advanced Usage

### Programmatic Configuration

```python
from lib.config import ConfigManager

# Load and validate
config = ConfigManager.load("config/user_config.yaml")

# Get module config
preproc_config = config.get_module_config("preprocessing")

# Export for workflow
workflow_config = config.export_for_workflow("esv")

# Validate
errors = config.validate()
if errors:
    for error in errors:
        print(f"Error: {error}")
```

---

### Environment Variables

Override config values with environment variables:

```bash
export METAWORKS_PIPELINE_PARALLEL_JOBS=8
export METAWORKS_CLASSIFICATION_MEMORY_GB=32

# Config manager will read these automatically
```

---

## Migration from Old Config

If you have old `config_ESV.yaml` files, migrate them:

```bash
python scripts/migrate_config.py \
  --input config/config_ESV.yaml \
  --output config/user_config.yaml
```

See [MIGRATION_GUIDE.md](MIGRATION_GUIDE.md) for details.

---

## Further Reading

- [CONFIGURATION_DESIGN.md](CONFIGURATION_DESIGN.md) - Architecture and design
- [CONFIGURATION_EXPLANATION.md](CONFIGURATION_EXPLANATION.md) - Detailed parameter explanations
- [MIGRATION_GUIDE.md](MIGRATION_GUIDE.md) - Migrating old configs
- [MODULE_STANDARDS.md](MODULE_STANDARDS.md) - Module development guide
