# MetaWorks Configuration Guide

This guide explains how to configure MetaWorks 2.0 pipeline runs using the profile-based configuration system. For detailed parameter references, see the [Module Parameters](MODULE_PARAMETERS.md) document.

---

## Table of Contents

1. [Quick Start](#quick-start)
2. [Configuration Overview](#configuration-overview)
3. [Profile System](#profile-system)
4. [Pipeline Data Flow](#pipeline-data-flow)
5. [Troubleshooting](#troubleshooting)
6. [Programmatic Configuration](#programmatic-configuration)

---

## Quick Start

### Option A: Use a Profile (Recommended)

```python
from lib.config.config_manager import ConfigManager

config = ConfigManager.load(
    profile="coi",
    workflow="esv",
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
cat > my_run.yaml << EOF
input:
  fastq_dir: "tests/testing_data"
  sample_source: "folder"
EOF

snakemake --configfile config/defaults.yaml config/presets/coi.yaml my_run.yaml
```

---

## Configuration Overview

The configuration system uses three layers, merged in order:

```
┌──────────────────────────────────────┐
│  User Overrides                      │  ← Your minimal config
└──────────────────────────────────────┘
            ↓ applied on top of
┌──────────────────────────────────────┐
│  Profile Configuration               │  ← Marker-specific settings
└──────────────────────────────────────┘
            ↓ applied on top of
┌──────────────────────────────────────┐
│  Pipeline Defaults                   │  ← config/defaults.yaml
└──────────────────────────────────────┘
```

Later layers override earlier ones. Most users only need to provide input data and choose a profile.

### Layer Responsibilities

| Layer | File(s) | Purpose | Edit? |
|--------|---------|---------|-------|
| **User Overrides** | Your YAML or API request | Your specific choices | Yes |
| **Profile** | `config/presets/*.yaml` | Marker-specific presets | Rarely |
| **Defaults** | `config/defaults.yaml` | Pipeline-wide defaults | Rarely |

---

## Profile System

Profiles are pre-configured settings for common marker genes. They dramatically reduce the amount of configuration needed.

### Available Profiles

#### Ready-to-use (validated with test data)

| Profile | Description | Marker |
|---------|-------------|--------|
| `coi` | COI arthropods | COI |
| `coi_vertebrate` | COI vertebrates | COI |
| `16s` | 16S rRNA bacteria/archaea | 16S |
| `28s` | 28S rRNA | 28S |
| `its` | ITS fungi | ITS |
| `12s` | 12S vertebrates | 12S |

#### Starter templates

| Profile | Description | Marker |
|---------|-------------|--------|
| `12s_fish` | 12S fish eDNA | 12S |
| `16s_vertebrate` | 16S vertebrates | 16S |
| `18s` | 18S eukaryotes | 18S |
| `18s_diatom` | 18S diatoms | 18S |
| `its_plants` | ITS plants | ITS |
| `rbcl` | rbcL (general) | rbcL |
| `rbcl_diatom` | rbcL diatoms | rbcL |
| `rbcl_landplant` | rbcL land plants | rbcL |

### Profile Contents

Each profile contains:

1. **Marker information** — The genetic marker type
2. **Classification settings** — Pre-configured classifier options
3. **Pseudogene filtering** — Appropriate genetic code and HMM settings

Example: `config/presets/coi.yaml`:

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

The Web UI includes a profile selector. Simply:
1. Choose your profile from the dropdown
2. Fill in input directory and run name
3. Submit!

#### API

```python
from lib.config.config_manager import ConfigManager

config = ConfigManager.load(
    profile="coi",
    workflow="esv",
    repo_root="/path/to/MetaWorks-2.0"
)
```

For API/UI use with dictionary overrides:

```python
config = ConfigManager.load_from_dict(
    profile="coi",
    workflow="esv",
    user_overrides={"input": {"fastq_dir": "/data/my_samples"}},
    repo_root="/path/to/MetaWorks-2.0"
)
```

#### CLI

```bash
snakemake \
  --configfile config/defaults.yaml \
  --configfile config/presets/coi.yaml \
  --configfile my_run.yaml
```

### Creating Custom Profiles

Create a new file in `config/presets/`:

```yaml
# config/presets/my_custom.yaml
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

## Pipeline Data Flow

MetaWorks 2.0 supports several related workflow variants. The Mermaid source files below are intended for use in Mermaid Chart or another Mermaid editor so they can be styled and exported separately.

- [Paired-end input mode](diagrams/workflow_input_paired.mmd)
- [Dual-indexed input mode](diagrams/workflow_input_dual_indexed.mmd)
- [Single-read input mode](diagrams/workflow_input_single_read.mmd)
- [Standard ESV workflow](diagrams/workflow_esv.mmd)
- [OTU workflow](diagrams/workflow_otu.mmd)
- [ITS-specific workflow](diagrams/workflow_its.mmd)
- [Global ESV post-analysis workflow](diagrams/workflow_global_esv.mmd)
- [Global OTU post-analysis workflow](diagrams/workflow_global_otu.mmd)
- [RDP classification engine — parallel_rdp.py](diagrams/workflow_rdp_wrapper.mmd)

**Input modes:** Paired-end reads are merged with SeqPrep, then primers are trimmed with Cutadapt using linked adapters from an `adapters.fasta`. Dual-indexed runs also pair reads with SeqPrep but derive per-sample adapter files from an adapter CSV (containing SampleID, Amplicon, Forward, Reverse columns). Single-read mode skips pairing — R1 reads are trimmed directly with Cutadapt `-g`; R2 reads are reverse-complemented first, then trimmed with Cutadapt `-a`.

**ESV workflow:** The core pipeline. Trimmed reads are reformatted, concatenated, dereplicated (VSEARCH `fastx_uniques`), denoised with UNOISE3 (VSEARCH `cluster_unoise`), chimera-checked (VSEARCH `uchime3_denovo`), and assembled into an ESV abundance table (VSEARCH `search_exact`). Taxonomic classification (RDP or SINTAX) produces `results.csv`.

**OTU workflow:** Extends the ESV pipeline. After chimera removal, denoised sequences are clustered at 97% similarity (VSEARCH `cluster_smallmem`) to produce centroid sequences, then an OTU table is built (VSEARCH `usearch_global`). Classification runs on centroid sequences instead of ESVs.

**ITS workflow:** Like the ESV pipeline but inserts ITSx extraction after chimera removal. ITSx extracts the configured ITS region (e.g. ITS2) from denoised sequences, stripping flanking rRNA gene portions. Classification then runs on the ITS-region sequences only.

**Classification engine:** Set by `modules.classification_engine`. RDP uses `parallel_rdp.py` to split the input FASTA into chunks and run concurrent RDP classifier JVM processes, then merges output into `rdp.out.tmp`. SINTAX uses VSEARCH `--sintax` and converts output to the same RDP-like format.

**Global ESV / OTU:** Post-analysis workflows that consume `results.csv` from multiple completed trial directories. Denoised sequences are concatenated across trials, globally dereplicated, then either exact-matched (Global ESV, id=1.0) or clustered at 97% (Global OTU). Trial sequences are mapped to the global database and a `GlobalESV` or `GlobalOTU` column is appended to each trial's `global_results.csv`.

The older PNG images in `docs/images/` are retained as historical MetaWorks v1 figures unless they are removed separately.

---
## Troubleshooting

### "Module X not found"

**Cause:** The module name does not match a key in the module registry.

**Solution:** Ensure the module name matches a key in `lib/config/module_registry.py`. Valid module names are: `trimming`, `denoising`, `clustering`, `itsx_extraction`, `classification`, `pseudogene_filtering`, `stats`, `global_esv`, `global_otu`.

### "No sequences after trimming"

**Cause:** Adapter file does not match your primers, or FASTQ files are in the wrong directory.

**Solution:**
1. Check `trimming.adapters` points to the correct adapter FASTA for your primers
2. Verify FASTQ files are in the directory specified by `input.fastq_dir`
3. Try decreasing `trimming.min_adapter_overlap` or increasing `trimming.error_rate`

### "Classifier path not found"

**Cause:** Custom classifier path does not point to an existing file.

**Solution:** If `classification.rdp.use_custom_classifier: true`, ensure `classification.rdp.classifier_path` points to an existing `.properties` file.

### Validation Errors

**Error:** `preprocessing.quality_score must be between 0 and 40`

**Solution:** Check parameter is within the valid range. See the parameter tables in [Module Parameters](MODULE_PARAMETERS.md) for valid ranges.

### "Pipeline runs but produces no results"

**Checklist:**
1. Adapter file matches your primers
2. FASTQ files are in the correct directory
3. Classifier path is correct (if using custom classifier)
4. Quality thresholds are appropriate for your data
5. At least `modules.trimming`, `modules.denoising`, and `modules.classification` are enabled

### "Profile 'X' not found"

**Cause:** The profile name does not match any file in `config/presets/`.

**Solution:** List available profiles:

```python
from lib.config.config_manager import ConfigManager
manager = ConfigManager("/path/to/MetaWorks-2.0")
profiles = manager.list_available_profiles()
for p in profiles:
    print(f"  {p['name']}: {p['description']}")
```

---

## Programmatic Configuration

```python
from lib.config.config_manager import ConfigManager, load_config

# Full control with overrides
config = ConfigManager.load(
    profile="coi",
    workflow="esv",
    overrides={"input": {"fastq_dir": "/data/samples"}},
    repo_root="/path/to/MetaWorks-2.0"
)

print(config.merged["classification"]["min_confidence"])

# Export for workflow
config_dict = config.export_for_workflow()
```

Or using the convenience function:

```python
config = load_config(profile="16s", workflow="esv")
```

For API/UI use with dictionary input:

```python
config = ConfigManager.load_from_dict(
    profile="its",
    workflow="esv",
    user_overrides={
        "input": {"fastq_dir": "/data/samples"},
        "modules": {"itsx_extraction": True},
    },
    repo_root="/path/to/MetaWorks-2.0"
)
```
