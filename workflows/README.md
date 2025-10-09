# MetaWorks Workflows

This directory contains workflow definitions that compose multiple modules into complete analysis pipelines.

## Available Workflows

### esv_basic.smk (Planned)
Basic ESV (Exact Sequence Variant) generation workflow without pseudogene filtering.

**Modules used:**
- Input validation
- Preprocessing (read pairing)
- Adapter trimming
- Denoising
- Taxonomic classification
- Results formatting

**Use case:** Standard metabarcoding analysis for most markers (ITS, 16S, 18S, rbcL)

### esv_pseudogene.smk (Planned)
ESV workflow with pseudogene filtering for protein-coding genes.

**Modules used:**
- All modules from esv_basic
- Pseudogene filtering (ORF detection + HMM or longest ORF)
- Taxon-specific filtering

**Use case:** COI and other protein-coding markers where pseudogene filtering is needed

## Workflow Structure

Each workflow file:
1. Loads default configurations
2. Imports required modules
3. Defines module compositions
4. Specifies the final output targets

Example:
```python
# workflows/esv_basic.smk

configfile: "../config/pipeline_defaults.yaml"

# Import modules
module preprocessing:
    snakefile: "../modules/preprocessing/main.smk"
    config: config

module trimming:
    snakefile: "../modules/trimming/main.smk"
    config: config

# Use rules from modules
use rule * from preprocessing
use rule * from trimming

# Define workflow completion
rule workflow_complete:
    input:
        results = config["dir"] + "/results.csv"
    output:
        touch(config["dir"] + "/workflow.complete")
```

## Creating Custom Workflows

To create a custom workflow:

1. Create a new `.smk` file in this directory
2. Import the modules you need
3. Configure module parameters
4. Define the workflow output
5. Update the main `Snakefile_ESV` to reference your workflow

See `docs/MODULE_STANDARDS.md` for detailed guidelines.

## Workflow Selection

Workflows are selected in the main Snakefile based on configuration:

```yaml
# config/config_ESV.yaml
workflow: "esv_pseudogene"  # or "esv_basic"
```

## Future Workflows

Planned workflow additions:
- `otu_clustering.smk` - OTU-based analysis
- `esv_denovo.smk` - De novo taxonomy assignment
- `esv_hybrid.smk` - Combined ESV/OTU approach
