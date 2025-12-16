# Denoising Module

## Overview

The denoising module generates Exact Sequence Variants (ESVs) using VSEARCH. It performs dereplication, denoising with UNOISE3 algorithm, chimera removal, and ESV table creation. The module supports both global pooling (all samples together) and per-sample processing approaches.

## Features

- **Global or per-sample processing** based on pooling configuration
- **Dereplication** to identify unique sequences with abundance information
- **UNOISE3 denoising** to separate true biological variants from sequencing errors
- **Chimera removal** using de novo detection
- **ESV table creation** mapping sequences to samples
- **Comprehensive logging** and benchmarking
- **Checkpoint support** for resumability

## Requirements

### Software Dependencies
- VSEARCH >= 2.15
- Python >= 3.8
- sed (for header renaming)
- gzip (for compression)

### Input Files
- Trimmed FASTA files from trimming module: `{dir}/{sample}.fasta.tmp`
- Requires `SAMPLES_UNIQUE` list from preprocessing module

## Configuration

### Required Parameters

```yaml
dir: "output/directory"
SAMPLES_UNIQUE: ["sample1", "sample2", ...]  # from preprocessing module
pooling: "Yes"  # or "No"
```

### Optional Parameters (with defaults)

```yaml
VSEARCH_DENOISE:
  minsize: 8    # Minimum abundance to retain cluster

VSEARCH_TABLE:
  t: 4          # Number of threads for table creation
```

## Module-Specific Configuration

Override defaults using the `modules` config namespace:

```yaml
modules:
  denoising:
    vsearch_denoise:
      minsize: 4    # Lower minimum abundance threshold
    vsearch_table:
      t: 8          # More threads for table creation
```

## Outputs

### Primary Outputs
- `{dir}/ESV.table.tmp` - ESV abundance table (sample × ESV matrix)
- `{dir}/cat.denoised.nonchimeras` - Final denoised, non-chimeric sequences

### Intermediate Files
- `{dir}/cat.fasta.gz` - Concatenated and compressed input sequences
- `{dir}/cat.uniques` - Dereplicated unique sequences with abundance tags
- `{dir}/cat.denoised` - Denoised centroid sequences

### Checkpoints
- `{dir}/checkpoints/denoising_complete.done` - Signals completion

### Logs and Benchmarks
- `{dir}/dereplication.log` - VSEARCH dereplication output
- `{dir}/denoising.log` - VSEARCH denoising output
- `{dir}/chimeraRemoval.log` - VSEARCH chimera detection output
- `{dir}/table.log` - ESV table creation output
- `{dir}/logs/{sample}.denoising.log` - Per-sample denoising logs (per-sample mode)
- `{dir}/benchmarks/denoising/{sample}.txt` - Performance metrics

## Processing Modes

### Global Pooling (pooling: "Yes") - Recommended
1. Concatenate all samples into a single file
2. Rename FASTA headers to replace hyphens with underscores
3. Dereplicate globally to identify unique sequences
4. Denoise using UNOISE3 algorithm with minimum abundance threshold
5. Remove chimeras using de novo detection
6. Create ESV table by mapping all original sequences to final ESVs

**Advantages:**
- Better error correction by pooling information across samples
- More accurate ESV identification
- Single ESV set for all samples

**Disadvantages:**
- Higher memory usage
- Longer processing time

### Per-Sample Processing (pooling: "No")
1. Process each sample separately through dereplication and denoising
2. Concatenate denoised sequences from all samples
3. Dereplicate the combined denoised sequences
4. Remove chimeras from the combined set
5. Create ESV table by mapping each sample to the final ESV set
6. Merge per-sample tables into final ESV table

**Advantages:**
- Lower memory usage per step
- Sample-specific analysis preserved

**Disadvantages:**
- Less accurate error correction
- Potential for ESV fragmentation across samples

## Usage Examples

### Example 1: Global Pooling (Recommended)

```yaml
# config.yaml
dir: "COI"
pooling: "Yes"
VSEARCH_DENOISE:
  minsize: 8
VSEARCH_TABLE:
  t: 4
```

### Example 2: Per-Sample Processing

```yaml
# config.yaml
dir: "results"
pooling: "No"
VSEARCH_DENOISE:
  minsize: 4      # Lower threshold for per-sample
VSEARCH_TABLE:
  t: 8            # More threads for table creation
```

### Example 3: Using in a Workflow

```python
# workflow.smk
module preprocessing:
    snakefile: "modules/preprocessing/main.smk"
    config: config

module trimming:
    snakefile: "modules/trimming/main.smk"
    config: config

module denoising:
    snakefile: "modules/denoising/main.smk"
    config: config

use rule * from preprocessing
use rule * from trimming
use rule * from denoising
```

## Module Import

To use this module in a Snakemake workflow:

```python
module denoising:
    snakefile: "modules/denoising/main.smk"
    config: config

# Use all rules from the module
use rule * from denoising

# Or use specific rules with prefixes
use rule denoise_pooled from denoising as denoise_esvs
```

## Validation

The module performs automatic validation:

1. **Required parameters**: Checks that `dir` and `SAMPLES_UNIQUE` are provided
2. **Pooling validation**: Ensures `pooling` is "Yes" or "No"
3. **Minsize validation**: Confirms minimum abundance threshold is positive
4. **Dependency check**: Verifies preprocessing and trimming modules ran first

## Error Handling

Common errors and solutions:

### "minsize must be positive, got: {value}"
**Solution**: Set `VSEARCH_DENOISE.minsize` to a positive integer (typically 2-8)

### "pooling must be 'Yes' or 'No', got: {value}"
**Solution**: Use either "Yes" or "No" for the pooling parameter

### "Memory allocation error" during dereplication
**Solution**: 
- Use per-sample mode instead of pooling
- Increase available memory
- Reduce input data size through pre-filtering

### "Too many chimeras detected"
**Solution**: 
- Check for contamination in samples
- Verify primer specificity
- Consider adjusting experimental conditions

## Integration with Other Modules

### Dependencies
The denoising module requires:
- `SAMPLES_UNIQUE` from preprocessing module
- `{dir}/{sample}.fasta.tmp` files from trimming module

### Provides
The denoising module creates:
- `{dir}/ESV.table.tmp` for results module
- `{dir}/cat.denoised.nonchimeras` for classification module

## Performance

Typical resource usage (global pooling mode):
- **Memory**: 4-8 GB (peak during dereplication)
- **Time**: 30-120 minutes (depending on data size)
- **Threads**: 1-4 (varies by step)
- **Disk**: 2-5x input size (temporary files)

### Per-sample mode:
- **Memory**: 2-4 GB per step (lower peak usage)
- **Time**: May be longer due to multiple processing rounds
- **Better for large datasets** where global pooling is not feasible

## VSEARCH Parameters Explained

### Dereplication (`--fastx_uniques`)
- Combines identical sequences and adds abundance information as `;size=N` tags
- Maintains sequence quality for downstream processing

### Denoising (`--cluster_unoise`)
- Implements UNOISE3 algorithm to separate true variants from errors
- Uses abundance information to identify likely true sequences
- `--minsize`: Minimum abundance threshold for cluster retention

### Chimera Detection (`--uchime3_denovo`)
- De novo chimera detection without reference database
- `--nonchimeras`: Output only non-chimeric sequences
- `--relabel 'Zotu'`: Relabel sequences as Zotu1, Zotu2, etc.

### ESV Table Creation (`--search_exact`)
- Maps all original sequences to final ESVs
- Creates abundance matrix for downstream analysis
- `--otutabout`: Output in tab-separated format

## Testing

Test the module independently:

```bash
# For global pooling
snakemake \
    --snakefile modules/denoising/main.smk \
    --configfile config/config_ESV.yaml \
    --cores 4 \
    --config pooling="Yes" \
    COI/checkpoints/denoising_complete.done

# For per-sample processing
snakemake \
    --snakefile modules/denoising/main.smk \
    --configfile config/config_ESV.yaml \
    --cores 8 \
    --config pooling="No" \
    COI/checkpoints/denoising_complete.done
```

## Troubleshooting

### High memory usage during dereplication
- Ensure sufficient RAM (4GB+ recommended)
- Consider per-sample mode for large datasets
- Check for very large input files

### Long processing time
- Increase threads for table creation (`VSEARCH_TABLE.t`)
- Verify input file quality and size
- Consider using SSD storage for temporary files

### Low ESV recovery
- Lower `minsize` parameter (but not below 2)
- Check for proper adapter removal in previous steps
- Verify sequence quality

### Chimera detection issues
- High chimera rates may indicate contamination
- Low chimera rates are expected for well-prepared samples
- Consider sample-specific chimera detection if needed

## Version History

- **2.0.0** (2025-10-28)
  - Initial modular version
  - Added both pooling and per-sample modes
  - Comprehensive validation and error handling
  - Added checkpoint support
  - Added logging and benchmarking

## See Also

- [Module Standards](../../docs/MODULE_STANDARDS.md) - Development guidelines
- [Trimming Module](../trimming/README.md) - Previous step
- [Classification Module](../classification/README.md) - Next step (planned)
- [VSEARCH Documentation](https://github.com/torognes/vsearch) - Tool documentation
