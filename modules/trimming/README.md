# Trimming Module

## Overview

The trimming module handles adapter trimming from paired reads using Cutadapt. It removes linked adapter sequences and performs quality trimming, preparing reads for downstream denoising and analysis.

## Features

- **Linked adapter trimming** using Cutadapt's file-based adapter input
- **Quality trimming** from both ends of sequences
- **Minimum length filtering** to remove short sequences
- **N-base filtering** to remove low-quality sequences
- **Automatic compression** of output FASTA files
- **Comprehensive logging** and benchmarking
- **Checkpoint support** for resumability

## Requirements

### Software Dependencies
- Cutadapt >= 4.0
- Python >= 3.8
- gzip (for compression)

### Input Files
- Paired reads from preprocessing module: `{dir}/paired/{sample}.fastq.gz`
- Linked adapter FASTA file with format:
  ```
  >AmpliconName;
  ^FwdPrimerSeq...ReverseComplementedRevPrimerSeq$
  ```

## Configuration

### Required Parameters

```yaml
dir: "output/directory"
CUTADAPT:
  fasta: "path/to/adapters_anchored.fasta"
```

### Optional Parameters (with defaults)

```yaml
CUTADAPT:
  m: 150         # Minimum sequence length after trimming
  q: "20,20"     # Quality cutoffs (5' and 3' ends)
  e: 0.1         # Maximum error rate (10%)
  O: 3           # Minimum adapter overlap
  mn: 3          # Maximum number of N's allowed
  rc: "No"       # Enable reverse complement mode
```

## Module-Specific Configuration

Override defaults using the `modules` config namespace:

```yaml
modules:
  trimming:
    cutadapt:
      m: 200      # Stricter minimum length
      q: "25,25"  # Higher quality threshold
```

## Outputs

### Primary Outputs
- `{dir}/trimmed/{sample}.fasta.gz` - Trimmed and compressed sequences

### Intermediate Files
- `{dir}/trimmed/{sample}.fasta` - Uncompressed trimmed sequences (temporary)

### Checkpoints
- `{dir}/checkpoints/trimming_complete.done` - Signals completion

### Logs and Benchmarks
- `{dir}/logs/trimming/{sample}.log` - Cutadapt output
- `{dir}/logs/trimming/{sample}.gzip.log` - Compression logs
- `{dir}/benchmarks/trimming/{sample}.txt` - Performance metrics

## Adapter File Format

The adapter file should contain linked adapters in anchored format:

```
>COI_BE;
^GGTCAACAAATCATAAAGATATTGG...TANACYTCNGGRTGNCCRAARAAYCA$

>COI_F230R;
^GGDACWGGWTGAACWGTWTAYCCHCC...CCNGAYATRGCNTTYCCNCG$
```

Where:
- `^` indicates 5' anchor (must be at start)
- `$` indicates 3' anchor (must be at end)
- `...` represents the variable region between primers
- Primers should be in the correct orientation for the merged reads

## Usage Examples

### Example 1: Basic Usage

```yaml
# config.yaml
dir: "COI"
CUTADAPT:
  fasta: "testing/adapters_anchored.fasta"
  m: 150
  q: "20,20"
```

### Example 2: Stricter Quality Filtering

```yaml
# config.yaml
dir: "results"
CUTADAPT:
  fasta: "adapters/primers.fasta"
  m: 200      # Longer minimum length
  q: "25,25"  # Higher quality cutoff
  mn: 0       # No N's allowed
  e: 0.05     # Lower error rate
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

use rule * from preprocessing
use rule * from trimming
```

## Module Import

To use this module in a Snakemake workflow:

```python
module trimming:
    snakefile: "modules/trimming/main.smk"
    config: config

# Use all rules from the module
use rule * from trimming

# Or use specific rules with prefixes
use rule trim_linked_adapters from trimming as trimming_trim
```

## Validation

The module performs automatic validation:

1. **Directory existence**: Checks that output directory is specified
2. **Adapter file validation**: Verifies adapter file exists
3. **Sample list validation**: Confirms SAMPLES_UNIQUE is available
4. **Dependency check**: Ensures preprocessing module ran first

## Error Handling

Common errors and solutions:

### "Adapter file not found: {path}"
**Solution**: Check the path to your adapter file in the configuration

### "SAMPLES_UNIQUE not found in config"
**Solution**: Ensure the preprocessing module runs before trimming

### "No sequences retained after trimming"
**Solution**: Check that:
- Adapter sequences are correct
- Quality thresholds aren't too strict
- Minimum length requirement is appropriate for your data

### Cutadapt error: "Adapter not found"
**Solution**: Verify adapter file format with anchored sequences

## Integration with Other Modules

### Dependencies
The trimming module requires:
- `SAMPLES_UNIQUE` from preprocessing module

### Provides
The trimming module creates:
- `{dir}/trimmed/{sample}.fasta.gz` for denoising module

## Performance

Typical resource usage per sample:
- **Memory**: ~4 GB
- **Time**: 10-60 minutes (depending on file size)
- **Threads**: 1 (Cutadapt is single-threaded)
- **Disk**: Output is ~10-50% of input size (after quality filtering)

## Cutadapt Parameters Explained

- **`-m`**: Minimum length - sequences shorter than this are discarded
- **`-q`**: Quality trimming - trim low-quality bases from ends
- **`-e`**: Error rate - maximum allowed mismatches in adapter alignment
- **`-O`**: Minimum overlap - adapter must overlap by at least this many bases
- **`--max-n`**: Maximum N's - sequences with more N's are discarded
- **`--prefix {name}`**: Use adapter name as prefix in output
- **`--discard-untrimmed`**: Remove sequences without adapters

## Testing

Test the module independently:

```bash
snakemake \
    --snakefile modules/trimming/main.smk \
    --configfile config/config_ESV.yaml \
    --cores 4 \
    COI/checkpoints/trimming_complete.done
```

## Troubleshooting

### Low sequence retention rate
- Check adapter sequences are correct
- Verify primer orientation in adapter file
- Consider lowering quality threshold (`q` parameter)
- Check if minimum length is too strict

### High N-content sequences being retained
- Lower the `mn` parameter
- Check quality of input sequences

### Slow processing
- Increase memory allocation
- Check disk I/O performance
- Consider processing in batches

## Version History

- **2.0.0** (2025-10-09)
  - Initial modular version
  - Added comprehensive validation
  - Added checkpoint support
  - Improved error messages
  - Added logging and benchmarking
  - Separated gzip compression step

## See Also

- [Module Standards](../../docs/MODULE_STANDARDS.md) - Development guidelines
- [Preprocessing Module](../preprocessing/README.md) - Previous step
- [Denoising Module](../denoising/README.md) - Next step (planned)
- [Cutadapt Documentation](https://cutadapt.readthedocs.io/) - Tool documentation
