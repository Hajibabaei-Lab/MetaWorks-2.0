# Preprocessing Module

## Overview

The preprocessing module handles the initial processing of paired-end Illumina sequencing data. It pairs forward and reverse reads using SeqPrep, preparing them for downstream adapter trimming and analysis.

## Features

- **Automatic sample detection** from FASTQ directory or CSV file
- **Quality-based read pairing** with configurable overlap parameters
- **Input validation** with informative error messages
- **Checkpoint support** for workflow resumability
- **Comprehensive logging** and benchmarking

## Requirements

### Software Dependencies
- SeqPrep >= 1.3
- Python >= 3.8
- pandas (for CSV-based sample input)

### Input Files
- Paired-end FASTQ files (gzipped)
- File naming pattern: `{sample}_L001_R1_001.fastq.gz` and `{sample}_L001_R2_001.fastq.gz`

## Configuration

### Required Parameters

```yaml
fastq_dir: "path/to/fastq/files"
dir: "output/directory"
sample_source: "folder"  # or "csv"
```

### Optional Parameters

```yaml
# If sample_source is "csv"
samples_csv: "path/to/samples.csv"

# SeqPrep parameters (with defaults)
SEQPREP:
  q: 13      # Phred quality cutoff
  o: 25      # Minimum overlap length (bp)
  m: 0.02    # Maximum fraction of mismatches in overlap
  n: 0.90    # Minimum fraction of matching overlap
```

## Module-Specific Configuration

You can override defaults using the `modules` config namespace:

```yaml
modules:
  preprocessing:
    seqprep:
      q: 20      # Stricter quality threshold
      o: 30      # Longer minimum overlap
```

## Outputs

### Primary Outputs
- `{dir}/paired/{sample}.fastq.gz` - Merged paired-end reads

### Intermediate Files
- `{sample}_R1.out` - Unpaired forward reads (temporary)
- `{sample}_R2.out` - Unpaired reverse reads (temporary)

### Checkpoints
- `{dir}/checkpoints/pairing_complete.done` - Signals completion

### Logs and Benchmarks
- `{dir}/logs/pairing/{sample}.log` - SeqPrep output
- `{dir}/benchmarks/pairing/{sample}.txt` - Performance metrics

## Usage Examples

### Example 1: Folder-based Input

```yaml
# config.yaml
fastq_dir: "testing/testing_data"
dir: "COI"
sample_source: "folder"

SEQPREP:
  q: 13
  o: 25
  m: 0.02
  n: 0.90
```

### Example 2: CSV-based Input

```yaml
# config.yaml
fastq_dir: "data/fastq"
dir: "results"
sample_source: "csv"
samples_csv: "samples.csv"
```

CSV format:
```csv
sample,path
sample1,/path/to/sample1_L001_R{read}_001.fastq.gz
sample2,/path/to/sample2_L001_R{read}_001.fastq.gz
```

### Example 3: Using in a Workflow

```python
# workflow.smk
module preprocessing:
    snakefile: "modules/preprocessing/main.smk"
    config: config

use rule * from preprocessing
```

## Module Import

To use this module in a Snakemake workflow:

```python
module preprocessing:
    snakefile: "modules/preprocessing/main.smk"
    config: config

# Use all rules from the module
use rule * from preprocessing

# Or use specific rules with prefixes
use rule pair_reads from preprocessing as preprocessing_pair_reads
```

## Validation

The module performs automatic validation:

1. **Directory existence**: Checks that `fastq_dir` exists
2. **Sample source validation**: Ensures valid `sample_source` value
3. **CSV validation**: Verifies CSV file exists if using CSV input
4. **Sample detection**: Confirms FASTQ files are found

## Error Handling

Common errors and solutions:

### "Missing required config parameter: fastq_dir"
**Solution**: Add `fastq_dir` to your configuration file

### "No FASTQ files found in {directory}"
**Solution**: Check that files match pattern `*_R1_001.fastq.gz`

### "FASTQ directory not found: {path}"
**Solution**: Verify the directory path exists and is correct

## Integration with Other Modules

The preprocessing module exports:
- `SAMPLES_UNIQUE` - List of detected samples (added to config)

This can be used by downstream modules:

```python
# In another module
expand("{dir}/output/{sample}.txt", sample=config["SAMPLES_UNIQUE"])
```

## Performance

Typical resource usage per sample:
- **Memory**: ~2 GB
- **Time**: 5-30 minutes (depending on file size)
- **Threads**: 1 (SeqPrep is single-threaded)

## Testing

Test the module with provided test data:

```bash
snakemake \
    --snakefile modules/preprocessing/main.smk \
    --configfile config/config_ESV.yaml \
    --cores 4 \
    COI/checkpoints/pairing_complete.done
```

## Version History

- **2.0.0** (2025-10-09)
  - Initial modular version
  - Added comprehensive validation
  - Added checkpoint support
  - Improved error messages
  - Added logging and benchmarking

## See Also

- [Module Standards](../../docs/MODULE_STANDARDS.md) - Development guidelines
- [Trimming Module](../trimming/README.md) - Next step in the pipeline
- [SeqPrep Documentation](https://github.com/jstjohn/SeqPrep) - Tool documentation
