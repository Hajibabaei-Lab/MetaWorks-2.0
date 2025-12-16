# Classification Module

## Overview

The classification module performs taxonomic assignment of ESVs using the RDP Classifier. It leverages parallel processing to efficiently classify large numbers of sequences against reference databases, supporting both custom-trained and built-in classifiers.

## Features

- **Parallel processing** using the `parallel_rdp.py` script for efficient classification
- **Custom and built-in classifier support** for different marker genes
- **Memory management** with configurable memory allocation
- **Comprehensive logging** and benchmarking
- **Checkpoint support** for resumability
- **Input validation** with informative error messages

## Requirements

### Software Dependencies
- RDP Classifier >= 2.12
- Python >= 3.8
- BioPython >= 1.70
- Java (required by RDP Classifier)

### Input Files
- Denoised sequences from denoising module: `{dir}/cat.denoised.nonchimeras`

## Configuration

### Required Parameters

```yaml
dir: "output/directory"
RDP:
  custom: "yes" # or "no"
```

### Optional Parameters (with defaults)

```yaml
RDP:
  memory: "20g"              # Memory allocation for RDP classifier
  t: "path/to/rRNAClassifier.properties"  # Properties file for custom classifier
  c: 0                      # Parameter for 16S built-in classifier
  f: "fixrank"              # Parameter for 16S built-in classifier
 g: "fungallsu"            # Parameter for fungal built-in classifier
marker: "COI"               # Marker gene type (for classifier options)
```

## Module-Specific Configuration

Override defaults using the `modules` config namespace:

```yaml
modules:
  classification:
    rdp:
      memory: "10g"         # Lower memory allocation
      custom: "yes"         # Use custom classifier
      t: "/path/to/custom/classifier.properties"
```

## Outputs

### Primary Outputs
- `{dir}/rdp.out.tmp` - Raw RDP classifier output with taxonomic assignments

### Checkpoints
- `{dir}/checkpoints/classification_complete.done` - Signals completion

### Logs and Benchmarks
- `{dir}/logs/classification.log` - Parallel RDP execution log
- `{dir}/benchmarks/classification.txt` - Performance metrics

## Classifier Types

### Custom Classifier (custom: "yes")
Uses a custom-trained RDP classifier with a specific properties file:

```yaml
RDP:
  custom: "yes"
  t: "/path/to/custom/rRNAClassifier.properties"
```

**Use case:** Marker-specific classifiers like COI, 12S, rbcL, etc.

### Built-in Classifiers (custom: "no")

#### 16S Classifier
```yaml
marker: "16S"
RDP:
  custom: "no"
  c: 0
  f: "fixrank"
```

#### Fungal 28S Classifier
```yaml
marker: "28S_fungi"
RDP:
  custom: "no"
 g: "fungallsu"
```

## Parallel Processing

The module uses the `parallel_rdp.py` script to:

1. **Split input FASTA** into chunks based on thread count
2. **Run RDP Classifier in parallel** on each chunk
3. **Concatenate results** into final output file
4. **Clean up temporary files**

### Performance Parameters
- `threads`: Number of parallel processes (default: 4)
- `memory`: Memory allocation per process (default: 20g)
- `timeout`: Maximum time per chunk (default: 12 hours)

## Usage Examples

### Example 1: Custom COI Classifier

```yaml
# config.yaml
dir: "COI"
marker: "COI"
RDP:
  custom: "yes"
  memory: "20g"
  t: "/path/to/coi/classifier/rRNAClassifier.properties"
```

### Example 2: Built-in 16S Classifier

```yaml
# config.yaml
dir: "16S_analysis"
marker: "16S"
RDP:
  custom: "no"
  memory: "10g"
  c: 0
  f: "fixrank"
```

### Example 3: Built-in Fungal Classifier

```yaml
# config.yaml
dir: "fungi"
marker: "28S_fungi"
RDP:
  custom: "no"
  memory: "15g"
  g: "fungallsu"
```

### Example 4: Using in a Workflow

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

module classification:
    snakefile: "modules/classification/main.smk"
    config: config

use rule * from preprocessing
use rule * from trimming
use rule * from denoising
use rule * from classification
```

## Module Import

To use this module in a Snakemake workflow:

```python
module classification:
    snakefile: "modules/classification/main.smk"
    config: config

# Use all rules from the module
use rule * from classification

# Or use specific rules with prefixes
use rule taxonomic_assignment from classification as classify_sequences
```

## Validation

The module performs automatic validation:

1. **Required parameters**: Checks that `dir` is provided
2. **Classifier validation**: Ensures `RDP.custom` is "yes" or "no"
3. **Properties file validation**: Verifies custom classifier file exists (when custom="yes")
4. **Dependency check**: Confirms denoising module ran first

## Error Handling

Common errors and solutions:

### "RDP classifier properties file not found: {path}"
**Solution**: 
- For custom classifiers: Verify path to `.properties` file is correct
- For built-in classifiers: Set `RDP.custom` to "no"

### "Memory allocation error" during classification
**Solution**: 
- Reduce memory allocation (e.g., "10g" instead of "20g")
- Increase available system memory
- Reduce thread count

### "Command timed out for {chunk_file}"
**Solution**:
- Increase timeout in `parallel_rdp.py` script
- Reduce number of threads
- Check system resources

### "Error processing {chunk_file}: {error}"
**Solution**:
- Verify FASTA file format is correct
- Check for sequence length issues
- Ensure sufficient disk space

## Integration with Other Modules

### Dependencies
The classification module requires:
- `{dir}/cat.denoised.nonchimeras` from denoising module

### Provides
The classification module creates:
- `{dir}/rdp.out.tmp` for pseudogene filtering or results modules

## Performance

Typical resource usage:
- **Memory**: 10-20 GB (configurable)
- **Time**: 30-240 minutes (depending on sequence count and classifier)
- **Threads**: 4 (configurable)
- **Disk**: 2-3x input size (for temporary chunks)

### Optimization Tips
- Use appropriate memory allocation for your dataset size
- Adjust thread count based on available CPU cores
- For large datasets, consider using SSD storage
- Monitor memory usage to avoid system swapping

## RDP Classifier Parameters Explained

### Custom Classifier (`-t` option)
- Points to `.properties` file containing trained classifier
- Includes reference sequences and taxonomic hierarchy
- Marker-specific training for better accuracy

### Built-in 16S Classifier (`-c` and `-f` options)
- `c`: Confidence threshold (0-1)
- `f`: Fix rank option for consistent taxonomic ranks

### Built-in Fungal Classifier (`-g` option)
- `g`: Classifier type (fungallsu, fungalits_unite, fungalits_warcup)

## Testing

Test the module independently:

```bash
# For custom classifier
snakemake \
    --snakefile modules/classification/main.smk \
    --configfile config/config_ESV.yaml \
    --cores 4 \
    --config RDP.custom="yes" \
    COI/checkpoints/classification_complete.done

# For built-in classifier
snakemake \
    --snakefile modules/classification/main.smk \
    --configfile config/config_ESV.yaml \
    --cores 4 \
    --config RDP.custom="no" marker="16S" \
    COI/checkpoints/classification_complete.done
```

## Troubleshooting

### Slow processing
- Check if enough memory is allocated
- Verify disk I/O performance
- Consider reducing thread count if system is overloaded

### Java memory errors
- Increase Java heap size in RDP configuration
- Ensure sufficient system memory
- Check for memory leaks in parallel processes

### Inconsistent taxonomic assignments
- Verify classifier training data quality
- Check for chimera contamination in input
- Consider using different confidence thresholds

### Failed parallel execution
- Check for proper thread management
- Verify temporary directory permissions
- Ensure BioPython is properly installed

## Version History

- **2.0.0** (2025-10-28)
  - Initial modular version
  - Added parallel processing support
  - Comprehensive validation and error handling
  - Added checkpoint support
  - Added logging and benchmarking

## See Also

- [Module Standards](../../docs/MODULE_STANDARDS.md) - Development guidelines
- [Denoising Module](../denoising/README.md) - Previous step
- [Pseudogene Module](../pseudogene/README.md) - Next step (planned)
- [RDP Classifier Documentation](https://rdp.cme.msu.edu/) - Tool documentation
