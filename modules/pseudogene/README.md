# Pseudogene Module

## Overview

The pseudogene module filters putative pseudogenes using ORF (Open Reading Frame) analysis. It implements two strategies for pseudogene detection: longest ORF analysis and HMM-based ORF analysis. The module also includes taxonomic subsetting capabilities to focus on specific groups.

## Features

- **Two pseudogene detection strategies** (longest ORF and HMM-based)
- **Taxonomic subsetting** using grep-based filtering
- **Flexible genetic code support** for different organisms
- **Configurable ORF parameters** (minimum length, start codons, etc.)
- **HMM profile integration** for protein-coding gene analysis
- **Comprehensive logging** and benchmarking
- **Checkpoint support** for resumability
- **Conditional execution** based on configuration

## Requirements

### Software Dependencies
- ORFfinder
- HMMER >= 3.0 (for HMM-based method)
- Python >= 3.8
- grep and awk (for taxonomic subsetting)

### Input Files
- RDP output from classification module: `{dir}/rdp.out.tmp`
- ESV table from denoising module: `{dir}/ESV.table.tmp`

## Configuration

### Required Parameters

```yaml
dir: "output/directory"
pseudogene_filtering: "yes"  # or "no" to skip filtering
```

### Optional Parameters (with defaults)

```yaml
pseudogene_filtering: "no"    # Whether to run pseudogene filtering
removal_type: 1              # 1 for longest ORF, 2 for HMM-based

# ORFfinder parameters
ORFFINDER:
  g: 5                       # Genetic code (5 = invertebrate mitochondrial)
  s: 2                       # Start codon (0=ATG only, 1=ATG+alt, 2=any sense)
  ml: 30                     # Minimum length (nt)
  n: 'true'                  # Ignore nested ORFs
  strand: 'plus'             # Strand to search

# HMM parameters (for removal_type: 2)
hmm: "hmm/bold.hmm"          # Path to HMM profile

# Taxonomic filtering parameters
grep_type: 1                 # 1 for simple grep, 2 for compound grep
taxon1: '-e Arthropoda'      # First taxonomic filter
taxon2: '-v Chordata'        # Second taxonomic filter (for grep_type 2)
```

## Module-Specific Configuration

Override defaults using the `modules` config namespace:

```yaml
modules:
  pseudogene:
    orffinder:
      g: 2                   # Vertebrate mitochondrial code
      ml: 50                 # Longer minimum length
    hmm:
      path: "/path/to/custom.hmm"
    grep:
      type: 1
      taxon1: '-e Chordata'
```

## Outputs

### Primary Outputs
- `{dir}/taxonomy.csv` - Filtered taxonomic assignments
- `{dir}/ESV.table` - Filtered ESV abundance table

### Intermediate Files
- `{dir}/taxon.zotus` - Taxonomic subset identifiers
- `{dir}/chimera.denoised.nonchimeras.taxon` - Taxon-subset sequences
- `{dir}/longest.orfs.fasta` - Longest ORFs (for longest ORF method)
- `{dir}/orfs.fasta.nt.hmm` - Nucleotide ORFs (for HMM method)
- `{dir}/orfs.fasta.aa.hmm` - Amino acid ORFs (for HMM method)
- `{dir}/hmm.txt` - HMM scan results (for HMM method)

### Checkpoints
- `{dir}/checkpoints/pseudogene_complete.done` - Signals completion

### Logs and Benchmarks
- `{dir}/logs/pseudogene.log` - Execution log
- `{dir}/benchmarks/pseudogene.txt` - Performance metrics

## Pseudogene Detection Strategies

### Strategy 1: Longest ORFs (removal_type: 1)
1. Identify ORFs in taxon-subset sequences
2. Extract the longest ORF for each sequence
3. Filter sequences based on ORF length and quality
4. Generate filtered taxonomy and ESV table

**Use case:** General pseudogene filtering for protein-coding genes

### Strategy 2: HMM-based ORFs (removal_type: 2)
1. Identify nucleotide and amino acid ORFs
2. Scan amino acid ORFs against HMM profile
3. Filter based on HMM match quality
4. Generate filtered taxonomy and ESV table

**Use case:** COI marker with BOLD HMM profile (requires hmm/bold.hmm)

## Taxonomic Subsetting

The module supports two grep-based taxonomic filtering approaches:

### Simple Grep (grep_type: 1)
```bash
grep {taxon1} rdp.out.tmp | awk '{print $1}' > taxon.zotus
```

### Compound Grep (grep_type: 2)
```bash
grep {taxon1} rdp.out.tmp | grep {taxon2} | awk '{print $1}' > taxon.zotus
```

## Usage Examples

### Example 1: Longest ORF Method

```yaml
# config.yaml
dir: "COI"
pseudogene_filtering: "yes"
removal_type: 1
ORFFINDER:
  g: 5        # Invertebrate mitochondrial
  ml: 30      # Minimum 30nt
  s: 2        # Any sense start codon
grep_type: 1
taxon1: '-e Arthropoda'
```

### Example 2: HMM-based Method (COI)

```yaml
# config.yaml
dir: "COI"
pseudogene_filtering: "yes"
removal_type: 2
marker: "COI"
ORFFINDER:
  g: 5        # Invertebrate mitochondrial
  ml: 30      # Minimum 30nt
hmm: "hmm/bold.hmm"
```

### Example 3: No Pseudogene Filtering

```yaml
# config.yaml
dir: "results"
pseudogene_filtering: "no"
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

module pseudogene:
    snakefile: "modules/pseudogene/main.smk"
    config: config

use rule * from preprocessing
use rule * from trimming
use rule * from denoising
use rule * from classification
use rule * from pseudogene
```

## Module Import

To use this module in a Snakemake workflow:

```python
module pseudogene:
    snakefile: "modules/pseudogene/main.smk"
    config: config

# Use all rules from the module
use rule * from pseudogene

# Or use specific rules with prefixes
use rule filter_rdp_longest from pseudogene as filter_taxonomy
```

## Validation

The module performs automatic validation:

1. **Required parameters**: Checks that `dir` is provided
2. **Pseudogene filtering validation**: Ensures value is "yes" or "no"
3. **Removal type validation**: Ensures value is 1 or 2 (when filtering enabled)
4. **HMM file validation**: Verifies HMM profile exists (for HMM method)
5. **Dependency check**: Confirms classification and denoising modules ran first

## Error Handling

Common errors and solutions:

### "HMM profile file not found: {path}"
**Solution**: 
- Verify path to HMM profile is correct
- Download required HMM profiles (e.g., BOLD COI HMM)
- Set correct path in configuration

### "Invalid or missing 'removal_type'"
**Solution**:
- Set `removal_type` to 1 (longest ORF) or 2 (HMM-based)
- For HMM method, ensure `marker` is "COI"

### "grep command failed"
**Solution**:
- Check taxonomic filter strings are valid
- Verify RDP output file format
- Ensure grep and awk are available

### ORFfinder errors
**Solution**:
- Verify ORFfinder is properly installed
- Check genetic code parameter is valid
- Ensure input FASTA format is correct

## Integration with Other Modules

### Dependencies
The pseudogene module requires:
- `{dir}/rdp.out.tmp` from classification module
- `{dir}/ESV.table.tmp` from denoising module

### Provides
The pseudogene module creates:
- `{dir}/taxonomy.csv` and `{dir}/ESV.table` for results module

## Performance

Typical resource usage:
- **Memory**: 2-8 GB (depending on method and data size)
- **Time**: 30-180 minutes (HMM method is typically slower)
- **Threads**: 1-4 (most steps are single-threaded)
- **Disk**: 2-5x input size (for intermediate ORF files)

### Method Comparison
- **Longest ORF**: Faster, less memory, good for general filtering
- **HMM-based**: More accurate for specific markers, requires HMM profile

## ORFfinder Parameters Explained

### Genetic Code (`-g`)
- 1 = Standard code (for rbcL)
- 2 = Vertebrate mitochondrial (for COI with vertebrates)
- 5 = Invertebrate mitochondrial (for COI with invertebrates)

### Start Codon (`-s`)
- 0 = ATG only
- 1 = ATG and alternative initiation codons
- 2 = Any sense codon (default)

### Minimum Length (`-ml`)
- Minimum ORF length in nucleotides
- Default: 30 nt (minimum 30)

### Nested ORFs (`-n`)
- 'true' = ignore nested ORFs
- 'false' = include nested ORFs

### Strand (`-strand`)
- 'plus' = forward strand only
- 'minus' = reverse strand only
- 'both' = both strands

## Testing

Test the module independently:

```bash
# For longest ORF method
snakemake \
    --snakefile modules/pseudogene/main.smk \
    --configfile config/config_ESV.yaml \
    --cores 4 \
    --config pseudogene_filtering="yes" removal_type=1 \
    COI/checkpoints/pseudogene_complete.done

# For HMM-based method
snakemake \
    --snakefile modules/pseudogene/main.smk \
    --configfile config/config_ESV.yaml \
    --cores 4 \
    --config pseudogene_filtering="yes" removal_type=2 marker="COI" \
    COI/checkpoints/pseudogene_complete.done

# Without pseudogene filtering
snakemake \
    --snakefile modules/pseudogene/main.smk \
    --configfile config/config_ESV.yaml \
    --cores 4 \
    --config pseudogene_filtering="no" \
    COI/checkpoints/pseudogene_complete.done
```

## Troubleshooting

### Slow HMM processing
- Check HMM profile quality and size
- Verify sufficient memory allocation
- Consider using longest ORF method for initial analysis

### Poor pseudogene detection
- Adjust ORFfinder parameters (minimum length, genetic code)
- Verify taxonomic filters are appropriate
- Consider using different removal method

### Memory issues
- Reduce input dataset size through pre-filtering
- Use single-threaded execution
- Increase available system memory

### Inconsistent results
- Verify genetic code matches target organisms
- Check HMM profile matches marker gene
- Ensure taxonomic filters are specific enough

## Version History

- **2.0.0** (2025-10-28)
  - Initial modular version
  - Added both longest ORF and HMM-based methods
  - Comprehensive validation and error handling
  - Added conditional execution based on configuration
 - Added checkpoint support
  - Added logging and benchmarking

## See Also

- [Module Standards](../../docs/MODULE_STANDARDS.md) - Development guidelines
- [Classification Module](../classification/README.md) - Previous step
- [ORFfinder Documentation](https://www.ncbi.nlm.nih.gov/orffinder/) - ORF prediction tool
- [HMMER Documentation](http://hmmer.org/) - HMM analysis tool
