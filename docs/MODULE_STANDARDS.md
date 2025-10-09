# MetaWorks Module Standards

## Overview

This document defines the standards and conventions for creating modules in the MetaWorks pipeline. Following these standards ensures consistency, maintainability, and interoperability between modules.

## Module Structure

Each module should be organized in its own directory under `modules/` with the following structure:

```
modules/
└── module_name/
    ├── module.yaml           # Module metadata and configuration
    ├── main.smk              # Main Snakemake rules
    ├── README.md             # Module documentation (optional)
    └── envs/                 # Conda environments (optional)
        └── environment.yaml
```

## Module Metadata (module.yaml)

Every module must include a `module.yaml` file that describes:

- Module identification (name, version, author)
- Input/output contracts
- Dependencies
- Parameters
- Resource requirements
- Validation rules

See `modules/module_template.yaml` for the complete schema.

## Snakemake Rules

### Naming Conventions

1. **Rule names**: Use descriptive, lowercase names with underscores
   ```python
   rule pair_reads:
   rule quality_filter:
   rule detect_samples:
   ```

2. **File paths**: Use config-based paths, not hardcoded
   ```python
   # Good
   output: "{dir}/processed/{sample}.fastq"
   
   # Bad
   output: "COI/processed/{sample}.fastq"
   ```

3. **Module-specific config**: Access via module config namespace
   ```python
   MODULE_CONFIG = config.get("modules", {}).get("preprocessing", {})
   quality_threshold = MODULE_CONFIG.get("quality", 20)
   ```

### Input Validation

All modules should validate their inputs:

```python
def validate_inputs():
    required = ["fastq_dir", "dir"]
    for req in required:
        if req not in config:
            raise ValueError(f"Missing required config: {req}")

validate_inputs()
```

### Checkpoints

Use Snakemake checkpoints for dynamic workflows:

```python
checkpoint detect_samples:
    """Detect all sample pairs in fastq_dir"""
    input:
        fastq_dir = config["fastq_dir"]
    output:
        sample_list = "{dir}/metadata/samples.txt"
    run:
        # Discovery logic
        pass
```

### Logging and Benchmarking

Always include log and benchmark directives:

```python
rule process_sample:
    input: ...
    output: ...
    log: "{dir}/logs/processing/{sample}.log"
    benchmark: "{dir}/benchmarks/processing/{sample}.txt"
    shell: """
        command ... 2>&1 | tee {log}
    """
```

## Python Library Integration

Modules should use the `lib/` Python package for complex logic:

```python
# In Snakemake rule
shell: """
    python -m lib.preprocessing.sequence_stats {input} > {output.stats}
"""

# Or in run block
run:
    from lib.preprocessing import SequenceStats
    stats = SequenceStats(input[0])
    result = stats.calculate()
```

### Creating Library Modules

1. Place in appropriate subpackage (`lib/preprocessing/`, `lib/taxonomy/`, etc.)
2. Include docstrings for all classes and functions
3. Provide CLI interface via `if __name__ == "__main__"`
4. Add unit tests

Example structure:

```python
"""Module description"""

from pathlib import Path
from typing import Dict, Union
import sys

class MyProcessor:
    """Processor class with clear purpose"""
    
    def __init__(self, filepath: Union[str, Path]):
        self.filepath = Path(filepath)
    
    def process(self) -> Dict:
        """Process data and return results"""
        pass

def main():
    """CLI interface"""
    if len(sys.argv) < 2:
        sys.exit("Usage: ...")
    
    # Implementation
    pass

if __name__ == "__main__":
    main()
```

## Module Independence

Modules should be:

1. **Self-contained**: Don't depend on specific directory structures outside their control
2. **Configurable**: Accept parameters via config, not hardcoded
3. **Testable**: Can be run independently with test data
4. **Documented**: Clear inputs, outputs, and parameters

## Workflow Composition

Modules are composed into workflows:

```python
# workflows/esv_basic.smk

module preprocessing:
    snakefile: "../modules/preprocessing/main.smk"
    config: config

module trimming:
    snakefile: "../modules/trimming/main.smk"
    config: config

# Use rules from modules
use rule * from preprocessing
use rule * from trimming
```

## Error Handling

1. **Validate inputs** before processing
2. **Use informative error messages** that guide users to solutions
3. **Handle edge cases** gracefully
4. **Clean up** temporary files on failure

```python
try:
    result = process_data(input_file)
except FileNotFoundError as e:
    raise FileNotFoundError(
        f"Input file not found: {input_file}. "
        f"Please check the 'fastq_dir' configuration."
    ) from e
```

## Testing

Each module should have:

1. **Unit tests** for Python code in `lib/`
2. **Integration tests** with sample data
3. **Documented test cases** in module README

Example test structure:

```
tests/
├── unit/
│   └── test_sequence_stats.py
├── integration/
│   └── test_preprocessing_module.py
└── data/
    └── sample.fastq.gz
```

## Documentation

Module README should include:

1. **Purpose**: What the module does
2. **Requirements**: Dependencies and prerequisites
3. **Configuration**: Required and optional parameters
4. **Examples**: Usage examples
5. **Outputs**: Description of output files

## Version Control

1. Follow semantic versioning (MAJOR.MINOR.PATCH)
2. Update version in `module.yaml` when making changes
3. Document breaking changes in module README

## Best Practices

1. **Use native Snakemake features** (checkpoints, temp files, protected files)
2. **Minimize external dependencies** - prefer stdlib when possible
3. **Make it API-friendly** - modules should be callable programmatically
4. **Think about resumability** - use checkpoints for long-running steps
5. **Resource awareness** - specify appropriate threads/memory requirements
6. **Logging is crucial** - capture all output for debugging

## Migration from Legacy Code

When refactoring existing code into modules:

1. Identify logical boundaries (preprocessing, classification, etc.)
2. Extract hardcoded values to config
3. Consolidate duplicate functionality
4. Add input validation
5. Improve error messages
6. Add logging and benchmarking
7. Create unit tests
8. Document the module

## Example: Complete Module

See `modules/preprocessing/` for a complete example module following all standards.
