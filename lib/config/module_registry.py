"""Python-backed registry for MetaWorks pipeline modules."""

from __future__ import annotations

from copy import deepcopy
from typing import Any, Dict, Iterable, Optional


MODULE_ALIASES = {
    "pseudogene": "pseudogene_filtering",
}


MODULE_REGISTRY: Dict[str, Dict[str, Any]] = {
    "preprocessing": {
        "module": {
            "name": "preprocessing",
            "version": "2.0.0",
            "description": "Preprocess raw reads with SeqPrep.",
            "author": "Hajibabaei Lab",
        },
        "directory": "preprocessing",
        "snakefile": "modules/preprocessing/rules/preprocessing.smk",
        "config_section": "preprocessing",
        "compatible_workflows": ["esv"],
        "validation": [
            {
                "parameter": "quality_score",
                "type": "integer",
                "min": 0,
                "max": 40,
                "description": "Phred quality score cutoff.",
            },
            {
                "parameter": "min_overlap",
                "type": "integer",
                "min": 10,
                "max": 100,
                "description": "Minimum overlap length.",
            },
            {
                "parameter": "max_mismatch",
                "type": "float",
                "min": 0.0,
                "max": 0.5,
                "description": "Maximum mismatch fraction.",
            },
            {
                "parameter": "min_match",
                "type": "float",
                "min": 0.0,
                "max": 1.0,
                "description": "Minimum matching overlap fraction.",
            },
        ],
        "resources": {
            "default": {"threads": 1, "memory_mb": 2000, "time_minutes": 30},
            "high_load": {"threads": 2, "memory_mb": 4000, "time_minutes": 45},
        },
        "depends_on": [],
        "checkpoints": [
            {
                "name": "preprocessing_complete",
                "trigger": "Read preprocessing complete",
                "output": "{output_dir}/checkpoints/preprocessing_complete.done",
            }
        ],
    },
    "trimming": {
        "module": {
            "name": "trimming",
            "version": "2.0.0",
            "description": "Trim linked adapters from paired reads with Cutadapt.",
            "author": "Hajibabaei Lab",
        },
        "directory": "trimming",
        "snakefile": "modules/trimming/rules/adapter_trimming.smk",
        "config_section": "trimming",
        "compatible_workflows": ["esv"],
        "validation": [
            {
                "parameter": "min_length",
                "type": "integer",
                "min": 10,
                "max": 10000,
                "description": "Minimum sequence length after trimming.",
            },
            {
                "parameter": "error_rate",
                "type": "float",
                "min": 0.0,
                "max": 1.0,
                "description": "Maximum allowed error rate.",
            },
            {
                "parameter": "min_adapter_overlap",
                "type": "integer",
                "min": 1,
                "max": 50,
                "description": "Minimum adapter overlap.",
            },
            {
                "parameter": "max_n_bases",
                "type": "integer",
                "min": 0,
                "max": 100,
                "description": "Maximum number of N bases allowed.",
            },
        ],
        "resources": {
            "default": {"threads": 1, "memory_mb": 4000, "time_minutes": 60},
            "high_load": {"threads": 2, "memory_mb": 8000, "time_minutes": 90},
        },
        "depends_on": ["preprocessing"],
        "checkpoints": [
            {
                "name": "trimming_complete",
                "trigger": "Adapter trimming complete",
                "output": "{output_dir}/checkpoints/trimming_complete.done",
            }
        ],
    },
    "denoising": {
        "module": {
            "name": "denoising",
            "version": "2.0.0",
            "description": "Generate ESVs with VSEARCH denoising and chimera removal.",
            "author": "Hajibabaei Lab",
        },
        "directory": "denoising",
        "snakefile": "modules/denoising/rules/denoising.smk",
        "config_section": "denoising",
        "compatible_workflows": ["esv"],
        "validation": [
            {
                "parameter": "min_cluster_size",
                "type": "integer",
                "min": 1,
                "max": 1000,
                "description": "Minimum reads per ESV cluster.",
            },
            {
                "parameter": "threads",
                "type": "integer",
                "min": 1,
                "max": 32,
                "description": "Number of threads for VSEARCH.",
            },
        ],
        "resources": {
            "default": {"threads": 4, "memory_mb": 8000, "time_minutes": 120},
            "high_load": {"threads": 8, "memory_mb": 16000, "time_minutes": 180},
        },
        "depends_on": ["trimming"],
        "checkpoints": [
            {
                "name": "denoising_complete",
                "trigger": "ESV generation complete",
                "output": "{output_dir}/checkpoints/denoising_complete.done",
            }
        ],
    },
    "classification": {
        "module": {
            "name": "classification",
            "version": "2.0.0",
            "description": "Assign taxonomy to ESVs.",
            "author": "Hajibabaei Lab",
        },
        "directory": "classification",
        "snakefile": "modules/classification/rules/classifier.smk",
        "config_section": "classification",
        "compatible_workflows": ["esv"],
        "backend_variants": ["rdp", "sintax"],
        "validation": [
            {
                "parameter": "marker",
                "type": "string",
                "allowed": [
                    "COI",
                    "16S",
                    "16S_vertebrate",
                    "18S_eukaryota",
                    "18S_diatom",
                    "12S_fish",
                    "12S_vertebrate",
                    "ITS_fungi",
                    "28S_fungi",
                    "rbcL_eukaryota",
                    "rbcL_diatom",
                    "rbcL_landPlant",
                    "ITS_plants",
                ],
                "description": "Marker gene type.",
            },
            {
                "parameter": "memory_gb",
                "type": "integer",
                "min": 1,
                "max": 100,
                "description": "RDP classifier memory allocation in GB.",
            },
            {
                "parameter": "min_confidence",
                "type": "float",
                "min": 0.0,
                "max": 1.0,
                "description": "Minimum confidence threshold.",
            },
            {
                "parameter": "threads",
                "type": "integer",
                "min": 1,
                "max": 32,
                "description": "Backend thread count.",
            },
        ],
        "resources": {
            "default": {"threads": 4, "memory_mb": 20000, "time_minutes": 240},
            "high_load": {"threads": 8, "memory_mb": 40000, "time_minutes": 360},
        },
        "depends_on": ["denoising"],
        "checkpoints": [
            {
                "name": "classification_complete",
                "trigger": "Taxonomic assignment complete",
                "output": "{output_dir}/checkpoints/classification_complete.done",
            }
        ],
    },
    "pseudogene_filtering": {
        "module": {
            "name": "pseudogene_filtering",
            "version": "2.0.0",
            "description": "Filter putative pseudogenes using ORF-based methods.",
            "author": "Hajibabaei Lab",
        },
        "directory": "pseudogene",
        "snakefile": "modules/pseudogene/rules/pseudogene.smk",
        "config_section": "pseudogene_filtering",
        "compatible_workflows": ["esv"],
        "validation": [
            {
                "parameter": "genetic_code",
                "type": "integer",
                "allowed": [
                    1,
                    2,
                    3,
                    4,
                    5,
                    6,
                    9,
                    10,
                    11,
                    12,
                    13,
                    14,
                    21,
                    22,
                    23,
                    24,
                    25,
                    26,
                    27,
                    28,
                    29,
                    30,
                    31,
                    33,
                ],
                "description": "NCBI genetic code table number.",
            },
            {
                "parameter": "min_orf_length",
                "type": "integer",
                "min": 30,
                "max": 10000,
                "description": "Minimum ORF length in nucleotides.",
            },
            {
                "parameter": "method",
                "type": "string",
                "allowed": ["hmm", "orf"],
                "description": "Pseudogene filtering method.",
            },
        ],
        "resources": {
            "default": {"threads": 4, "memory_mb": 8000, "time_minutes": 180},
            "high_load": {"threads": 8, "memory_mb": 16000, "time_minutes": 270},
        },
        "depends_on": ["classification", "denoising"],
        "checkpoints": [
            {
                "name": "pseudogene_complete",
                "trigger": "Pseudogene filtering complete",
                "output": "{output_dir}/checkpoints/pseudogene_complete.done",
            }
        ],
    },
    "stats": {
        "module": {
            "name": "stats",
            "version": "2.0.0",
            "description": "Generate pipeline statistics and summary reports.",
            "author": "Hajibabaei Lab",
        },
        "directory": "stats",
        "snakefile": "modules/stats/rules/stats.smk",
        "config_section": "stats",
        "compatible_workflows": ["esv"],
        "validation": [],
        "resources": {
            "default": {"threads": 1, "memory_mb": 4000, "time_minutes": 60},
            "high_load": {"threads": 2, "memory_mb": 8000, "time_minutes": 90},
        },
        "depends_on": ["trimming"],
        "checkpoints": [
            {
                "name": "stats_complete",
                "trigger": "Statistics generation complete",
                "output": "{output_dir}/checkpoints/stats_complete.done",
            }
        ],
    },
    "utils": {
        "module": {
            "name": "utils",
            "version": "2.0.0",
            "description": "Shared utilities for results consolidation.",
            "author": "Hajibabaei Lab",
        },
        "directory": "utils",
        "snakefile": "modules/utils/rules/results.smk",
        "config_section": "utils",
        "compatible_workflows": ["esv"],
        "internal": True,
        "validation": [],
        "resources": {
            "default": {"threads": 1, "memory_mb": 2000, "time_minutes": 30},
            "high_load": {"threads": 2, "memory_mb": 4000, "time_minutes": 45},
        },
        "depends_on": ["classification"],
        "checkpoints": [
            {
                "name": "utils_complete",
                "trigger": "Utilities processing complete",
                "output": "{output_dir}/checkpoints/utils_complete.done",
            }
        ],
    },
}


def normalize_module_name(module_name: str) -> str:
    """Resolve module aliases to canonical registry keys."""
    return MODULE_ALIASES.get(module_name, module_name)


def get_registered_module_names(include_internal: bool = True) -> list[str]:
    """List registered module names."""
    names = []
    for module_name, module_config in MODULE_REGISTRY.items():
        if not include_internal and module_config.get("internal", False):
            continue
        names.append(module_name)
    return names


def get_module_registry_entry(module_name: str) -> Dict[str, Any]:
    """Return a copy of one module definition."""
    canonical_name = normalize_module_name(module_name)
    if canonical_name not in MODULE_REGISTRY:
        raise KeyError(canonical_name)
    return deepcopy(MODULE_REGISTRY[canonical_name])


def load_module_registry(modules: Optional[Iterable[str]] = None) -> Dict[str, Dict[str, Any]]:
    """Return selected module definitions from the registry."""
    if modules is None:
        modules = get_registered_module_names()

    loaded: Dict[str, Dict[str, Any]] = {}
    for module_name in modules:
        canonical_name = normalize_module_name(module_name)
        loaded[canonical_name] = get_module_registry_entry(canonical_name)
    return loaded
