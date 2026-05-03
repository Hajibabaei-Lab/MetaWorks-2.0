"""Python-backed registry for MetaWorks pipeline modules."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, Iterable, Optional


MODULE_ALIASES = {
    "pseudogene": "pseudogene_filtering",
}


MODULE_REGISTRY: Dict[str, Dict[str, Any]] = {
    "shared_utils": {
        "module": {
            "name": "shared_utils",
            "version": "2.0.0",
            "description": "Shared Snakemake runtime helpers.",
            "author": "Hajibabaei Lab",
            "enabled_by_default": True,
        },
        "directory": "utils",
        "snakefile": "workflow/rules/utils/utils.smk",
        "compatible_workflows": ["esv"],
        "internal": True,
        "activation": "always",
        "validation": [],
        "resources": None,
        "depends_on": [],
        "checkpoints": [],
    },
    "trimming": {
        "module": {
            "name": "trimming",
            "version": "2.0.0",
            "description": "Trim linked adapters from paired reads with Cutadapt.",
            "author": "Hajibabaei Lab",
            "enabled_by_default": True,
        },
        "directory": "trimming",
        "snakefile": "workflow/rules/trimming/adapter_trimming.smk",
        "config_section": "trimming",
        "compatible_workflows": ["esv"],
        "activation": "enabled",
        "stage": "trimming",
        "terminal_outputs": ["trimmed/{sample}.fasta.gz"],
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
        "depends_on": [],
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
            "enabled_by_default": True,
        },
        "directory": "denoising",
        "snakefile": "workflow/rules/denoising/denoising.smk",
        "config_section": "denoising",
        "compatible_workflows": ["esv"],
        "activation": "enabled",
        "stage": "denoising",
        "terminal_outputs": ["cat.denoised.nonchimeras"],
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
    "clustering": {
        "module": {
            "name": "clustering",
            "version": "2.0.0",
            "description": "OTU clustering at 97% similarity using VSEARCH cluster_smallmem.",
            "author": "Hajibabaei Lab",
            "enabled_by_default": False,
        },
        "directory": "clustering",
        "snakefile": "workflow/rules/clustering/clustering.smk",
        "config_section": "clustering",
        "compatible_workflows": ["esv"],
        "activation": "enabled",
        "stage": "clustering",
        "terminal_outputs": ["centroids.fasta", "OTU.table"],
        "validation": [
            {
                "parameter": "cluster_id",
                "type": "float",
                "min": 0.5,
                "max": 1.0,
                "description": "Clustering similarity threshold (OTU radius).",
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
        "depends_on": ["denoising"],
        "checkpoints": [
            {
                "name": "clustering_complete",
                "trigger": "OTU clustering complete",
                "output": "{output_dir}/checkpoints/clustering_complete.done",
            }
        ],
    },
    "classification": {
        "module": {
            "name": "classification",
            "version": "2.0.0",
            "description": "Assign taxonomy to ESVs.",
            "author": "Hajibabaei Lab",
            "enabled_by_default": True,
        },
        "directory": "classification",
        "snakefile": "workflow/rules/classification/classifier.smk",
        "config_section": "classification",
        "compatible_workflows": ["esv"],
        "activation": "enabled",
        "stage": "classification",
        "terminal_outputs": [],
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
            "enabled_by_default": False,
        },
        "directory": "pseudogene",
        "snakefile": "workflow/rules/pseudogene/pseudogene.smk",
        "config_section": "pseudogene_filtering",
        "compatible_workflows": ["esv"],
        "activation": "classification_stage",
        "stage": "classification",
        "terminal_outputs": ["taxonomy.csv", "ESV.table"],
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
            "enabled_by_default": True,
        },
        "directory": "stats",
        "snakefile": "workflow/rules/stats/stats.smk",
        "config_section": "stats",
        "compatible_workflows": ["esv"],
        "activation": "enabled",
        "additive": True,
        "terminal_outputs": [
            "stats/R1.stats",
            "stats/R2.stats",
            "stats/paired.stats",
            "stats/trimmed.stats",
        ],
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
            "enabled_by_default": True,
        },
        "directory": "utils",
        "snakefile": "workflow/rules/utils/results.smk",
        "config_section": "utils",
        "compatible_workflows": ["esv"],
        "internal": True,
        "activation": "classification_stage",
        "stage": "classification",
        "terminal_outputs": ["results.csv"],
        "validation": [],
        "resources": {
            "default": {"threads": 1, "memory_mb": 2000, "time_minutes": 30},
            "high_load": {"threads": 2, "memory_mb": 4000, "time_minutes": 45},
        },
        "depends_on": ["pseudogene_filtering"],
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


def is_module_enabled(config: Dict[str, Any], module_name: str) -> bool:
    """Resolve whether a module toggle is enabled for the current config."""
    canonical_name = normalize_module_name(module_name)
    entry = MODULE_REGISTRY.get(canonical_name, {})
    default = (
        entry.get("module", {}).get("enabled_by_default", False)
        if isinstance(entry.get("module"), dict)
        else False
    )
    modules = config.get("modules", {})
    if not isinstance(modules, dict):
        return bool(default)
    value = modules.get(canonical_name, default)
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.strip().lower() in ("1", "true", "yes", "on")
    return bool(value)


def classification_stage_enabled(config: Dict[str, Any]) -> bool:
    """Return whether the classification stage can run."""
    return is_module_enabled(config, "denoising") and is_module_enabled(config, "classification")


def should_include_module(config: Dict[str, Any], module_name: str) -> bool:
    """Return whether a registry module should be included in the workflow."""
    canonical_name = normalize_module_name(module_name)
    entry = MODULE_REGISTRY[canonical_name]
    activation = entry.get("activation", "enabled")

    if activation == "always":
        return True
    if activation == "enabled":
        return is_module_enabled(config, canonical_name)
    if activation == "classification_stage":
        return classification_stage_enabled(config)
    raise ValueError(f"Unknown module activation mode for {canonical_name}: {activation!r}")


def validate_module_selection(config: Dict[str, Any]) -> list[str]:
    """Validate enabled user-facing modules against registry dependencies."""
    errors: list[str] = []

    for module_name, entry in MODULE_REGISTRY.items():
        if entry.get("internal", False):
            continue
        if not is_module_enabled(config, module_name):
            continue
        for dependency in entry.get("depends_on", []):
            dependency_name = normalize_module_name(dependency)
            dependency_entry = MODULE_REGISTRY.get(dependency_name, {})
            if dependency_entry.get("internal", False):
                continue
            if not is_module_enabled(config, dependency_name):
                errors.append(
                    f"modules.{module_name}=true requires modules.{dependency_name}=true"
                )

    return errors


def _resolve_included_module_names(config: Dict[str, Any]) -> list[str]:
    """Resolve included modules in dependency order."""
    included = {
        normalize_module_name(name)
        for name in MODULE_REGISTRY
        if should_include_module(config, name)
    }
    ordered: list[str] = []
    visiting: set[str] = set()
    visited: set[str] = set()

    def visit(module_name: str) -> None:
        if module_name in visited:
            return
        if module_name in visiting:
            raise ValueError(f"Cyclic module dependency detected at {module_name}")

        visiting.add(module_name)
        for dependency in MODULE_REGISTRY[module_name].get("depends_on", []):
            dependency_name = normalize_module_name(dependency)
            if dependency_name in included:
                visit(dependency_name)
        visiting.remove(module_name)
        visited.add(module_name)
        ordered.append(module_name)

    for module_name in MODULE_REGISTRY:
        if module_name not in included:
            continue
        visit(module_name)

    return ordered


def resolve_module_include_paths(
    config: Dict[str, Any],
    repo_root: Optional[str] = None,
) -> list[str]:
    """Resolve the snakefile include paths for the active workflow."""
    base_dir = Path(repo_root).resolve() if repo_root else None
    include_paths: list[str] = []

    for module_name in _resolve_included_module_names(config):
        snakefile = MODULE_REGISTRY[module_name].get("snakefile")
        if not snakefile:
            continue
        if base_dir is not None:
            include_paths.append(str((base_dir / str(snakefile)).resolve()))
        else:
            include_paths.append(str(snakefile))

    return include_paths


def _expand_terminal_output(
    output_dir: str,
    pattern: str,
    samples: Iterable[str],
) -> list[str]:
    """Expand a terminal-output pattern into concrete paths."""
    full_pattern = str(Path(output_dir) / pattern)
    if "{sample}" not in pattern:
        return [full_pattern]
    return [full_pattern.format(sample=sample) for sample in samples]


def resolve_terminal_targets(config: Dict[str, Any], samples: Iterable[str]) -> list[str]:
    """Resolve durable workflow targets for rule all."""
    output_dir = config["pipeline"]["output_dir"]
    samples_list = list(samples)
    stage_rank = {"trimming": 0, "denoising": 1, "clustering": 2, "classification": 3}

    included = _resolve_included_module_names(config)
    highest_stage: Optional[str] = None
    for module_name in included:
        entry = MODULE_REGISTRY[module_name]
        if entry.get("additive", False):
            continue
        stage = entry.get("stage")
        if stage not in stage_rank:
            continue
        if highest_stage is None or stage_rank[stage] > stage_rank[highest_stage]:
            highest_stage = stage

    targets: list[str] = []
    for module_name in included:
        entry = MODULE_REGISTRY[module_name]
        include_outputs = entry.get("additive", False) or entry.get("stage") == highest_stage
        if not include_outputs:
            continue
        for pattern in entry.get("terminal_outputs", []):
            targets.extend(_expand_terminal_output(output_dir, pattern, samples_list))

    deduped: list[str] = []
    seen: set[str] = set()
    for target in targets:
        if target in seen:
            continue
        seen.add(target)
        deduped.append(target)

    if not deduped:
        raise ValueError(
            "No final targets selected. Enable at least one runnable module set "
            "(e.g., modules.stats=true, modules.trimming=true, or "
            "modules.classification=true with modules.denoising=true)."
        )

    return deduped
