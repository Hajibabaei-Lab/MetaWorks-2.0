"""Build a Galaxy-style JSON schema describing configurable pipeline parameters."""

from __future__ import annotations

import re
import typing
from typing import Any, Dict, List, Optional, get_args, get_origin

from .loader import merge_configs
from .module_registry import MODULE_REGISTRY
from .schemas import (
    ClassificationConfig,
    ClusteringConfig,
    DenoisingConfig,
    ITSxConfig,
    InputConfig,
    ModuleSelection,
    OutputConfig,
    PipelineConfig,
    PseudogeneConfig,
    StatsConfig,
    TrimmingConfig,
)

SECTION_MODELS = {
    "pipeline": PipelineConfig,
    "input": InputConfig,
    "modules": ModuleSelection,
    "trimming": TrimmingConfig,
    "denoising": DenoisingConfig,
    "clustering": ClusteringConfig,
    "itsx_extraction": ITSxConfig,
    "classification": ClassificationConfig,
    "pseudogene_filtering": PseudogeneConfig,
    "stats": StatsConfig,
    "output": OutputConfig,
}

SECTION_LABELS = {
    "pipeline": "Pipeline Settings",
    "input": "Input Settings",
    "modules": "Module Toggles",
    "trimming": "Trimming Parameters (Cutadapt)",
    "denoising": "Denoising Parameters (VSEARCH)",
    "clustering": "OTU Clustering Parameters (VSEARCH)",
    "itsx_extraction": "ITSx Extraction Parameters",
    "classification": "Classification Parameters",
    "pseudogene_filtering": "Pseudogene Filtering",
    "stats": "Statistics",
    "output": "Output Settings",
}

FILE_REF_FIELDS = frozenset(
    {
        "adapters",
        "classifier_path",
        "hmm_profile",
        "db_fasta",
        "fastq_dir",
        "samples_csv",
    }
)

SECTION_ENABLED_BY = {
    "trimming": "modules.trimming",
    "denoising": "modules.denoising",
    "clustering": "modules.clustering",
    "itsx_extraction": "modules.itsx_extraction",
    "classification": "modules.classification",
    "pseudogene_filtering": "modules.pseudogene_filtering",
    "stats": "modules.stats",
}

SYSTEM_MANAGED_FIELDS = frozenset(
    {
        ("pipeline", "name"),
        ("pipeline", "output_dir"),
    }
)

CLASSIFICATION_GROUPS = {
    "rdp": {
        "label": "RDP Classifier",
        "visible_when": "modules.classification_engine == 'rdp'",
    },
    "sintax": {
        "label": "SINTAX Classifier",
        "visible_when": "modules.classification_engine == 'sintax'",
    },
}


def _resolve_python_type(annotation: Any) -> str:
    origin = get_origin(annotation)
    if origin is typing.Union:
        args = [a for a in get_args(annotation) if a is not type(None)]
        if args:
            return _resolve_python_type(args[0])
    if annotation is bool:
        return "boolean"
    if annotation is int:
        return "integer"
    if annotation is float:
        return "float"
    return "string"


def _extract_constraints(field_info: Any) -> Optional[Dict[str, Any]]:
    constraints: Dict[str, Any] = {}
    for meta in getattr(field_info, "metadata", []) or []:
        for attr in ("ge", "le", "gt", "lt", "min_length", "max_length", "multiple_of"):
            val = getattr(meta, attr, None)
            if val is not None:
                constraints[attr] = val
        pattern = getattr(meta, "pattern", None)
        if pattern is not None:
            constraints["pattern"] = str(pattern)
    return constraints or None


def _get_registry_allowed(section_name: str, field_name: str) -> Optional[List[Any]]:
    entry = MODULE_REGISTRY.get(section_name, {})
    for rule in entry.get("validation", []):
        if rule.get("parameter") == field_name:
            allowed = rule.get("allowed")
            if allowed is not None:
                return list(allowed)
    return None


def _infer_type_from_value(value: Any) -> str:
    if isinstance(value, bool):
        return "boolean"
    if isinstance(value, int):
        return "integer"
    if isinstance(value, float):
        return "float"
    return "string"


def _field_type(
    annotation: Any, field_name: str, allowed: Optional[List[Any]]
) -> str:
    if field_name in FILE_REF_FIELDS:
        return "file_ref"
    if allowed is not None:
        return "select"
    return _resolve_python_type(annotation)


def _pattern_to_options(pattern: str) -> Optional[List[str]]:
    match = re.match(r"^\^\((.+)\)\$$", pattern)
    if match:
        return [v.strip() for v in match.group(1).split("|")]
    return None


def _build_scalar_field(
    field_name: str,
    field_info: Any,
    section_name: str,
    effective_default: Any,
) -> Dict[str, Any]:
    annotation = field_info.annotation
    description = getattr(field_info, "description", None)
    constraints = _extract_constraints(field_info)
    allowed = _get_registry_allowed(section_name, field_name)

    if (
        constraints
        and "pattern" in constraints
        and allowed is None
    ):
        options = _pattern_to_options(constraints["pattern"])
        if options:
            allowed = options
            del constraints["pattern"]

    schema_type = _field_type(annotation, field_name, allowed)

    result: Dict[str, Any] = {"type": schema_type}
    if description:
        result["description"] = description
    if constraints:
        result["constraints"] = constraints
    if allowed is not None:
        result["options"] = allowed
    if effective_default is not None or not field_info.is_required():
        result["default"] = effective_default

    return result


def _build_inferred_field(field_name: str, value: Any) -> Dict[str, Any]:
    result: Dict[str, Any] = {"type": _infer_type_from_value(value)}
    result["default"] = value
    return result


def _build_modules_section(merged: Dict[str, Any]) -> Dict[str, Any]:
    model = ModuleSelection
    modules_data = merged.get("modules", {})
    fields: Dict[str, Any] = {}

    for key, value in modules_data.items():
        if key in model.model_fields:
            field_info = model.model_fields[key]
            if key == "classification_engine":
                classification_entry = MODULE_REGISTRY.get("classification", {})
                engine_options = classification_entry.get("backend_variants", ["rdp", "sintax"])
                fields[key] = {
                    "type": "select",
                    "options": engine_options,
                    "default": value,
                    "description": getattr(field_info, "description", None),
                }
            else:
                fields[key] = {
                    "type": "boolean",
                    "default": value,
                }
        else:
            fields[key] = _build_inferred_field(key, value)

    return fields


def _build_classification_group(
    group_name: str,
    group_data: Dict[str, Any],
    group_def: Dict[str, Any],
) -> Dict[str, Any]:
    fields: Dict[str, Any] = {}

    for sub_name, sub_value in group_data.items():
        allowed = _get_registry_allowed("classification", sub_name)

        sub_schema: Dict[str, Any] = {}
        if sub_name in FILE_REF_FIELDS:
            sub_schema["type"] = "file_ref"
        elif allowed is not None:
            sub_schema["type"] = "select"
            sub_schema["options"] = allowed
        elif sub_value is not None:
            sub_schema["type"] = _infer_type_from_value(sub_value)
        else:
            sub_schema["type"] = "string"

        for rule in MODULE_REGISTRY.get("classification", {}).get("validation", []):
            if rule.get("parameter") == sub_name:
                rule_constraints: Dict[str, Any] = {}
                if rule.get("min") is not None:
                    rule_constraints["ge"] = rule["min"]
                if rule.get("max") is not None:
                    rule_constraints["le"] = rule["max"]
                if rule_constraints:
                    sub_schema["constraints"] = rule_constraints
                if rule.get("description"):
                    sub_schema["description"] = rule["description"]
                break

        sub_schema["default"] = sub_value
        fields[sub_name] = sub_schema

    return {
        "type": "group",
        "label": group_def["label"],
        "visible_when": group_def["visible_when"],
        "fields": fields,
    }


def _build_classification_section(merged: Dict[str, Any]) -> Dict[str, Any]:
    model = ClassificationConfig
    class_data = merged.get("classification", {})
    fields: Dict[str, Any] = {}

    for key, value in class_data.items():
        if key in CLASSIFICATION_GROUPS and isinstance(value, dict):
            fields[key] = _build_classification_group(
                key, value, CLASSIFICATION_GROUPS[key]
            )
        elif key in model.model_fields and not isinstance(value, dict):
            field_info = model.model_fields[key]
            fields[key] = _build_scalar_field(
                key, field_info, "classification", value
            )
        elif not isinstance(value, dict):
            fields[key] = _build_inferred_field(key, value)

    return fields


def _build_standard_section(
    section_name: str,
    section_data: Dict[str, Any],
    model_cls: Any,
) -> Dict[str, Any]:
    fields: Dict[str, Any] = {}

    for key, value in section_data.items():
        if (section_name, key) in SYSTEM_MANAGED_FIELDS:
            continue
        if isinstance(value, dict):
            continue
        if key in model_cls.model_fields:
            field_info = model_cls.model_fields[key]
            fields[key] = _build_scalar_field(key, field_info, section_name, value)
        else:
            fields[key] = _build_inferred_field(key, value)

    return fields


def build_config_schema(
    defaults_config: Dict[str, Any],
    profile_config: Dict[str, Any],
    profile_name: str,
    workflow: str = "esv",
) -> Dict[str, Any]:
    profile_to_merge = {k: v for k, v in profile_config.items() if k != "profile"}
    profile_meta = profile_config.get("profile", {})
    marker = profile_meta.get("marker", "")

    merged = merge_configs(dict(defaults_config), profile_to_merge)

    sections: Dict[str, Any] = {}

    for section_name, section_data in merged.items():
        if section_name not in SECTION_MODELS:
            continue
        if not isinstance(section_data, dict):
            continue

        model_cls = SECTION_MODELS[section_name]
        section_schema: Dict[str, Any] = {
            "label": SECTION_LABELS.get(section_name, section_name),
        }

        if section_name in SECTION_ENABLED_BY:
            section_schema["enabled_by"] = SECTION_ENABLED_BY[section_name]
            section_schema["collapsed"] = True

        if section_name == "modules":
            section_schema["fields"] = _build_modules_section(merged)
        elif section_name == "classification":
            section_schema["fields"] = _build_classification_section(merged)
        else:
            section_schema["fields"] = _build_standard_section(
                section_name, section_data, model_cls
            )

        sections[section_name] = section_schema

    return {
        "profile": profile_name,
        "marker": marker,
        "workflow": workflow,
        "sections": sections,
    }
