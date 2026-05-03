"""Configuration manager for loading, merging, and validating configs."""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional

from .loader import ConfigError, load_yaml, merge_configs, validate_parameter
from .module_registry import (
    get_module_registry_entry,
    get_registered_module_names,
    is_module_enabled,
    load_module_registry,
    normalize_module_name,
)
from .schemas import ModuleConfig, UserConfig, ValidationError


@dataclass(frozen=True)
class ResolvedConfig:
    """Immutable, fully-merged configuration ready for Snakemake or API consumption."""

    merged: Dict[str, Any]
    profile: Optional[str]
    workflow: str
    module_configs: Dict[str, Dict[str, Any]]

    def get_module_config(self, module_name: str) -> Dict[str, Any]:
        """Get the merged configuration for a specific module."""
        canonical_name = normalize_module_name(module_name)

        if canonical_name not in self.module_configs:
            raise ConfigError(f"Module not found: {module_name}")

        config_section = self.module_configs[canonical_name].get("config_section", canonical_name)
        module_section = self.merged.get(config_section, {})
        return dict(module_section) if isinstance(module_section, dict) else {}

    def export_for_workflow(self, workflow_type: Optional[str] = None) -> Dict[str, Any]:
        """Export configuration in the canonical nested schema for Snakemake."""
        if workflow_type is not None and workflow_type not in ("esv", "otu"):
            raise ConfigError(f"Invalid workflow type: {workflow_type}")
        return dict(self.merged)

    def validate(self) -> List[str]:
        """Validate the merged configuration. Returns list of errors (empty if valid)."""
        errors: List[str] = []

        for module_name, module_config in self.module_configs.items():
            validation_rules = module_config.get("validation", [])
            for rule in validation_rules:
                if hasattr(rule, "model_dump"):
                    rule_data = rule.model_dump()
                else:
                    rule_data = dict(rule)

                param_path = rule_data.get("parameter", "")
                param_type = rule_data.get("type", "string")
                min_value = rule_data.get("min")
                max_value = rule_data.get("max")
                allowed_values = rule_data.get("allowed")

                param_value = self._get_nested_param(module_name, param_path)

                if param_value is not None:
                    is_valid, error_msg = validate_parameter(
                        param_value, param_type, min_value, max_value, allowed_values
                    )
                    if not is_valid:
                        errors.append(f"{module_name}.{param_path}: {error_msg}")

            try:
                ModuleConfig(**module_config)
            except Exception as e:
                errors.append(f"Module {module_name} config validation: {str(e)}")

        for module_name, module_config in self.module_configs.items():
            if module_config.get("internal", False):
                continue
            if not is_module_enabled(self.merged, module_name):
                continue
            for dependency in module_config.get("depends_on", []):
                if dependency not in self.module_configs:
                    errors.append(
                        f"Module {module_name} depends on unknown module {dependency}"
                    )
                    continue
                dependency_config = self.module_configs[dependency]
                if dependency_config.get("internal", False):
                    continue
                if not is_module_enabled(self.merged, dependency):
                    errors.append(
                        f"modules.{module_name}=true requires modules.{dependency}=true"
                    )

        try:
            UserConfig(**self.merged)
        except Exception as e:
            errors.append(f"Workflow config validation: {str(e)}")

        return errors

    def _get_nested_param(self, module_name: str, param_path: str) -> Any:
        canonical_name = normalize_module_name(module_name)
        config_key = canonical_name
        if canonical_name in self.module_configs:
            config_key = self.module_configs[canonical_name].get(
                "config_section", canonical_name
            )

        if config_key in self.merged:
            module_section = self.merged[config_key]
            if isinstance(module_section, dict):
                if param_path in module_section:
                    return module_section[param_path]
                for _, params in module_section.items():
                    if isinstance(params, dict) and param_path in params:
                        return params[param_path]

        return None


class ConfigManager:
    """
    Builder for resolved pipeline configurations.

    The configuration is loaded in layers:
    1. Defaults configuration (pipeline-wide parameter defaults)
    2. Profile configuration (marker-specific settings)
    3. User configuration (user's choices and overrides)

    Merge order: defaults → profile → user_config

    Call ``merge()`` to produce an immutable ``ResolvedConfig``.
    """

    def __init__(self, repo_root: Optional[str] = None):
        if repo_root is None:
            repo_root = os.getcwd()

        self.repo_root = Path(repo_root)
        self.defaults_config: Optional[Dict[str, Any]] = None
        self.profile_config: Optional[Dict[str, Any]] = None
        self.module_configs: Dict[str, Dict[str, Any]] = {}
        self.user_config: Optional[Dict[str, Any]] = None
        self.current_profile: Optional[str] = None

        self.defaults_config_path = self.repo_root / "config" / "defaults.yaml"
        self.profiles_dir = self.repo_root / "config" / "presets"
        self.modules_dir = self.repo_root / "modules"

    def load_defaults_config(self, defaults_config_path: Optional[str] = None) -> Dict[str, Any]:
        if defaults_config_path is None:
            defaults_config_path = str(self.defaults_config_path)

        self.defaults_config = load_yaml(defaults_config_path)
        return self.defaults_config

    def list_available_profiles(self) -> List[Dict[str, str]]:
        profiles = []
        
        if not self.profiles_dir.exists():
            return profiles

        for profile_file in self.profiles_dir.glob("*.yaml"):
            try:
                profile_data = load_yaml(str(profile_file))
                profile_info = profile_data.get("profile", {})
                profiles.append({
                    "name": profile_info.get("name", profile_file.stem),
                    "description": profile_info.get("description", ""),
                    "marker": profile_info.get("marker", ""),
                    "file": profile_file.stem,
                })
            except Exception:
                continue

        return sorted(profiles, key=lambda p: p["name"])

    def load_profile(self, profile_name: str) -> Dict[str, Any]:
        profile_path = self.profiles_dir / f"{profile_name}.yaml"
        
        if not profile_path.exists():
            available = [p["name"] for p in self.list_available_profiles()]
            raise ConfigError(
                f"Profile '{profile_name}' not found. "
                f"Available profiles: {', '.join(available) if available else 'none'}"
            )

        self.profile_config = load_yaml(str(profile_path))
        self.current_profile = profile_name
        return self.profile_config

    def load_module_configs(self, modules: Optional[List[str]] = None) -> Dict[str, Dict[str, Any]]:
        try:
            if modules is None:
                self.module_configs = load_module_registry()
            else:
                self.module_configs = {}
                for module_name in modules:
                    canonical_name = normalize_module_name(module_name)
                    self.module_configs[canonical_name] = get_module_registry_entry(canonical_name)
        except KeyError as e:
            requested = str(e).strip("'")
            available = ", ".join(get_registered_module_names())
            raise ConfigError(
                f"Unknown module '{requested}'. Available modules: {available}"
            ) from e

        return self.module_configs

    def load_user_config(self, user_config_path: str) -> Dict[str, Any]:
        self.user_config = load_yaml(user_config_path)
        return self.user_config

    def merge(
        self, 
        workflow: str = "esv",
        output_dir: Optional[str] = None,
    ) -> ResolvedConfig:
        """
        Merge all configuration layers and return an immutable ResolvedConfig.

        Merge order: Defaults → Profile → User.
        Later layers override earlier ones.

        Raises:
            ConfigError: If required configs haven't been loaded
        """
        if self.defaults_config is None:
            raise ConfigError("Defaults config not loaded. Call load_defaults_config() first.")

        if not self.module_configs:
            raise ConfigError("Module configs not loaded. Call load_module_configs() first.")

        merged = dict(self.defaults_config)

        if self.profile_config:
            profile_to_merge = {k: v for k, v in self.profile_config.items() if k != "profile"}
            merged = merge_configs(merged, profile_to_merge)

        if self.user_config:
            merged = merge_configs(merged, self.user_config)

        if "pipeline" not in merged:
            merged["pipeline"] = {}
        merged["pipeline"]["name"] = workflow
        
        if output_dir:
            merged["pipeline"]["output_dir"] = output_dir
        elif "pipeline" not in merged or "output_dir" not in merged.get("pipeline", {}):
            merged["pipeline"]["output_dir"] = f"{workflow.upper()}_results"

        return ResolvedConfig(
            merged=merged,
            profile=self.current_profile,
            workflow=workflow,
            module_configs=dict(self.module_configs),
        )

    @classmethod
    def load(
        cls,
        profile: str = "coi",
        workflow: str = "esv",
        user_config_path: Optional[str] = None,
        overrides: Optional[Dict[str, Any]] = None,
        repo_root: Optional[str] = None,
    ) -> ResolvedConfig:
        """Load all configs, merge, validate, and return a ResolvedConfig."""
        manager = cls(repo_root)

        manager.load_defaults_config()
        manager.load_module_configs()
        manager.load_profile(profile)

        if user_config_path:
            manager.load_user_config(user_config_path)

        if overrides:
            if manager.user_config is None:
                manager.user_config = {}
            manager.user_config = merge_configs(manager.user_config, overrides)

        resolved = manager.merge(workflow=workflow)

        errors = resolved.validate()
        if errors:
            raise ValidationError(errors)

        return resolved

    @classmethod
    def load_from_dict(
        cls,
        profile: str,
        workflow: str,
        user_overrides: Optional[Dict[str, Any]] = None,
        repo_root: Optional[str] = None,
    ) -> ResolvedConfig:
        """Load configuration from a dictionary of overrides (for API/UI use)."""
        manager = cls(repo_root)

        manager.load_defaults_config()
        manager.load_module_configs()
        manager.load_profile(profile)

        if user_overrides:
            manager.user_config = dict(user_overrides)

        return manager.merge(workflow=workflow)


def load_config(
    profile: str = "coi",
    workflow: str = "esv",
    user_config_path: Optional[str] = None,
    overrides: Optional[Dict[str, Any]] = None,
    repo_root: Optional[str] = None,
) -> ResolvedConfig:
    """Load and validate all configurations."""
    return ConfigManager.load(
        profile=profile,
        workflow=workflow,
        user_config_path=user_config_path,
        overrides=overrides,
        repo_root=repo_root,
    )
