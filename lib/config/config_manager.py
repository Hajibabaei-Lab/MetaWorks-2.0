"""Configuration manager for loading, merging, and validating configs."""

import os
from pathlib import Path
from typing import Any, Dict, List, Optional

from .loader import ConfigError, load_yaml, merge_configs, validate_parameter
from .module_registry import (
    get_module_registry_entry,
    get_registered_module_names,
    load_module_registry,
    normalize_module_name,
)
from .schemas import ModuleConfig, SystemConfig, UserConfig, ValidationError


class ConfigManager:
    """
    Manages configuration loading, merging, and validation.

    The configuration is loaded in layers:
    1. System configuration (immutable infrastructure settings)
    2. Defaults configuration (pipeline-wide parameter defaults)
    3. Profile configuration (marker-specific settings)
    4. User configuration (user's choices and overrides)

    Merge order: defaults → profile → user_config
    """

    def __init__(self, repo_root: Optional[str] = None):
        """
        Initialize the configuration manager.

        Args:
            repo_root: Path to MetaWorks repository root.
                      If None, uses current working directory.
        """
        if repo_root is None:
            repo_root = os.getcwd()

        self.repo_root = Path(repo_root)
        self.system_config: Optional[Dict[str, Any]] = None
        self.defaults_config: Optional[Dict[str, Any]] = None
        self.profile_config: Optional[Dict[str, Any]] = None
        self.module_configs: Dict[str, Dict[str, Any]] = {}
        self.user_config: Optional[Dict[str, Any]] = None
        self.merged_config: Optional[Dict[str, Any]] = None
        self.current_workflow: str = "esv"
        self.current_profile: Optional[str] = None

        # Default configuration paths
        self.system_config_path = self.repo_root / "config" / "system_config.yaml"
        self.defaults_config_path = self.repo_root / "config" / "defaults.yaml"
        self.profiles_dir = self.repo_root / "config" / "profiles"
        self.modules_dir = self.repo_root / "modules"

    def load_system_config(self, system_config_path: Optional[str] = None) -> Dict[str, Any]:
        """
        Load the system configuration.

        Args:
            system_config_path: Path to system config file.
                               If None, uses default path.

        Returns:
            System configuration dictionary

        Raises:
            ConfigError: If system config cannot be loaded
        """
        if system_config_path is None:
            system_config_path = str(self.system_config_path)

        self.system_config = load_yaml(system_config_path)

        # Validate system config structure
        required_sections = [
            "system",
            "paths",
            "runtime",
            "schedulers",
            "resources",
            "logging",
            "security",
            "monitoring",
            "maintenance",
        ]
        for section in required_sections:
            if section not in self.system_config:
                raise ConfigError(f"System config missing required section: {section}")

        return self.system_config

    def load_defaults_config(self, defaults_config_path: Optional[str] = None) -> Dict[str, Any]:
        """
        Load the pipeline defaults configuration.

        Args:
            defaults_config_path: Path to defaults config file.
                                 If None, uses default path.

        Returns:
            Defaults configuration dictionary

        Raises:
            ConfigError: If defaults config cannot be loaded
        """
        if defaults_config_path is None:
            defaults_config_path = str(self.defaults_config_path)

        self.defaults_config = load_yaml(defaults_config_path)
        return self.defaults_config

    def list_available_profiles(self) -> List[Dict[str, str]]:
        """
        List all available profiles.

        Returns:
            List of profile info dictionaries with 'name' and 'description' keys
        """
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
                # If profile can't be loaded, skip it
                continue

        return sorted(profiles, key=lambda p: p["name"])

    def load_profile(self, profile_name: str) -> Dict[str, Any]:
        """
        Load a profile configuration by name.

        Args:
            profile_name: Name of the profile (e.g., "coi", "16s")

        Returns:
            Profile configuration dictionary

        Raises:
            ConfigError: If profile cannot be loaded
        """
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
        """
        Load module registry entries.

        Args:
            modules: List of module names to load. If None, loads all registered modules.

        Returns:
            Dictionary mapping module names to their configurations

        Raises:
            ConfigError: If a requested module is not registered
        """
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
        """
        Load the user configuration.

        Args:
            user_config_path: Path to user config file

        Returns:
            User configuration dictionary

        Raises:
            ConfigError: If user config cannot be loaded
        """
        self.user_config = load_yaml(user_config_path)
        return self.user_config

    def merge(
        self, 
        workflow: str = "esv",
        output_dir: Optional[str] = None,
    ) -> Dict[str, Any]:
        """
        Merge all configuration layers.

        Merge order: Defaults → Profile → User
        Later layers override earlier ones.

        Args:
            workflow: Workflow type ("esv" or "otu")
            output_dir: Optional output directory override

        Returns:
            Merged configuration dictionary

        Raises:
            ConfigError: If required configs haven't been loaded
        """
        if self.defaults_config is None:
            raise ConfigError("Defaults config not loaded. Call load_defaults_config() first.")

        if not self.module_configs:
            raise ConfigError("Module configs not loaded. Call load_module_configs() first.")

        self.current_workflow = workflow

        # Start with defaults
        merged = dict(self.defaults_config)

        # Merge profile config on top (if loaded)
        if self.profile_config:
            # Remove the profile metadata section before merging
            profile_to_merge = {k: v for k, v in self.profile_config.items() if k != "profile"}
            merged = merge_configs(merged, profile_to_merge)

        # Merge user config on top (if loaded)
        if self.user_config:
            merged = merge_configs(merged, self.user_config)

        # Set workflow-specific values
        if "pipeline" not in merged:
            merged["pipeline"] = {}
        merged["pipeline"]["name"] = workflow
        
        # Set output directory
        if output_dir:
            merged["pipeline"]["output_dir"] = output_dir
        elif "pipeline" not in merged or "output_dir" not in merged.get("pipeline", {}):
            # Default output dir based on workflow type
            merged["pipeline"]["output_dir"] = f"{workflow.upper()}_results"

        self.merged_config = merged
        return merged

    def validate(self) -> List[str]:
        """
        Validate the merged configuration.

        Validates:
        - System config structure
        - Module configs against validation rules
        - User config against schemas
        - Parameter consistency

        Returns:
            List of validation errors (empty if valid)
        """
        errors = []

        if self.merged_config is None:
            errors.append("Configuration not merged. Call merge() first.")
            return errors

        # Validate system config
        if self.system_config:
            try:
                SystemConfig(**self.system_config)
            except Exception as e:
                errors.append(f"System config validation: {str(e)}")

        # Validate module configs
        for module_name, module_config in self.module_configs.items():
            # Check validation rules
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

                # Get parameter value from user config or module default
                param_value = self._get_nested_param(self.merged_config, module_name, param_path)

                if param_value is not None:
                    is_valid, error_msg = validate_parameter(
                        param_value, param_type, min_value, max_value, allowed_values
                    )
                    if not is_valid:
                        errors.append(f"{module_name}.{param_path}: {error_msg}")

            # Validate module config structure
            try:
                ModuleConfig(**module_config)
            except Exception as e:
                errors.append(f"Module {module_name} config validation: {str(e)}")

        modules_section = self.merged_config.get("modules", {})
        if isinstance(modules_section, dict):
            for module_name, module_config in self.module_configs.items():
                if module_config.get("internal", False):
                    continue
                if not self._is_module_enabled(module_name, modules_section):
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
                    if not self._is_module_enabled(dependency, modules_section):
                        errors.append(
                            f"modules.{module_name}=true requires modules.{dependency}=true"
                        )

        # Validate user config
        if self.merged_config is not None:
            try:
                UserConfig(**self.merged_config)
            except Exception as e:
                errors.append(f"Workflow config validation: {str(e)}")

        return errors

    def get_module_config(self, module_name: str) -> Dict[str, Any]:
        """
        Get the merged configuration for a specific module.

        Args:
            module_name: Name of the module

        Returns:
            Merged configuration for the module

        Raises:
            ConfigError: If module doesn't exist or config not merged
        """
        if self.merged_config is None:
            raise ConfigError("Configuration not merged. Call merge() first.")

        canonical_name = normalize_module_name(module_name)

        if canonical_name not in self.module_configs:
            raise ConfigError(f"Module not found: {module_name}")

        config_section = self.module_configs[canonical_name].get("config_section", canonical_name)
        module_section = self.merged_config.get(config_section, {})
        return dict(module_section) if isinstance(module_section, dict) else {}

    def export_for_workflow(self, workflow_type: str = "esv") -> Dict[str, Any]:
        """
        Export configuration in the canonical nested schema for Snakemake.

        The returned dict matches the structure expected by Snakefile_ESV:
        pipeline.*, input.*, modules.* (toggles), and per-module top-level
        sections (trimming.*, denoising.*, classification.*, etc.).

        Args:
            workflow_type: Type of workflow ("esv" or "otu")

        Returns:
            Canonical configuration dictionary for Snakemake

        Raises:
            ConfigError: If config not merged or workflow type invalid
        """
        if self.merged_config is None:
            raise ConfigError("Configuration not merged. Call merge() first.")

        if workflow_type not in ["esv", "otu"]:
            raise ConfigError(f"Invalid workflow type: {workflow_type}")

        return dict(self.merged_config)

    def _get_nested_param(self, config: Dict[str, Any], module_name: str, param_path: str) -> Any:
        """
        Get a nested parameter value from config.

        Args:
            config: Configuration dictionary
            module_name: Module name
            param_path: Parameter path (e.g., "quality_score")

        Returns:
            Parameter value or None if not found
        """
        canonical_name = normalize_module_name(module_name)
        config_key = canonical_name
        if canonical_name in self.module_configs:
            config_key = self.module_configs[canonical_name].get(
                "config_section", canonical_name
            )

        if config_key in config:
            module_section = config[config_key]
            if isinstance(module_section, dict):
                if param_path in module_section:
                    return module_section[param_path]
                for _, params in module_section.items():
                    if isinstance(params, dict) and param_path in params:
                        return params[param_path]

        return None

    def _is_module_enabled(self, module_name: str, modules_section: Dict[str, Any]) -> bool:
        """Resolve whether a module is enabled from the runtime modules block."""
        value = modules_section.get(module_name)
        if value is None:
            module_metadata = self.module_configs.get(module_name, {}).get("module", {})
            if isinstance(module_metadata, dict):
                return bool(module_metadata.get("enabled_by_default", False))
            return False
        if isinstance(value, bool):
            return value
        if isinstance(value, str):
            return value.strip().lower() in ("1", "true", "yes", "on")
        return bool(value)

    @classmethod
    def load(
        cls,
        profile: str = "coi",
        workflow: str = "esv",
        user_config_path: Optional[str] = None,
        overrides: Optional[Dict[str, Any]] = None,
        repo_root: Optional[str] = None,
    ) -> "ConfigManager":
        """
        Convenience method to load all configurations.

        Args:
            profile: Profile name (e.g., "coi", "16s")
            workflow: Workflow type ("esv" or "otu")
            user_config_path: Path to user config file (optional)
            overrides: Direct parameter overrides (optional)
            repo_root: Path to repository root (optional)

        Returns:
            ConfigManager instance with loaded configs

        Raises:
            ConfigError: If loading or validation fails
        """
        manager = cls(repo_root)

        # Load all configs
        manager.load_system_config()
        manager.load_defaults_config()
        manager.load_module_configs()
        manager.load_profile(profile)

        # Load user config if provided
        if user_config_path:
            manager.load_user_config(user_config_path)

        # Apply direct overrides if provided
        if overrides:
            if manager.user_config is None:
                manager.user_config = {}
            manager.user_config = merge_configs(manager.user_config, overrides)

        # Merge configs
        manager.merge(workflow=workflow)

        # Validate
        errors = manager.validate()
        if errors:
            raise ValidationError(errors)

        return manager

    @classmethod
    def load_from_dict(
        cls,
        profile: str,
        workflow: str,
        user_overrides: Optional[Dict[str, Any]] = None,
        repo_root: Optional[str] = None,
    ) -> "ConfigManager":
        """
        Load configuration from a dictionary of overrides (for API/UI use).

        Args:
            profile: Profile name (e.g., "coi", "16s")
            workflow: Workflow type ("esv" or "otu")
            user_overrides: Dictionary of user overrides
            repo_root: Path to repository root (optional)

        Returns:
            ConfigManager instance with loaded configs

        Raises:
            ConfigError: If loading fails
        """
        manager = cls(repo_root)

        # Load all configs
        manager.load_system_config()
        manager.load_defaults_config()
        manager.load_module_configs()
        
        # Load profile
        manager.load_profile(profile)

        # Set user overrides directly
        if user_overrides:
            manager.user_config = dict(user_overrides)

        # Merge configs
        manager.merge(workflow=workflow)

        return manager


# Convenience function for quick loading
def load_config(
    profile: str = "coi",
    workflow: str = "esv",
    user_config_path: Optional[str] = None,
    overrides: Optional[Dict[str, Any]] = None,
    repo_root: Optional[str] = None,
) -> ConfigManager:
    """
    Load and validate all configurations.

    Args:
        profile: Profile name (e.g., "coi", "16s")
        workflow: Workflow type ("esv" or "otu")
        user_config_path: Path to user config file (optional)
        overrides: Direct parameter overrides (optional)
        repo_root: Path to repository root (optional)

    Returns:
        ConfigManager instance with loaded and validated configs

    Raises:
        ConfigError: If loading or validation fails
        ValidationError: If validation fails
    """
    return ConfigManager.load(
        profile=profile,
        workflow=workflow,
        user_config_path=user_config_path,
        overrides=overrides,
        repo_root=repo_root,
    )
