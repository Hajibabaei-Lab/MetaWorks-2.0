"""Configuration manager for loading, merging, and validating configs."""

import os
from pathlib import Path
from typing import Any, Dict, List, Optional

from .loader import ConfigError, load_yaml, merge_configs, validate_parameter
from .schemas import ModuleConfig, SystemConfig, UserConfig, ValidationError


class ConfigManager:
    """
    Manages configuration loading, merging, and validation.

    The configuration is loaded in three layers:
    1. System configuration (immutable infrastructure settings)
    2. Module configurations (module defaults and metadata)
    3. User configuration (user's choices)

    User config overrides module config, which overrides system config.
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
        self.module_configs: Dict[str, Dict[str, Any]] = {}
        self.user_config: Optional[Dict[str, Any]] = None
        self.merged_config: Optional[Dict[str, Any]] = None

        # Default configuration paths
        self.system_config_path = self.repo_root / "config" / "system_config.yaml"
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

    def load_module_configs(self, modules: Optional[List[str]] = None) -> Dict[str, Dict[str, Any]]:
        """
        Load module configurations.

        Args:
            modules: List of module names to load. If None, loads all modules.

        Returns:
            Dictionary mapping module names to their configurations

        Raises:
            ConfigError: If module config cannot be loaded
        """
        if modules is None:
            # Auto-discover modules
            modules = []
            for module_dir in self.modules_dir.iterdir():
                module_config_file = module_dir / "module_config.yaml"
                if module_config_file.exists():
                    modules.append(module_dir.name)

        self.module_configs = {}
        for module_name in modules:
            module_config_path = self.modules_dir / module_name / "module_config.yaml"

            if not module_config_path.exists():
                raise ConfigError(f"Module config not found: {module_config_path}")

            try:
                config = load_yaml(str(module_config_path))
                self.module_configs[module_name] = config
            except Exception as e:
                raise ConfigError(f"Error loading module config {module_name}: {str(e)}")

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

        # Validate required sections
        required_sections = ["pipeline", "input"]
        for section in required_sections:
            if section not in self.user_config:
                raise ConfigError(f"User config missing required section: {section}")

        return self.user_config

    def merge(self) -> Dict[str, Any]:
        """
        Merge all configuration layers.

        Merge order: System → Module → User
        User config overrides module config, which overrides system config.

        Returns:
            Merged configuration dictionary

        Raises:
            ConfigError: If configs haven't been loaded
        """
        if self.system_config is None:
            raise ConfigError("System config not loaded. Call load_system_config() first.")

        if not self.module_configs:
            raise ConfigError("Module configs not loaded. Call load_module_configs() first.")

        if self.user_config is None:
            raise ConfigError("User config not loaded. Call load_user_config() first.")

        # Start with system config
        merged = self.system_config.copy()

        # Add module configs
        merged["modules"] = {}
        for module_name, module_config in self.module_configs.items():
            merged["modules"][module_name] = module_config

        # Merge user config on top
        merged = merge_configs(merged, self.user_config)

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
        try:
            SystemConfig(**self.merged_config)
        except Exception as e:
            errors.append(f"System config validation: {str(e)}")

        # Validate module configs
        for module_name, module_config in self.module_configs.items():
            # Check validation rules
            validation_rules = module_config.get("validation", [])
            for rule in validation_rules:
                param_path = rule.get("parameter", "")
                param_type = rule.get("type", "string")
                min_value = rule.get("min")
                max_value = rule.get("max")
                allowed_values = rule.get("allowed")

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

        # Validate user config
        try:
            UserConfig(**self.user_config)
        except Exception as e:
            errors.append(f"User config validation: {str(e)}")

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

        if module_name not in self.module_configs:
            raise ConfigError(f"Module not found: {module_name}")

        # Get module defaults
        module_defaults = self.module_configs[module_name].get("parameters", {})

        # Get user overrides
        user_overrides = self.user_config.get(module_name, {})

        # Merge them
        return merge_configs(module_defaults, user_overrides)

    def export_for_workflow(self, workflow_type: str = "esv") -> Dict[str, Any]:
        """
        Export configuration for Snakemake workflow.

        Args:
            workflow_type: Type of workflow ("esv" or "otu")

        Returns:
            Flattened configuration dictionary for Snakemake

        Raises:
            ConfigError: If config not merged or workflow type invalid
        """
        if self.merged_config is None:
            raise ConfigError("Configuration not merged. Call merge() first.")

        if workflow_type not in ["esv", "otu"]:
            raise ConfigError(f"Invalid workflow type: {workflow_type}")

        # Flatten configuration for Snakemake
        # Convert module configs to old-style format for compatibility
        exported = {}

        # Add system settings
        exported.update(self.system_config)

        # Add module parameters in old format
        for module_name, module_config in self.module_configs.items():
            merged_module_config = self.get_module_config(module_name)

            # Convert nested parameters to old flat format
            for param_group, params in merged_module_config.items():
                for param_name, param_value in params.items():
                    old_key = f"{module_name}_{param_group}_{param_name}"
                    exported[old_key] = param_value

        # Add user settings
        exported.update(self.user_config)

        return exported

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
        # Try user config first
        if module_name in config:
            module_section = config[module_name]
            if isinstance(module_section, dict):
                # Flatten nested parameters
                for param_group, params in module_section.items():
                    if isinstance(params, dict) and param_path in params:
                        return params[param_path]

        # Try module defaults
        if module_name in self.module_configs:
            module_config = self.module_configs[module_name]
            parameters = module_config.get("parameters", {})
            for param_group, params in parameters.items():
                if isinstance(params, dict) and param_path in params:
                    return params[param_path]

        return None

    @classmethod
    def load(cls, user_config_path: str, repo_root: Optional[str] = None) -> "ConfigManager":
        """
        Convenience method to load all configurations.

        Args:
            user_config_path: Path to user config file
            repo_root: Path to repository root (optional)

        Returns:
            ConfigManager instance with loaded configs

        Raises:
            ConfigError: If loading or validation fails
        """
        manager = cls(repo_root)

        # Load all configs
        manager.load_system_config()
        manager.load_module_configs()
        manager.load_user_config(user_config_path)

        # Merge configs
        manager.merge()

        # Validate
        errors = manager.validate()
        if errors:
            raise ValidationError(errors)

        return manager


# Convenience function for quick loading
def load_config(user_config_path: str, repo_root: Optional[str] = None) -> ConfigManager:
    """
    Load and validate all configurations.

    Args:
        user_config_path: Path to user config file
        repo_root: Path to repository root (optional)

    Returns:
        ConfigManager instance with loaded and validated configs

    Raises:
        ConfigError: If loading or validation fails
        ValidationError: If validation fails
    """
    return ConfigManager.load(user_config_path, repo_root)
