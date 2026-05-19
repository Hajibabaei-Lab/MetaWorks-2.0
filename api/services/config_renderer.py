"""Config rendering and writing for pipeline runs."""

import logging
from pathlib import Path
from typing import Any, Dict

import yaml

from lib.config import ConfigManager
from lib.exceptions import ConfigurationError

from ..schemas import RenderConfigResponse, WorkflowType

logger = logging.getLogger(__name__)


class ConfigRenderer:
    """Renders pipeline configurations from profiles and overrides."""

    def __init__(self, repo_root: str):
        self._repo_root = repo_root

    def _get_config_manager(self) -> ConfigManager:
        return ConfigManager(repo_root=self._repo_root)

    def render_config(
        self,
        profile: str,
        workflow: WorkflowType,
        overrides: Dict[str, Any],
    ) -> RenderConfigResponse:
        """Render configuration using profile, workflow, and user overrides."""
        try:
            config_manager = self._get_config_manager()
            config_manager.load_defaults_config()
            config_manager.load_module_configs()
            config_manager.load_profile(profile)

            if overrides:
                config_manager.user_config = dict(overrides)

            resolved = config_manager.merge(workflow=workflow.value)
            merged_dict = resolved.export_for_workflow()
            rendered = yaml.safe_dump(merged_dict, sort_keys=False)

            return RenderConfigResponse(workflow=workflow, merged=rendered)
        except Exception as exc:
            raise ConfigurationError(
                f"Failed to render config: {str(exc)}",
                config_key="render",
                suggestion="Check profile name and overrides",
            ) from exc

    def default_config_sections(self, workflow: WorkflowType) -> Dict[str, Any]:
        """Get default config sections (legacy method)."""
        config_manager = self._get_config_manager()
        config_manager.load_defaults_config()
        config_manager.load_module_configs()
        resolved = config_manager.merge(workflow=workflow.value)
        return resolved.export_for_workflow()

    def write_rendered_config(
        self,
        run_dir: Path,
        profile: str,
        workflow: WorkflowType,
        overrides: Dict[str, Any],
    ) -> Path:
        """Write rendered config to run directory. Returns config file path."""
        rendered = self.render_config(profile, workflow, overrides)
        config_path = run_dir / "rendered_config.yaml"
        config_path.write_text(rendered.merged)
        return config_path
