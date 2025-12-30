"""
Configuration management routes.

This module handles all configuration-related API endpoints.
"""

from pathlib import Path

from fastapi import HTTPException
from fastapi.responses import PlainTextResponse

from lib.exceptions import ConfigurationError, ValidationError

from ..schemas import (
    ConfigSectionsResponse,
    RenderConfigRequest,
    RenderConfigResponse,
    WorkflowType,
)


def register_config_routes(app, manager, settings):
    """Register all configuration-related routes with the FastAPI app."""

    @app.get("/configs/defaults/{workflow}", response_class=PlainTextResponse)
    def get_default_config(workflow: WorkflowType) -> str:
        """Get default configuration for a workflow."""
        config_path = settings.default_configs.get(workflow.value)
        if not config_path or not Path(config_path).exists():
            raise HTTPException(
                status_code=404,
                detail=f"Default config not found for {workflow.value}",
            )
        try:
            return Path(config_path).read_text()
        except IOError as exc:
            raise HTTPException(
                status_code=500,
                detail=f"Failed to read config file: {str(exc)}",
            ) from exc

    @app.get("/configs/defaults/{workflow}/sections", response_model=ConfigSectionsResponse)
    def get_default_config_sections(workflow: WorkflowType) -> ConfigSectionsResponse:
        """Get configuration sections for a workflow."""
        try:
            sections = manager.default_config_sections(workflow)
        except ConfigurationError as exc:
            raise HTTPException(
                status_code=404,
                detail=f"Configuration error: {exc.message}",
            ) from exc
        except Exception as exc:
            raise HTTPException(
                status_code=500,
                detail=f"Unexpected error: {str(exc)}",
            ) from exc
        return ConfigSectionsResponse(workflow=workflow, sections=sections)

    @app.post("/configs/render", response_model=RenderConfigResponse)
    def render_config(payload: RenderConfigRequest) -> RenderConfigResponse:
        """Render configuration with overrides."""
        try:
            return manager.render_config(payload.workflow, payload.overrides)
        except ConfigurationError as exc:
            raise HTTPException(
                status_code=400,
                detail=f"Configuration error: {exc.message}",
            ) from exc
        except ValidationError as exc:
            raise HTTPException(
                status_code=400,
                detail=f"Validation error: {exc.message}",
            ) from exc
        except Exception as exc:
            raise HTTPException(
                status_code=500,
                detail=f"Unexpected error: {str(exc)}",
            ) from exc

    @app.get("/settings/paths")
    def settings_paths():
        """Get system paths."""
        from ..schemas import PathsResponse

        return PathsResponse(
            repo_root=str(settings.repo_root), runtime_cache=str(settings.singularity_cache)
        )
