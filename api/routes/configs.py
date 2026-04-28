"""Configuration management routes."""

from typing import TYPE_CHECKING

import yaml
from fastapi import APIRouter, HTTPException
from fastapi.responses import PlainTextResponse

from lib.config import ConfigManager
from lib.exceptions import ConfigurationError, ValidationError

from ..schemas import (
    ConfigSectionsResponse,
    ProfileInfo,
    ProfileListResponse,
    RenderConfigResponse,
    RenderConfigWithProfileRequest,
    WorkflowType,
)

if TYPE_CHECKING:
    from ..config import Settings
    from ..job_manager import JobManager


def _build_config_router(manager: "JobManager", settings: "Settings") -> APIRouter:
    """Build the configuration router."""

    router = APIRouter(tags=["config"])

    def _get_config_manager() -> ConfigManager:
        """Get a ConfigManager instance with repo root from settings."""
        return ConfigManager(repo_root=str(settings.repo_root))

    @router.get("/profiles", response_model=ProfileListResponse)
    def list_profiles() -> ProfileListResponse:
        """List all available configuration profiles."""
        try:
            config_manager = _get_config_manager()
            profiles = config_manager.list_available_profiles()
            return ProfileListResponse(
                profiles=[ProfileInfo(**p) for p in profiles]
            )
        except Exception as exc:
            raise HTTPException(
                status_code=500,
                detail=f"Failed to list profiles: {str(exc)}",
            ) from exc

    @router.get("/profiles/{profile_name}")
    def get_profile(profile_name: str) -> dict:
        """Get a specific profile's configuration."""
        try:
            config_manager = _get_config_manager()
            config_manager.load_profile(profile_name)
            return config_manager.profile_config or {}
        except ConfigurationError as exc:
            raise HTTPException(
                status_code=404,
                detail=exc.message,
            ) from exc
        except Exception as exc:
            raise HTTPException(
                status_code=500,
                detail=f"Failed to load profile: {str(exc)}",
            ) from exc

    @router.get("/configs/defaults")
    def get_defaults() -> dict:
        """Get the pipeline defaults configuration."""
        try:
            config_manager = _get_config_manager()
            config_manager.load_defaults_config()
            return config_manager.defaults_config or {}
        except Exception as exc:
            raise HTTPException(
                status_code=500,
                detail=f"Failed to load defaults: {str(exc)}",
            ) from exc

    @router.post("/configs/render", response_model=RenderConfigResponse)
    def render_config(payload: RenderConfigWithProfileRequest) -> RenderConfigResponse:
        """Render configuration with profile and overrides."""
        try:
            config_manager = _get_config_manager()
            
            # Load defaults and profile
            config_manager.load_defaults_config()
            config_manager.load_module_configs()
            config_manager.load_profile(payload.profile)
            
            # Apply overrides
            if payload.overrides:
                config_manager.user_config = dict(payload.overrides)
            
            # Merge with workflow type
            config_manager.merge(workflow=payload.workflow.value)
            
            # Export as YAML
            merged_dict = config_manager.export_for_workflow(payload.workflow.value)
            merged_yaml = yaml.safe_dump(merged_dict, sort_keys=False)
            
            return RenderConfigResponse(workflow=payload.workflow, merged=merged_yaml)
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

    @router.get("/configs/defaults/{workflow}", response_class=PlainTextResponse)
    def get_default_config(workflow: WorkflowType) -> str:
        """Get default configuration for a workflow (legacy endpoint)."""
        # Load defaults and render with no profile
        try:
            config_manager = _get_config_manager()
            config_manager.load_defaults_config()
            config_manager.load_module_configs()
            config_manager.merge(workflow=workflow.value)
            merged_dict = config_manager.export_for_workflow(workflow.value)
            return yaml.safe_dump(merged_dict, sort_keys=False)
        except Exception as exc:
            raise HTTPException(
                status_code=500,
                detail=f"Failed to load defaults: {str(exc)}",
            ) from exc

    @router.get("/configs/defaults/{workflow}/sections", response_model=ConfigSectionsResponse)
    def get_default_config_sections(workflow: WorkflowType) -> ConfigSectionsResponse:
        """Get configuration sections for a workflow."""
        try:
            config_manager = _get_config_manager()
            config_manager.load_defaults_config()
            config_manager.load_module_configs()
            config_manager.merge(workflow=workflow.value)
            sections = config_manager.export_for_workflow(workflow.value)
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

    @router.get("/settings/paths")
    def settings_paths():
        """Get system paths."""
        from ..schemas import PathsResponse

        allowed_runtimes = [
            item.strip().lower()
            for item in str(settings.allowed_runtimes or "").split(",")
            if item.strip()
        ]
        if not allowed_runtimes:
            allowed_runtimes = ["docker", "apptainer"]

        return PathsResponse(
            repo_root=str(settings.repo_root),
            runtime_cache=str(settings.singularity_cache),
            allowed_runtimes=allowed_runtimes,
            retention_policy=str(settings.retention_policy or "until_download").lower(),
            default_runtime=str(settings.default_runtime or "docker").lower(),
            container_image=str(settings.container_image or "docker://metaworks:latest"),
        )

    return router


def register_config_routes(
    app,
    manager: "JobManager",
    settings: "Settings",
    *,
    prefix: str = "",
    include_in_schema: bool = True,
) -> None:
    """Register all configuration-related routes with the FastAPI app."""
    app.include_router(
        _build_config_router(manager, settings),
        prefix=prefix,
        include_in_schema=include_in_schema,
    )
