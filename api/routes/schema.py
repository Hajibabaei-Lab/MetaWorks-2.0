"""Config schema discovery route — returns Galaxy-style JSON schema for the Vue form builder."""

from typing import TYPE_CHECKING

from fastapi import APIRouter, HTTPException, Query

from lib.config import ConfigManager
from lib.config.loader import ConfigError
from lib.config.schema_builder import build_config_schema

from ..schemas import ConfigSchemaResponse

if TYPE_CHECKING:
    from ..config import Settings
    from ..job_manager import JobManager


def _build_schema_router(manager: "JobManager", settings: "Settings") -> APIRouter:
    router = APIRouter(tags=["config"])

    def _get_config_manager() -> ConfigManager:
        return ConfigManager(repo_root=str(settings.repo_root))

    @router.get("/config/schema", response_model=ConfigSchemaResponse)
    def get_config_schema(
        profile: str = Query(default="coi", description="Profile name (e.g. coi, 16s, its)"),
        workflow: str = Query(default="esv", description="Workflow type (esv or otu)"),
    ) -> ConfigSchemaResponse:
        try:
            config_manager = _get_config_manager()
            defaults = config_manager.load_defaults_config()
            profile_config = config_manager.load_profile(profile)
            schema_dict = build_config_schema(
                defaults_config=defaults,
                profile_config=profile_config,
                profile_name=profile,
                workflow=workflow,
            )
            return ConfigSchemaResponse(**schema_dict)
        except ConfigError as exc:
            raise HTTPException(
                status_code=404,
                detail=str(exc),
            ) from exc
        except Exception as exc:
            raise HTTPException(
                status_code=500,
                detail=f"Failed to build schema: {exc}",
            ) from exc

    return router


def register_schema_routes(
    app,
    manager: "JobManager",
    settings: "Settings",
    *,
    prefix: str = "",
    include_in_schema: bool = True,
) -> None:
    app.include_router(
        _build_schema_router(manager, settings),
        prefix=prefix,
        include_in_schema=include_in_schema,
    )
