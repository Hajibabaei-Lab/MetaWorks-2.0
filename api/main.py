"""
MetaWorks Control Node API.

Main FastAPI application for managing bioinformatics pipeline runs.
"""

import os

from fastapi import APIRouter, FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, PlainTextResponse, RedirectResponse
from fastapi.staticfiles import StaticFiles
from pydantic import ValidationError as PydanticValidationError

from .config import Settings, settings
from .job_manager import JobManager
from .routes import assets, configs, runs, schema

API_PREFIX = "/api"


def _parse_cors_origins(raw_origins: str) -> list[str]:
    origins = [item.strip() for item in str(raw_origins or "").split(",") if item.strip()]
    return origins or ["*"]


def _build_health_router() -> APIRouter:
    router = APIRouter(tags=["system"])

    @router.get("/health", response_class=PlainTextResponse)
    def health() -> str:
        """Health check endpoint."""
        return "ok"

    return router


def _mount_legacy_ui(app: FastAPI, app_settings: Settings) -> None:
    ui_dist = app_settings.repo_root / "ui" / "dist"
    ui_path = app_settings.repo_root / "ui"

    if ui_dist.exists():
        app.mount("/", StaticFiles(directory=ui_dist, html=True), name="ui")
    elif ui_path.exists():
        app.mount("/", StaticFiles(directory=ui_path, html=True), name="ui")


def create_app(app_settings: Settings | None = None) -> FastAPI:
    """Create the FastAPI application."""
    active_settings = app_settings or settings
    manager = JobManager()

    app = FastAPI(
        title="MetaWorks Control Node API",
        version="0.2.0",
        docs_url=f"{API_PREFIX}/docs",
        openapi_url=f"{API_PREFIX}/openapi.json",
        redoc_url=None,
    )

    origins = _parse_cors_origins(active_settings.cors_allowed_origins)
    allow_credentials = origins != ["*"]
    app.add_middleware(
        CORSMiddleware,
        allow_origins=origins,
        allow_credentials=allow_credentials,
        allow_methods=["*"],
        allow_headers=["*"],
    )

    @app.exception_handler(PydanticValidationError)
    async def validation_exception_handler(request, exc):
        """Handle Pydantic validation errors with detailed error messages."""
        errors = exc.errors()
        error_messages = []

        for error in errors:
            field = " -> ".join(str(loc) for loc in error["loc"])
            message = error["msg"]
            error_messages.append(f"{field}: {message}")

        return JSONResponse(
            status_code=400,
            content={"detail": f"Validation error: {'; '.join(error_messages)}"},
            headers={"X-Error": "Validation failed"},
        )

    health_router = _build_health_router()
    app.include_router(health_router, prefix=API_PREFIX)

    runs.register_run_routes(app, manager, prefix=API_PREFIX)
    configs.register_config_routes(app, manager, active_settings, prefix=API_PREFIX)
    assets.register_asset_routes(app, active_settings, prefix=API_PREFIX)
    schema.register_schema_routes(app, manager, active_settings, prefix=API_PREFIX)

    @app.get("/docs", include_in_schema=False)
    def docs_redirect() -> RedirectResponse:
        """Legacy docs alias."""
        return RedirectResponse(url=f"{API_PREFIX}/docs")

    @app.get("/openapi.json", include_in_schema=False)
    def openapi_redirect() -> RedirectResponse:
        """Legacy OpenAPI alias."""
        return RedirectResponse(url=f"{API_PREFIX}/openapi.json")

    if active_settings.serve_legacy_ui:
        _mount_legacy_ui(app, active_settings)

    for path in [
        active_settings.run_root,
        active_settings.artifact_root,
        active_settings.adapter_root,
        active_settings.classifier_root,
        active_settings.staging_root,
    ]:
        os.makedirs(path, exist_ok=True)

    return app


app = create_app()
