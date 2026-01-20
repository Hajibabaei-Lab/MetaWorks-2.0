"""
MetaWorks Control Node API.

Main FastAPI application for managing bioinformatics pipeline runs.
"""

import os
from pathlib import Path

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import PlainTextResponse
from fastapi.staticfiles import StaticFiles
from pydantic import ValidationError as PydanticValidationError

from .config import settings
from .job_manager import JobManager
from .routes import assets, configs, runs

# Create FastAPI app
app = FastAPI(title="MetaWorks Control Node API", version="0.1.0")

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Initialize job manager
manager = JobManager()


@app.exception_handler(PydanticValidationError)
async def validation_exception_handler(request, exc):
    """Handle Pydantic validation errors with detailed error messages."""
    errors = exc.errors()
    error_messages = []

    for error in errors:
        field = " -> ".join(str(loc) for loc in error["loc"])
        message = error["msg"]
        error_messages.append(f"{field}: {message}")

    error_detail = "; ".join(error_messages)

    from fastapi import HTTPException

    return HTTPException(
        status_code=400,
        detail=f"Validation error: {error_detail}",
        headers={"X-Error": "Validation failed"},
    )


@app.get("/health", response_class=PlainTextResponse)
def health() -> str:
    """Health check endpoint."""
    return "ok"


# Register route modules
runs.register_run_routes(app, manager)
configs.register_config_routes(app, manager, settings)
assets.register_asset_routes(app, settings)


# Mount UI static files
# Check for Vue build output first, then fall back to old UI
ui_dist = Path(__file__).resolve().parent.parent / "ui" / "dist"
ui_path = Path(__file__).resolve().parent.parent / "ui"

if ui_dist.exists():
    # Serve Vue build output (modern UI)
    app.mount("/", StaticFiles(directory=ui_dist, html=True), name="ui")
elif ui_path.exists():
    # Fall back to old UI if build output doesn't exist
    app.mount("/", StaticFiles(directory=ui_path, html=True), name="ui")


# Ensure runtime directories exist on startup
for path in [
    settings.run_root,
    settings.artifact_root,
    settings.adapter_root,
    settings.classifier_root,
    settings.staging_root,
]:
    os.makedirs(path, exist_ok=True)
