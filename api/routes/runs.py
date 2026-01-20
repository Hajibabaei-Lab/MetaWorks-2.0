"""
Run management routes.

This module handles all run-related API endpoints.
"""

from typing import TYPE_CHECKING, Optional

from fastapi import HTTPException
from fastapi.responses import FileResponse

from lib.exceptions import (
    ConfigurationError,
    FileProcessingError,
    StateManagementError,
    ValidationError,
)

from ..schemas import (
    ArtifactResponse,
    DeleteRunResponse,
    LogResponse,
    RunStatus,
    RunSubmissionRequest,
)

if TYPE_CHECKING:
    from fastapi import FastAPI

    from ..job_manager import JobManager


def register_run_routes(app: "FastAPI", manager: "JobManager") -> None:  # noqa: C901
    """Register all run-related routes with the FastAPI app."""

    @app.post("/runs", response_model=RunStatus)
    def submit_run(payload: RunSubmissionRequest) -> RunStatus:
        """Submit a new run for execution."""
        try:
            status = manager.submit_run(payload)
        except (ConfigurationError, ValidationError) as exc:
            raise HTTPException(
                status_code=400,
                detail=f"Configuration error: {exc.message}",
            ) from exc
        except StateManagementError as exc:
            raise HTTPException(
                status_code=409,
                detail=f"Run management error: {exc.message}",
            ) from exc
        except Exception as exc:
            raise HTTPException(
                status_code=500,
                detail=f"Unexpected error: {str(exc)}",
            ) from exc
        return status

    @app.get("/runs/{run_id}", response_model=RunStatus)
    def run_status(run_id: str) -> RunStatus:
        """Get the status of a specific run."""
        try:
            return manager.get_status(run_id)
        except StateManagementError as exc:
            raise HTTPException(
                status_code=404,
                detail=f"Run not found: {exc.message}",
            ) from exc
        except Exception as exc:
            raise HTTPException(
                status_code=500,
                detail=f"Unexpected error: {str(exc)}",
            ) from exc

    @app.get("/runs/{run_id}/logs", response_model=LogResponse)
    def run_logs(run_id: str, lines: Optional[int] = None) -> LogResponse:
        """Get recent log lines for a run."""
        try:
            data = manager.tail_logs(run_id, lines=lines)
            return LogResponse(**data)
        except StateManagementError as exc:
            raise HTTPException(
                status_code=404,
                detail=f"Run not found: {exc.message}",
            ) from exc
        except FileProcessingError as exc:
            raise HTTPException(
                status_code=500,
                detail=f"Log file error: {exc.message}",
            ) from exc
        except Exception as exc:
            raise HTTPException(
                status_code=500,
                detail=f"Unexpected error: {str(exc)}",
            ) from exc

    @app.get("/runs/{run_id}/logs/download")
    def download_log_file(run_id: str) -> FileResponse:
        """Download the full log file for a run."""
        from pathlib import Path

        try:
            status = manager.get_status(run_id)
        except StateManagementError as exc:
            raise HTTPException(
                status_code=404,
                detail=f"Run not found: {exc.message}",
            ) from exc
        log_path = status.log_path
        if not log_path or not Path(log_path).exists():
            raise HTTPException(
                status_code=404,
                detail="Log file not available yet",
            )
        return FileResponse(log_path, filename=f"{run_id}-snakemake.log", media_type="text/plain")

    @app.post("/runs/{run_id}/cancel", response_model=RunStatus)
    def cancel_run(run_id: str) -> RunStatus:
        """Cancel a running run."""
        try:
            return manager.cancel_run(run_id)
        except StateManagementError as exc:
            raise HTTPException(
                status_code=404,
                detail=f"Run not found: {exc.message}",
            ) from exc
        except Exception as exc:
            raise HTTPException(
                status_code=500,
                detail=f"Failed to cancel run: {str(exc)}",
            ) from exc

    @app.get("/runs/{run_id}/artifacts", response_model=ArtifactResponse)
    def package_artifacts(run_id: str) -> ArtifactResponse:
        """Package run artifacts into a tarball."""
        try:
            return manager.package_artifacts(run_id)
        except StateManagementError as exc:
            raise HTTPException(
                status_code=404,
                detail=f"Run not found: {exc.message}",
            ) from exc
        except FileProcessingError as exc:
            raise HTTPException(
                status_code=500,
                detail=f"Artifact packaging error: {exc.message}",
            ) from exc
        except Exception as exc:
            raise HTTPException(
                status_code=500,
                detail=f"Unexpected error: {str(exc)}",
            ) from exc

    @app.get("/runs/{run_id}/artifacts/download")
    def download_artifact(run_id: str) -> FileResponse:
        """Download the artifacts tarball for a run."""
        from pathlib import Path

        try:
            result = manager.package_artifacts(run_id)
        except StateManagementError as exc:
            raise HTTPException(
                status_code=404,
                detail=f"Run not found: {exc.message}",
            ) from exc
        except FileProcessingError as exc:
            raise HTTPException(
                status_code=500,
                detail=f"Artifact packaging error: {exc.message}",
            ) from exc
        except Exception as exc:
            raise HTTPException(
                status_code=500,
                detail=f"Unexpected error: {str(exc)}",
            ) from exc
        archive_path = Path(result.archive_path)
        if not archive_path.exists():
            raise HTTPException(
                status_code=404,
                detail="Artifact archive not found",
            )
        return FileResponse(
            archive_path, filename=f"{run_id}.tar.gz", media_type="application/gzip"
        )

    @app.delete("/runs/{run_id}", response_model=DeleteRunResponse)
    def delete_run(run_id: str) -> DeleteRunResponse:
        """Delete a run and its associated files."""
        try:
            result = manager.delete_run(run_id)
            return DeleteRunResponse(**result)
        except StateManagementError as exc:
            raise HTTPException(
                status_code=404,
                detail=f"Run not found: {exc.message}",
            ) from exc
        except FileProcessingError as exc:
            raise HTTPException(
                status_code=500,
                detail=f"File deletion error: {exc.message}",
            ) from exc
        except Exception as exc:
            raise HTTPException(
                status_code=500,
                detail=f"Unexpected error: {str(exc)}",
            ) from exc
