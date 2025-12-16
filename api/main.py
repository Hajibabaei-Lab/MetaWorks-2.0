import os
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, File, HTTPException, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, PlainTextResponse
from fastapi.staticfiles import StaticFiles

from .config import settings
from .job_manager import JobManager
from .schemas import (
    ArtifactResponse,
    ConfigSectionsResponse,
    DeleteRunResponse,
    ListAssetsResponse,
    LogResponse,
    PathsResponse,
    RenderConfigRequest,
    RenderConfigResponse,
    RunStatus,
    RunSubmissionRequest,
    UploadResponse,
    WorkflowType,
)

app = FastAPI(title="MetaWorks Control Node API", version="0.1.0")
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

manager = JobManager()


@app.get("/health", response_class=PlainTextResponse)
def health() -> str:
    return "ok"


@app.post("/runs", response_model=RunStatus)
def submit_run(payload: RunSubmissionRequest) -> RunStatus:
    try:
        status = manager.submit_run(payload)
    except Exception as exc:  # pylint: disable=broad-except
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    return status


@app.get("/runs/{run_id}", response_model=RunStatus)
def run_status(run_id: str) -> RunStatus:
    try:
        return manager.get_status(run_id)
    except Exception as exc:  # pylint: disable=broad-except
        raise HTTPException(status_code=404, detail=str(exc)) from exc


@app.get("/runs/{run_id}/logs", response_model=LogResponse)
def run_logs(run_id: str, lines: Optional[int] = None) -> LogResponse:
    try:
        data = manager.tail_logs(run_id, lines=lines)
        return LogResponse(**data)
    except Exception as exc:  # pylint: disable=broad-except
        raise HTTPException(status_code=404, detail=str(exc)) from exc


@app.get("/runs/{run_id}/logs/download")
def download_log_file(run_id: str) -> FileResponse:
    try:
        status = manager.get_status(run_id)
    except Exception as exc:  # pylint: disable=broad-except
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    log_path = status.log_path
    if not log_path or not Path(log_path).exists():
        raise HTTPException(status_code=404, detail="Log file not available yet")
    return FileResponse(log_path, filename=f"{run_id}-snakemake.log", media_type="text/plain")


@app.post("/runs/{run_id}/cancel", response_model=RunStatus)
def cancel_run(run_id: str) -> RunStatus:
    try:
        return manager.cancel_run(run_id)
    except Exception as exc:  # pylint: disable=broad-except
        raise HTTPException(status_code=404, detail=str(exc)) from exc


@app.get("/runs/{run_id}/artifacts", response_model=ArtifactResponse)
def package_artifacts(run_id: str) -> ArtifactResponse:
    try:
        return manager.package_artifacts(run_id)
    except Exception as exc:  # pylint: disable=broad-except
        raise HTTPException(status_code=404, detail=str(exc)) from exc


@app.get("/runs/{run_id}/artifacts/download")
def download_artifact(run_id: str) -> FileResponse:
    try:
        result = manager.package_artifacts(run_id)
    except Exception as exc:  # pylint: disable=broad-except
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    archive_path = Path(result.archive_path)
    if not archive_path.exists():
        raise HTTPException(status_code=404, detail="Artifact archive not found")
    return FileResponse(archive_path, filename=f"{run_id}.tar.gz", media_type="application/gzip")


@app.get("/configs/defaults/{workflow}", response_class=PlainTextResponse)
def get_default_config(workflow: WorkflowType) -> str:
    config_path = settings.default_configs.get(workflow.value)
    if not config_path or not Path(config_path).exists():
        raise HTTPException(status_code=404, detail=f"Default config not found for {workflow.value}")
    return Path(config_path).read_text()


@app.get("/configs/defaults/{workflow}/sections", response_model=ConfigSectionsResponse)
def get_default_config_sections(workflow: WorkflowType) -> ConfigSectionsResponse:
    try:
        sections = manager.default_config_sections(workflow)
    except Exception as exc:  # pylint: disable=broad-except
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    return ConfigSectionsResponse(workflow=workflow, sections=sections)


@app.delete("/runs/{run_id}", response_model=DeleteRunResponse)
def delete_run(run_id: str) -> DeleteRunResponse:
    try:
        result = manager.delete_run(run_id)
        return DeleteRunResponse(**result)
    except Exception as exc:  # pylint: disable=broad-except
        raise HTTPException(status_code=404, detail=str(exc)) from exc


@app.post("/configs/render", response_model=RenderConfigResponse)
def render_config(payload: RenderConfigRequest) -> RenderConfigResponse:
    try:
        return manager.render_config(payload.workflow, payload.overrides)
    except Exception as exc:  # pylint: disable=broad-except
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.get("/settings/paths", response_model=PathsResponse)
def settings_paths() -> PathsResponse:
    return PathsResponse(repo_root=str(settings.repo_root), runtime_cache=str(settings.singularity_cache))


def _save_upload(target_root: Path, file: UploadFile) -> UploadResponse:
    target_root.mkdir(parents=True, exist_ok=True)
    target_path = target_root / file.filename
    content = file.file.read()
    target_path.write_bytes(content)
    return UploadResponse(name=file.filename, path=str(target_path))


@app.post("/classifiers", response_model=UploadResponse)
def upload_classifier(file: UploadFile = File(...)) -> UploadResponse:
    try:
        return _save_upload(settings.classifier_root, file)
    except Exception as exc:  # pylint: disable=broad-except
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.get("/classifiers", response_model=ListAssetsResponse)
def list_classifiers() -> ListAssetsResponse:
    items = [p.name for p in settings.classifier_root.glob("*") if p.is_file()]
    return ListAssetsResponse(items=items)


@app.post("/adapters", response_model=UploadResponse)
def upload_adapter(file: UploadFile = File(...)) -> UploadResponse:
    try:
        return _save_upload(settings.adapter_root, file)
    except Exception as exc:  # pylint: disable=broad-except
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.get("/adapters", response_model=ListAssetsResponse)
def list_adapters() -> ListAssetsResponse:
    items = [p.name for p in settings.adapter_root.glob("*") if p.is_file()]
    return ListAssetsResponse(items=items)


@app.post("/delete/classifiers/{name}", response_model=ListAssetsResponse)
def delete_classifier(name: str) -> ListAssetsResponse:
    target = settings.classifier_root / name
    if target.exists():
        target.unlink()
    return list_classifiers()


@app.post("/delete/adapters/{name}", response_model=ListAssetsResponse)
def delete_adapter(name: str) -> ListAssetsResponse:
    target = settings.adapter_root / name
    if target.exists():
        target.unlink()
    return list_adapters()


ui_path = Path(__file__).resolve().parent.parent / "ui"
if ui_path.exists():
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
