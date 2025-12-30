from __future__ import annotations

from datetime import datetime
from enum import Enum
from typing import Any, Dict, List, Optional

from pydantic import BaseModel, Field, field_validator, model_validator


class WorkflowType(str, Enum):
    esv = "esv"
    otu = "otu"


class RuntimeType(str, Enum):
    conda = "conda"
    docker = "docker"
    apptainer = "apptainer"


class RunSubmissionRequest(BaseModel):
    workflow: WorkflowType = Field(..., description="Workflow to run (esv|otu).")
    run_name: str = Field(
        ..., min_length=1, max_length=255, description="Friendly name for the run."
    )
    runtime: RuntimeType = Field(
        default=RuntimeType.conda, description="Runtime executor (conda|docker|apptainer)."
    )
    container_image: Optional[str] = Field(
        default=None,
        min_length=1,
        description="Container URI or path (docker://image or .sif) for docker/apptainer runtimes.",
    )
    bind_paths: List[str] = Field(
        default_factory=list,
        description="Extra bind mounts (src:dest) when using docker/apptainer, run/input paths are auto-bound.",
    )
    cache_dir: Optional[str] = Field(
        default=None,
        description="Container cache directory for docker/apptainer runtimes.",
    )
    config_overrides: Dict[str, Any] = Field(
        default_factory=dict, description="Partial config overrides to merge with defaults."
    )
    input_dir: str = Field(
        ..., min_length=1, description="Path to input FASTQ folder on control node or shared FS."
    )
    samples_csv: Optional[str] = Field(
        default=None,
        description="Optional samples CSV path if using sample_source=csv.",
    )
    sample_source: str = Field(
        default="folder",
        pattern="^(folder|csv)$",
        description="folder|csv (matches Snakemake config).",
    )
    notes: Optional[str] = Field(
        default=None, max_length=2000, description="Optional free-form notes for the run."
    )
    cores: Optional[int] = Field(
        default=None, ge=1, le=256, description="Override default cores for scheduler job (1-256)."
    )
    mem_gb: Optional[int] = Field(
        default=None,
        ge=1,
        le=1024,
        description="Override default memory (GB) for scheduler job (1-1024 GB).",
    )
    dry_run: bool = Field(
        default=False, description="If true, do not submit to scheduler—only stage config."
    )

    @field_validator("input_dir")
    @classmethod
    def validate_input_dir(cls, v):
        """Validate input_dir is a valid path."""
        if not v or not v.strip():
            raise ValueError("input_dir cannot be empty")
        return v

    @model_validator(mode="after")
    @classmethod
    def validate_container_image(cls, data):
        """Ensure container_image is provided for container runtimes."""
        if data.runtime in (RuntimeType.docker, RuntimeType.apptainer) and not data.container_image:
            raise ValueError(f"container_image is required when using {data.runtime} runtime")
        return data


class RunStatus(BaseModel):
    run_id: str
    workflow: WorkflowType
    runtime: RuntimeType
    status: str
    scheduler_job_id: Optional[str] = None
    submitted_at: Optional[datetime] = None
    started_at: Optional[datetime] = None
    ended_at: Optional[datetime] = None
    pid: Optional[int] = None
    return_code: Optional[int] = None
    message: Optional[str] = None
    run_dir: Optional[str] = None
    log_path: Optional[str] = None
    command: Optional[str] = None
    artifact_path: Optional[str] = None


class LogResponse(BaseModel):
    run_id: str
    tail: List[str]


class RenderConfigRequest(BaseModel):
    workflow: WorkflowType
    overrides: Dict[str, Any] = Field(default_factory=dict)


class RenderConfigResponse(BaseModel):
    workflow: WorkflowType
    merged: str


class ConfigSectionsResponse(BaseModel):
    workflow: WorkflowType
    sections: Dict[str, Any]


class UploadResponse(BaseModel):
    name: str
    path: str


class ListAssetsResponse(BaseModel):
    items: List[str]


class ArtifactResponse(BaseModel):
    run_id: str
    archive_path: str


class DeleteRunResponse(BaseModel):
    run_id: str
    removed_paths: List[str]


class PathsResponse(BaseModel):
    repo_root: str
    runtime_cache: str
