from pathlib import Path
from typing import Dict, List

from pydantic import Field

try:  # Pydantic v2+
    from pydantic_settings import BaseSettings
except Exception:  # pragma: no cover
    from pydantic import BaseSettings


# Resolve the repository root relative to this file so defaults work on any host path.
DEFAULT_REPO_ROOT = Path(__file__).resolve().parent.parent


class Settings(BaseSettings):
    """Runtime settings for the API and runtime-aware Snakemake runner."""

    repo_root: Path = Field(
        default=DEFAULT_REPO_ROOT,
        description="Root of the MetaWorks checkout used for Snakemake execution.",
    )
    run_root: Path = Field(
        default=DEFAULT_REPO_ROOT / "runtime" / "runs",
        description="Root directory where per-run working directories are created.",
    )
    staging_root: Path = Field(
        default=DEFAULT_REPO_ROOT / "runtime" / "staging",
        description="Location for staging uploads before submission.",
    )
    artifact_root: Path = Field(
        default=DEFAULT_REPO_ROOT / "runtime" / "artifacts",
        description="Location where packaged archives are written.",
    )
    classifier_root: Path = Field(
        default=DEFAULT_REPO_ROOT / "runtime" / "classifiers",
        description="Managed classifier files live here.",
    )
    adapter_root: Path = Field(
        default=DEFAULT_REPO_ROOT / "runtime" / "adapters",
        description="Managed adapter FASTA files live here.",
    )
    state_file: Path = Field(
        default=DEFAULT_REPO_ROOT / "runtime" / "state" / "runs.json",
        description="Tracks run metadata, PIDs/return codes, and paths.",
    )
    snakefiles: Dict[str, Path] = Field(
        default={
            "esv": Path("workflow/Snakefile"),
        },
        description="Workflow entrypoints keyed by workflow name.",
    )
    scheduler: str = Field(
        default="slurm",
        description="Scheduler flavor (slurm|lsf|pbs).",
    )
    submit_enabled: bool = Field(
        default=False,
        description="If False, scheduler submit is dry-run and only records the command.",
    )
    default_runtime: str = Field(
        default="docker", description="Default runtime if not provided (conda|docker|apptainer)."
    )
    allowed_runtimes: str = Field(
        default="docker,apptainer",
        description="Comma-separated runtimes allowed for web submissions.",
    )
    retention_policy: str = Field(
        default="until_download",
        description="Run data retention mode: until_download|immediate|manual.",
    )
    container_image: str = Field(
        default="docker://metaworks:latest",
        description="Default container URI or path for docker/apptainer runs.",
    )
    bind_paths: List[str] = Field(
        default_factory=list,
        description="Optional additional bind mounts (src:dest) appended to auto binds for container runs.",
    )
    singularity_cache: Path = Field(
        default=DEFAULT_REPO_ROOT / "runtime" / "cache",
        description="Default cache/prefix for Apptainer/Singularity pulls.",
    )
    default_cores: int = Field(default=32, description="Cores requested per job.")
    default_mem_gb: int = Field(default=120, description="Memory requested per job.")
    log_tail_lines: int = Field(default=200, description="Default lines to return for log tail.")
    serve_legacy_ui: bool = Field(
        default=True,
        description="Serve the legacy static UI from the backend process.",
    )
    cors_allowed_origins: str = Field(
        default="http://localhost:5173,http://127.0.0.1:5173,http://localhost:8080,http://127.0.0.1:8080",
        description="Comma-separated origins allowed to access the API directly during development.",
    )

    class Config:
        env_prefix = "METAWORKS_"


settings = Settings()
