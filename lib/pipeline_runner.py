"""
Pipeline runner abstractions for MetaWorks workflows.

This module provides a light-weight orchestration layer that can be reused by
the forthcoming API as well as CLI utilities.  It does not deal with any
specific web framework concerns and is intentionally self-contained.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import threading
import uuid
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Union

from lib.exceptions import FileProcessingError, PipelineExecutionError
from lib.utils.file_utils import ensure_directory


class WorkflowType(str, Enum):
    """Supported workflows that can be triggered by the runner."""

    ESV = "esv"
    OTU = "otu"


class RunStatus(str, Enum):
    """Lifecycle states for a workflow execution."""

    QUEUED = "queued"
    STAGING = "staging"
    RUNNING = "running"
    CANCELLED = "cancelled"
    FAILED = "failed"
    COMPLETED = "completed"


@dataclass
class PipelineRunRequest:
    """
    Description of a run submission coming from an API or CLI client.

    Attributes
    ----------
    workflow:
        Workflow type to execute (ESV or OTU).
    config:
        Path to the base configuration that will be copied into the run space.
    config_overrides:
        Optional dictionary of overrides that should be merged into the config.
        These are stored alongside the run for reproducibility.
    input_paths:
        Paths to input files that need to be staged into the working directory.
    classifiers:
        Optional paths to classifier artifacts referenced by the config.
    adapters:
        Optional paths to adapter fasta files required by the config.
    """

    workflow: WorkflowType
    config: Path
    config_overrides: Optional[Dict[str, Any]] = None
    input_paths: Optional[Iterable[Path]] = None
    classifiers: Optional[Iterable[Path]] = None
    adapters: Optional[Iterable[Path]] = None


@dataclass
class PipelineRunInfo:
    """Metadata for an individual run, persisted on disk as JSON."""

    run_id: str
    workflow: WorkflowType
    status: RunStatus
    workdir: Path
    created_at: datetime
    started_at: Optional[datetime] = None
    finished_at: Optional[datetime] = None
    message: Optional[str] = None
    log_path: Optional[Path] = None
    exit_code: Optional[int] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_json(self) -> Dict[str, Any]:
        """Serialize the run info for storage."""
        return {
            "run_id": self.run_id,
            "workflow": self.workflow.value,
            "status": self.status.value,
            "workdir": str(self.workdir),
            "created_at": self.created_at.isoformat(),
            "started_at": self.started_at.isoformat() if self.started_at else None,
            "finished_at": self.finished_at.isoformat() if self.finished_at else None,
            "message": self.message,
            "log_path": str(self.log_path) if self.log_path else None,
            "exit_code": self.exit_code,
            "metadata": self.metadata,
        }

    @classmethod
    def from_json(cls, payload: Dict[str, Any]) -> "PipelineRunInfo":
        """Deserialize a run record that was saved to disk."""

        def _parse_time(value: Optional[str]) -> Optional[datetime]:
            if not value:
                return None
            return datetime.fromisoformat(value)

        return cls(
            run_id=payload["run_id"],
            workflow=WorkflowType(payload["workflow"]),
            status=RunStatus(payload["status"]),
            workdir=Path(payload["workdir"]),
            created_at=_parse_time(payload.get("created_at")) or datetime.utcnow(),
            started_at=_parse_time(payload.get("started_at")),
            finished_at=_parse_time(payload.get("finished_at")),
            message=payload.get("message"),
            log_path=Path(payload["log_path"]) if payload.get("log_path") else None,
            exit_code=payload.get("exit_code"),
            metadata=payload.get("metadata") or {},
        )


class PipelineRunner:
    """
    Coordinates staging and execution of MetaWorks workflows.

    This implementation focuses on filesystem-based coordination so it can run
    locally or be adapted to HPC schedulers later without changing API code.
    """

    METADATA_FILENAME = "run_metadata.json"

    def __init__(
        self,
        runs_root: Union[str, Path],
        snakefiles: Optional[Dict[WorkflowType, Union[str, Path]]] = None,
        env: Optional[Dict[str, str]] = None,
    ) -> None:
        """
        Parameters
        ----------
        runs_root:
            Directory that will contain per-run working directories.
        snakefiles:
            Mapping from workflow type to snakefile path. Defaults to
            Snakefile_ESV for ESV and workflows/otu_clustering.smk for OTU.
        env:
            Optional environment variables that should be passed to subprocesses.
        """
        self.runs_root = ensure_directory(runs_root)
        self.env = env
        default_snakefiles = {
            WorkflowType.ESV: Path("Snakefile_ESV"),
            WorkflowType.OTU: Path("workflows/otu_clustering.smk"),
        }
        if snakefiles:
            for key, value in snakefiles.items():
                default_snakefiles[key] = Path(value)
        self.snakefiles = default_snakefiles

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------
    def submit_run(self, request: PipelineRunRequest) -> PipelineRunInfo:
        """
        Stage inputs and launch a workflow asynchronously.

        Returns metadata immediately; callers can poll `get_run`.
        """
        run_id = uuid.uuid4().hex
        run_dir = ensure_directory(self.runs_root / run_id)
        staged_config = self._stage_config(run_dir, request)
        staged_inputs = self._stage_inputs(run_dir, request.input_paths)
        metadata = {
            "config": str(staged_config),
            "inputs": [str(path) for path in staged_inputs],
            "workflow": request.workflow.value,
        }
        run_info = PipelineRunInfo(
            run_id=run_id,
            workflow=request.workflow,
            status=RunStatus.QUEUED,
            workdir=run_dir,
            created_at=datetime.utcnow(),
            metadata=metadata,
        )
        self._persist_metadata(run_info)

        thread = threading.Thread(
            target=self._execute_run,
            args=(run_info, staged_config),
            kwargs={"extra_env": self.env},
            daemon=True,
        )
        thread.start()

        return run_info

    def get_run(self, run_id: str) -> Optional[PipelineRunInfo]:
        """Return stored run metadata if available."""
        metadata_path = self._metadata_path(self.runs_root / run_id)
        if not metadata_path.exists():
            return None
        payload = json.loads(metadata_path.read_text())
        return PipelineRunInfo.from_json(payload)

    def list_runs(self) -> Dict[str, PipelineRunInfo]:
        """Return all runs keyed by run_id."""
        runs: Dict[str, PipelineRunInfo] = {}
        for run_dir in self.runs_root.iterdir():
            metadata_path = self._metadata_path(run_dir)
            if metadata_path.exists():
                payload = json.loads(metadata_path.read_text())
                info = PipelineRunInfo.from_json(payload)
                runs[info.run_id] = info
        return runs

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    def _execute_run(
        self,
        run_info: PipelineRunInfo,
        config_path: Path,
        extra_env: Optional[Dict[str, str]] = None,
    ) -> None:
        """Launch Snakemake and update metadata as the job progresses."""
        run_info.status = RunStatus.STAGING
        self._persist_metadata(run_info)

        snakefile = self.snakefiles[run_info.workflow]
        log_path = run_info.workdir / "snakemake.log"
        run_info.log_path = log_path
        run_info.status = RunStatus.RUNNING
        run_info.started_at = datetime.utcnow()
        self._persist_metadata(run_info)

        cmd = [
            "snakemake",
            "--snakefile",
            str(snakefile),
            "--configfile",
            str(config_path),
            "--directory",
            str(run_info.workdir),
        ]

        env = dict(extra_env or {})
        try:
            with log_path.open("w") as log_file:
                process = subprocess.run(
                    cmd,
                    cwd=run_info.workdir,
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                    check=False,
                    env=env if env else None,
                )
        except FileNotFoundError as exc:
            run_info.status = RunStatus.FAILED
            run_info.message = "Snakemake executable not found"
            run_info.finished_at = datetime.utcnow()
            run_info.exit_code = -1
            self._persist_metadata(run_info)
            raise PipelineExecutionError(
                f"Snakemake executable not found: {exc}",
                stage="execution",
                suggestion="Install Snakemake: conda install -c bioconda snakemake",
            ) from exc

        run_info.exit_code = process.returncode
        run_info.finished_at = datetime.utcnow()
        if process.returncode == 0:
            run_info.status = RunStatus.COMPLETED
        else:
            run_info.status = RunStatus.FAILED
            run_info.message = f"Snakemake exited with code {process.returncode}"
        self._persist_metadata(run_info)

    def _stage_config(self, run_dir: Path, request: PipelineRunRequest) -> Path:
        """Copy base config into the run directory."""
        try:
            config_destination = run_dir / request.config.name
            shutil.copy2(request.config, config_destination)
        except (FileNotFoundError, IOError) as exc:
            raise FileProcessingError(
                f"Failed to copy config file: {request.config}",
                filepath=str(request.config),
                suggestion="Check that config file exists and is readable",
            ) from exc

        if request.config_overrides:
            try:
                overrides_path = run_dir / "config_overrides.json"
                overrides_path.write_text(json.dumps(request.config_overrides, indent=2))
            except IOError as exc:
                raise FileProcessingError(
                    f"Failed to write config overrides: {overrides_path}",
                    filepath=str(overrides_path),
                    suggestion="Check write permissions in run directory",
                ) from exc

        return config_destination

    def _stage_inputs(
        self,
        run_dir: Path,
        input_paths: Optional[Iterable[Path]],
    ) -> Iterable[Path]:
        """Copy inputs into a staging directory for reproducibility."""
        staged_inputs = []
        if not input_paths:
            return staged_inputs
        inputs_dir = ensure_directory(run_dir / "inputs")
        for path in input_paths:
            try:
                destination = inputs_dir / Path(path).name
                shutil.copy2(path, destination)
                staged_inputs.append(destination)
            except (FileNotFoundError, IOError) as exc:
                raise FileProcessingError(
                    f"Failed to stage input file: {path}",
                    filepath=str(path),
                    suggestion="Check that input file exists and is readable",
                ) from exc
        return staged_inputs

    def _metadata_path(self, run_dir: Path) -> Path:
        return run_dir / self.METADATA_FILENAME

    def _persist_metadata(self, run_info: PipelineRunInfo) -> None:
        """Write metadata JSON file for a run."""
        try:
            metadata_path = self._metadata_path(run_info.workdir)
            metadata_path.write_text(json.dumps(run_info.to_json(), indent=2))
        except IOError as exc:
            raise FileProcessingError(
                f"Failed to write metadata file: {metadata_path}",
                filepath=str(metadata_path),
                suggestion="Check write permissions in run directory",
            ) from exc
