"""Thin orchestration facade composing run lifecycle services."""

import logging
import os
import shlex
import subprocess
import threading
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

from lib.exceptions import (
    ConfigurationError,
    FileProcessingError,
    PipelineExecutionError,
    StateManagementError,
)
from lib.runtime import CondaRuntime, DockerRuntime, SingularityRuntime
from lib.state_management import RunMetadata, RunStateStore

from .config import settings
from .schemas import (
    ArtifactResponse,
    RenderConfigResponse,
    RunStatus,
    RunSubmissionRequest,
    RuntimeType,
    WorkflowType,
)
from .services.config_renderer import ConfigRenderer
from .services.process_monitor import ProcessMonitor
from .services.run_file_service import RunFileService

logger = logging.getLogger(__name__)


class SnakemakeRunner:
    """Runtime-aware Snakemake command builder and launcher."""

    def __init__(self):
        self.repo_root = settings.repo_root

    def _handler_for(self, runtime: RuntimeType):
        if runtime == RuntimeType.conda:
            return CondaRuntime(repo_root=settings.repo_root)
        if runtime == RuntimeType.docker:
            return DockerRuntime(repo_root=settings.repo_root)
        if runtime == RuntimeType.apptainer:
            return SingularityRuntime(repo_root=settings.repo_root)
        raise ConfigurationError(
            f"Unknown runtime type: {runtime}",
            config_key="runtime",
        )

    def _snakefile_for(self, workflow: WorkflowType) -> Path:
        snakefile = settings.snakefiles.get(workflow.value)
        if snakefile is None:
            raise ConfigurationError(
                f"No snakefile configured for workflow {workflow}",
                config_key=f"snakefiles.{workflow.value}",
                suggestion=f"Add snakefile path to settings for {workflow.value}",
            )
        snakefile_path = Path(snakefile)
        if not snakefile_path.is_absolute():
            snakefile_path = self.repo_root / snakefile_path
        return snakefile_path

    def _resolve_binds(self, run_dir: Path, submission: RunSubmissionRequest) -> List[str]:
        binds: List[str] = []
        for item in settings.bind_paths + (submission.bind_paths or []):
            if item and item not in binds:
                binds.append(item)
        return binds

    def build_command(
        self, submission: RunSubmissionRequest, config_path: Path, run_dir: Path
    ) -> List[str]:
        snakefile = self._snakefile_for(submission.workflow)
        cores = submission.cores or settings.default_cores
        handler = self._handler_for(submission.runtime)

        if submission.runtime == RuntimeType.conda:
            return handler.build_command(
                snakefile=snakefile, config_path=config_path, workdir=run_dir, cores=cores
            )

        container_image = submission.container_image or settings.container_image
        additional_binds = self._resolve_binds(run_dir, submission)

        if submission.runtime == RuntimeType.apptainer:
            cache_dir = submission.cache_dir or settings.singularity_cache
            return handler.build_command(
                snakefile=snakefile,
                config_path=config_path,
                workdir=run_dir,
                cores=cores,
                image=container_image,
                bind_paths=additional_binds,
                cache_dir=Path(cache_dir),
            )

        return handler.build_command(
            snakefile=snakefile,
            config_path=config_path,
            workdir=run_dir,
            cores=cores,
            image=container_image,
            bind_paths=additional_binds,
        )

    def launch(self, command: List[str], log_path: Path) -> subprocess.Popen:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_file = log_path.open("a")
        popen_kwargs: Dict[str, Any] = {
            "stdout": log_file,
            "stderr": subprocess.STDOUT,
            "cwd": self.repo_root,
        }
        if hasattr(os, "setsid"):
            popen_kwargs["preexec_fn"] = os.setsid
        elif os.name == "nt" and hasattr(subprocess, "CREATE_NEW_PROCESS_GROUP"):
            popen_kwargs["creationflags"] = subprocess.CREATE_NEW_PROCESS_GROUP
        return subprocess.Popen(command, **popen_kwargs)


class JobManager:
    """Facade that composes config rendering, file management, process monitoring, and state."""

    def __init__(self):
        self._state_store = RunStateStore(settings.state_file)
        self.runner = SnakemakeRunner()
        self._config_renderer = ConfigRenderer(str(settings.repo_root))
        self._file_service = RunFileService()
        self._process_monitor = ProcessMonitor()

        settings.run_root.mkdir(parents=True, exist_ok=True)
        settings.artifact_root.mkdir(parents=True, exist_ok=True)
        settings.staging_root.mkdir(parents=True, exist_ok=True)
        settings.singularity_cache.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # State store helpers
    # ------------------------------------------------------------------

    def _dict_to_metadata(self, run_id: str, record: Dict[str, Any]) -> RunMetadata:
        return RunMetadata(
            run_id=run_id,
            workflow=str(record["workflow"]),
            runtime=str(record["runtime"]),
            status=str(record["status"]),
            submitted_at=str(record["submitted_at"]),
            run_dir=str(record["run_dir"]),
            started_at=record.get("started_at"),
            ended_at=record.get("ended_at"),
            run_name=record.get("run_name"),
            sample_source=record.get("sample_source"),
            config_overrides=record.get("config_overrides"),
            input_dir=record.get("input_dir"),
            samples_csv=record.get("samples_csv"),
            notes=record.get("notes"),
            keep_outputs=record.get("keep_outputs"),
            pid=record.get("pid"),
            log_path=record.get("log_path"),
            command=record.get("command"),
            message=record.get("message"),
            return_code=record.get("return_code"),
            artifact_path=record.get("artifact_path"),
            scheduler_job_id=record.get("scheduler_job_id"),
        )

    def _metadata_to_dict(self, metadata: RunMetadata) -> Dict[str, Any]:
        return metadata.to_dict()

    def _get_metadata(self, run_id: str) -> Optional[RunMetadata]:
        return self._state_store.get_run(run_id)

    def _set_metadata(self, run_id: str, record: Dict[str, Any]) -> None:
        metadata = self._dict_to_metadata(run_id, record)
        existing = self._state_store.get_run(run_id)
        if existing:
            self._state_store.update_run(run_id, record)
        else:
            self._state_store.create_run(run_id, metadata)

    def _all_runs(self) -> Dict[str, Any]:
        runs = self._state_store.list_runs()
        return {run_id: self._metadata_to_dict(metadata) for run_id, metadata in runs.items()}

    def _delete_run_state(self, run_id: str) -> None:
        self._state_store.delete_run(run_id)

    def _require_metadata(self, run_id: str) -> RunMetadata:
        metadata = self._get_metadata(run_id)
        if not metadata:
            raise StateManagementError(
                f"Unknown run_id {run_id}",
                run_id=run_id,
                suggestion="Check run ID or list all runs",
            )
        return metadata

    # ------------------------------------------------------------------
    # Config rendering (delegates to ConfigRenderer)
    # ------------------------------------------------------------------

    def render_config_with_profile(
        self, profile: str, workflow: WorkflowType, overrides: Dict[str, Any]
    ) -> RenderConfigResponse:
        return self._config_renderer.render_config(profile, workflow, overrides)

    def render_config(
        self, workflow: WorkflowType, overrides: Dict[str, Any]
    ) -> RenderConfigResponse:
        return self.render_config_with_profile("coi", workflow, overrides)

    def default_config_sections(self, workflow: WorkflowType) -> Dict[str, Any]:
        return self._config_renderer.default_config_sections(workflow)

    def _write_rendered_config(
        self, run_dir: Path, profile: str, workflow: WorkflowType, overrides: Dict[str, Any]
    ) -> Path:
        return self._config_renderer.write_rendered_config(run_dir, profile, workflow, overrides)

    # ------------------------------------------------------------------
    # Run lifecycle
    # ------------------------------------------------------------------

    def submit_run(self, submission: RunSubmissionRequest) -> RunStatus:
        timestamp = datetime.utcnow().strftime("%Y%m%d%H%M%S")
        run_id = f"{submission.workflow.value}-{timestamp}"
        run_dir = settings.run_root / run_id

        try:
            run_dir.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            raise FileProcessingError(
                f"Failed to create run directory: {run_dir}",
                filepath=str(run_dir),
                suggestion="Check permissions and disk space",
            ) from exc

        try:
            self._file_service.stage_support_files(run_dir)
        except OSError as exc:
            raise FileProcessingError(
                f"Failed to stage run support files into: {run_dir}",
                filepath=str(run_dir),
                suggestion="Check permissions and available disk space",
            ) from exc

        try:
            runtime = submission.runtime or RuntimeType(settings.default_runtime)
        except ValueError as exc:
            raise ConfigurationError(
                f"Unsupported runtime '{settings.default_runtime}'",
                config_key="default_runtime",
                suggestion="Use one of: conda, docker, singularity",
            ) from exc
        allowed = self._allowed_runtime_names()
        if runtime.value not in allowed:
            raise ConfigurationError(
                f"Runtime '{runtime.value}' is disabled for web runs",
                config_key="runtime",
                suggestion=f"Use one of: {', '.join(allowed)}",
            )
        submission = submission.copy(update={"runtime": runtime})

        effective_overrides: Dict[str, Any] = dict(submission.config_overrides or {})
        input_block = effective_overrides.get("input")
        if not isinstance(input_block, dict):
            input_block = {}
        input_block["fastq_dir"] = submission.input_dir
        input_block["sample_source"] = submission.sample_source
        if submission.sample_source == "csv" and submission.samples_csv:
            input_block["samples_csv"] = submission.samples_csv
        effective_overrides["input"] = input_block

        record: Dict[str, Any] = {
            "run_id": run_id,
            "workflow": submission.workflow.value,
            "runtime": runtime.value,
            "status": "queued",
            "submitted_at": datetime.utcnow().isoformat(),
            "scheduler_job_id": None,
            "run_dir": str(run_dir),
            "run_name": submission.run_name,
            "sample_source": submission.sample_source,
            "config_overrides": effective_overrides,
            "input_dir": submission.input_dir,
            "samples_csv": submission.samples_csv,
            "notes": submission.notes,
            "keep_outputs": bool(submission.keep_outputs),
            "pid": None,
            "log_path": str(run_dir / "logs" / "snakemake.log"),
        }
        self._set_metadata(run_id, record)

        config_path = self._write_rendered_config(
            run_dir, submission.profile, submission.workflow, effective_overrides
        )

        command = self.runner.build_command(submission, config_path, run_dir)
        record["command"] = " ".join(shlex.quote(token) for token in command)

        if submission.dry_run:
            record["status"] = "staged"
            record["message"] = f"Dry run. Command: {record['command']}"
            self._set_metadata(run_id, record)
            return RunStatus(**record)

        log_path = Path(record["log_path"])
        process = self.runner.launch(command, log_path)
        record["pid"] = process.pid
        record["status"] = "running"
        record["started_at"] = datetime.utcnow().isoformat()
        record["message"] = f"Running with PID {process.pid}"
        self._set_metadata(run_id, record)

        monitor_thread = threading.Thread(
            target=self._watch_process, args=(run_id, process), daemon=True
        )
        monitor_thread.start()
        return RunStatus(**record)

    def _watch_process(self, run_id: str, process: subprocess.Popen) -> None:
        return_code = process.wait()
        metadata = self._get_metadata(run_id)
        if not metadata:
            return
        record = self._metadata_to_dict(metadata)
        finalization = self._process_monitor.finalize_completed(run_id, return_code)
        record.update(finalization)

        retention_policy = self._retention_policy()
        should_cleanup = (
            retention_policy == "immediate"
            or (retention_policy == "until_download" and not bool(record.get("keep_outputs", True)))
        )
        if should_cleanup:
            removed = self._file_service.cleanup_run_files(record)
            if removed:
                record["message"] = f"{record['message']} (cleanup policy: {retention_policy})"
            record["log_path"] = None
            record["artifact_path"] = None
        self._set_metadata(run_id, record)

    # ------------------------------------------------------------------
    # Run status + cancellation
    # ------------------------------------------------------------------

    def get_status(self, run_id: str) -> RunStatus:
        metadata = self._require_metadata(run_id)
        record = self._metadata_to_dict(metadata)
        record = self._process_monitor.resolve_status(record)
        if record.get("status") != metadata.status:
            self._set_metadata(run_id, record)
        return RunStatus(**record)

    def list_runs(self) -> List[RunStatus]:
        runs: List[RunStatus] = []
        for run_id in sorted(self._all_runs().keys(), reverse=True):
            try:
                runs.append(self.get_status(run_id))
            except StateManagementError:
                logger.debug("Skipping stale run metadata for %s during list_runs", run_id)
        return runs

    def cancel_run(self, run_id: str) -> RunStatus:
        metadata = self._require_metadata(run_id)
        record = self._metadata_to_dict(metadata)
        pid = record.get("pid")
        if pid:
            try:
                record["message"] = self._process_monitor.terminate(int(pid))
            except OSError as exc:
                raise PipelineExecutionError(
                    f"Could not terminate PID {pid}: {exc}",
                    stage="cancel",
                    suggestion="Check if process is still running",
                ) from exc
        record["status"] = "cancelled"
        record["ended_at"] = datetime.utcnow().isoformat()
        self._set_metadata(run_id, record)
        return RunStatus(**record)

    # ------------------------------------------------------------------
    # Log streaming
    # ------------------------------------------------------------------

    def tail_logs(self, run_id: str, lines: Optional[int] = None) -> Dict[str, Any]:
        metadata = self._require_metadata(run_id)
        record = self._metadata_to_dict(metadata)
        run_dir = Path(record["run_dir"])
        log_dir = run_dir / "logs"
        target_log = Path(record.get("log_path") or (log_dir / "snakemake.log"))
        if not target_log.exists():
            return {"run_id": run_id, "tail": ["log not yet available"]}

        try:
            num_lines = lines or settings.log_tail_lines
            content = target_log.read_text().splitlines()[-num_lines:]
            return {"run_id": run_id, "tail": content}
        except IOError as exc:
            raise FileProcessingError(
                f"Failed to read log file: {target_log}",
                filepath=str(target_log),
                suggestion="Check file permissions",
            ) from exc

    # ------------------------------------------------------------------
    # Artifact packaging + download cleanup
    # ------------------------------------------------------------------

    def package_artifacts(self, run_id: str) -> ArtifactResponse:
        metadata = self._require_metadata(run_id)
        record = self._metadata_to_dict(metadata)
        run_dir = Path(record["run_dir"])
        archive_path = self._file_service.package_artifacts(run_dir, run_id)
        record["artifact_path"] = str(archive_path)
        self._set_metadata(run_id, record)
        return ArtifactResponse(run_id=run_id, archive_path=str(archive_path))

    def cleanup_after_download(self, run_id: str, archive_path: str) -> None:
        metadata = self._get_metadata(run_id)
        if not metadata:
            return
        record = self._metadata_to_dict(metadata)
        retention_policy = self._retention_policy()
        removed = self._file_service.cleanup_after_download(archive_path, record, retention_policy)

        record["artifact_path"] = None
        if retention_policy == "until_download":
            record["log_path"] = None
            record["keep_outputs"] = False
        if removed:
            record["message"] = (
                (record.get("message") or "Completed")
                + " (artifact downloaded; run files cleaned)"
            )
        self._set_metadata(run_id, record)

    # ------------------------------------------------------------------
    # Run deletion
    # ------------------------------------------------------------------

    def delete_run(self, run_id: str) -> Dict[str, Any]:
        metadata = self._require_metadata(run_id)
        record = self._metadata_to_dict(metadata)

        try:
            self.cancel_run(run_id)
        except Exception as exc:
            logger.debug("Best-effort cancel failed for %s: %s", run_id, exc)

        removed = self._file_service.delete_run_files(record)

        try:
            if not self._all_runs():
                self._file_service.cleanup_repo_snakemake()
        except Exception as exc:
            logger.debug("Best-effort repo snakemake removal failed: %s", exc)

        self._delete_run_state(run_id)
        return {"run_id": run_id, "removed_paths": removed}

    # ------------------------------------------------------------------
    # Settings helpers
    # ------------------------------------------------------------------

    def _allowed_runtime_names(self) -> List[str]:
        values = [
            item.strip().lower()
            for item in str(settings.allowed_runtimes or "").split(",")
            if item.strip()
        ]
        return values or ["docker", "apptainer"]

    def _retention_policy(self) -> str:
        policy = str(settings.retention_policy or "until_download").strip().lower()
        if policy not in {"until_download", "immediate", "manual"}:
            return "until_download"
        return policy
