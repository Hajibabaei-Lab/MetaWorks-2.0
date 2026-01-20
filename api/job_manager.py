import logging
import os
import shlex
import shutil
import signal
import subprocess
import threading
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, cast

import yaml

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

logger = logging.getLogger(__name__)


def _deep_merge(base: Dict[str, Any], overrides: Dict[str, Any]) -> Dict[str, Any]:
    """Recursively merge overrides into base."""
    merged = dict(base)
    for key, value in overrides.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            merged[key] = _deep_merge(base[key], value)
        else:
            merged[key] = value
    return merged


class SnakemakeRunner:
    """Runtime-aware Snakemake command builder and launcher."""

    def __init__(self):
        self.repo_root = settings.repo_root
        self.conda_handler = CondaRuntime(repo_root=settings.repo_root)
        self.docker_handler = DockerRuntime(repo_root=settings.repo_root)
        self.singularity_handler = SingularityRuntime(repo_root=settings.repo_root)

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
        """Collect explicit bind strings (src:dest) supplied via settings or submission."""
        binds: List[str] = []
        for item in settings.bind_paths + (submission.bind_paths or []):
            if item and item not in binds:
                binds.append(item)
        return binds

    def build_command(
        self, submission: RunSubmissionRequest, config_path: Path, run_dir: Path
    ) -> List[str]:
        """Build Snakemake command using appropriate runtime handler."""
        snakefile = self._snakefile_for(submission.workflow)
        cores = submission.cores or settings.default_cores

        # Select appropriate runtime handler
        if submission.runtime == RuntimeType.conda:
            return self.conda_handler.build_command(
                snakefile=snakefile, config_path=config_path, workdir=run_dir, cores=cores
            )

        # For container runtimes, get image and additional bind paths
        container_image = submission.container_image or settings.container_image
        additional_binds = self._resolve_binds(run_dir, submission)

        if submission.runtime == RuntimeType.docker:
            return self.docker_handler.build_command(
                snakefile=snakefile,
                config_path=config_path,
                workdir=run_dir,
                cores=cores,
                image=container_image,
                additional_binds=additional_binds,
            )

        # Singularity runtime
        cache_dir = submission.cache_dir or settings.singularity_cache
        return self.singularity_handler.build_command(
            snakefile=snakefile,
            config_path=config_path,
            workdir=run_dir,
            cores=cores,
            image=container_image,
            additional_binds=additional_binds,
            cache_dir=Path(cache_dir),
        )

    def launch(self, command: List[str], log_path: Path) -> subprocess.Popen:
        """Launch a subprocess with proper configuration."""
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_file = log_path.open("a")
        # Run from repo root to keep Snakefile relative paths stable.
        popen_kwargs: Dict[str, Any] = {
            "stdout": log_file,
            "stderr": subprocess.STDOUT,
            "cwd": self.repo_root,
        }
        if hasattr(os, "setsid"):
            popen_kwargs["preexec_fn"] = os.setsid  # POSIX
        elif os.name == "nt" and hasattr(subprocess, "CREATE_NEW_PROCESS_GROUP"):
            popen_kwargs["creationflags"] = subprocess.CREATE_NEW_PROCESS_GROUP  # Windows
        return subprocess.Popen(command, **popen_kwargs)


class JobManager:
    def __init__(self):
        self._state_store = RunStateStore(settings.state_file)
        self.runner = SnakemakeRunner()
        settings.run_root.mkdir(parents=True, exist_ok=True)
        settings.artifact_root.mkdir(parents=True, exist_ok=True)
        settings.staging_root.mkdir(parents=True, exist_ok=True)
        settings.singularity_cache.mkdir(parents=True, exist_ok=True)

    def _dict_to_metadata(self, run_id: str, record: Dict[str, Any]) -> RunMetadata:
        """Convert dict record to RunMetadata for storage."""
        workflow = str(record["workflow"])
        runtime = str(record["runtime"])
        status = str(record["status"])
        submitted_at = str(record["submitted_at"])
        run_dir = str(record["run_dir"])
        return RunMetadata(
            run_id=run_id,
            workflow=workflow,
            runtime=runtime,
            status=status,
            submitted_at=submitted_at,
            started_at=record.get("started_at"),
            ended_at=record.get("ended_at"),
            run_dir=run_dir,
            run_name=record.get("run_name"),
            sample_source=record.get("sample_source"),
            config_overrides=record.get("config_overrides"),
            input_dir=record.get("input_dir"),
            samples_csv=record.get("samples_csv"),
            notes=record.get("notes"),
            pid=record.get("pid"),
            log_path=record.get("log_path"),
            command=record.get("command"),
            message=record.get("message"),
            return_code=record.get("return_code"),
            artifact_path=record.get("artifact_path"),
            scheduler_job_id=record.get("scheduler_job_id"),
        )

    def _metadata_to_dict(self, metadata: RunMetadata) -> Dict[str, Any]:
        """Convert RunMetadata back to dict for API compatibility."""
        record = metadata.to_dict()
        # Ensure backward compatibility with field name changes
        if "ended_at" in record and "return_code" in record:
            pass  # Already in correct format
        return record

    def _get_metadata(self, run_id: str) -> Optional[RunMetadata]:
        """Get run metadata from state store."""
        return self._state_store.get_run(run_id)

    def _set_metadata(self, run_id: str, record: Dict[str, Any]) -> None:
        """Set run metadata in state store."""
        metadata = self._dict_to_metadata(run_id, record)
        existing = self._state_store.get_run(run_id)
        if existing:
            self._state_store.update_run(run_id, record)
        else:
            self._state_store.create_run(run_id, metadata)

    def _all_runs(self) -> Dict[str, Any]:
        """Get all runs as dict for compatibility."""
        runs = self._state_store.list_runs()
        return {run_id: self._metadata_to_dict(metadata) for run_id, metadata in runs.items()}

    def _delete_run_state(self, run_id: str) -> None:
        """Delete run from state store."""
        self._state_store.delete_run(run_id)

    def _load_default_config(self, workflow: WorkflowType) -> Dict[str, Any]:
        path = settings.default_configs.get(workflow.value)
        if not path:
            raise ConfigurationError(
                f"No default config registered for workflow {workflow}",
                config_key=f"default_configs.{workflow.value}",
                suggestion="Check settings.default_configs configuration",
            )
        config_path = Path(path)
        if not config_path.is_absolute():
            config_path = settings.repo_root / config_path
        try:
            with config_path.open() as fh:
                loaded = yaml.safe_load(fh)
                if loaded is None:
                    return {}
                if not isinstance(loaded, dict):
                    raise ConfigurationError(
                        f"Config file did not parse to a mapping: {config_path}",
                        filepath=str(config_path),
                        suggestion="Ensure the top-level YAML value is a mapping/object",
                    )
                return cast(Dict[str, Any], loaded)
        except yaml.YAMLError as exc:
            raise ConfigurationError(
                f"Invalid YAML in config file: {path}",
                filepath=str(config_path),
                suggestion=f"Fix YAML syntax in {config_path}",
            ) from exc
        except FileNotFoundError as exc:
            raise ConfigurationError(
                f"Config file not found: {path}",
                filepath=str(config_path),
                suggestion="Verify config file path in settings",
            ) from exc

    def render_config(
        self, workflow: WorkflowType, overrides: Dict[str, Any]
    ) -> RenderConfigResponse:
        base_config = self._load_default_config(workflow)
        merged = _deep_merge(base_config, overrides or {})
        rendered = yaml.safe_dump(merged, sort_keys=False)
        return RenderConfigResponse(workflow=workflow, merged=rendered)

    def default_config_sections(self, workflow: WorkflowType) -> Dict[str, Any]:
        return self._load_default_config(workflow)

    def _write_rendered_config(
        self, run_dir: Path, workflow: WorkflowType, overrides: Dict[str, Any]
    ) -> Path:
        rendered = self.render_config(workflow, overrides)
        config_path = run_dir / "rendered_config.yaml"
        config_path.write_text(rendered.merged)
        return config_path

    def _stage_run_support_files(self, run_dir: Path) -> None:
        """
        Stage repo-relative resources into the run directory.

        The Snakemake runner uses `--directory <run_dir>`, so many workflow paths
        like `python_scripts/...`, `config/...`, or `tests/...` must exist
        relative to `run_dir`.
        """

        repo_root = settings.repo_root

        def _copytree_if_missing(src: Path, dest: Path) -> None:
            if dest.exists() or dest.is_symlink():
                return
            if not src.exists():
                return
            shutil.copytree(src, dest)

        def _copyfile_if_missing(src: Path, dest: Path) -> None:
            if dest.exists() or dest.is_symlink():
                return
            if not src.exists():
                return
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dest)

        # Scripts referenced as `python3 python_scripts/<script>.py` in many rules.
        _copytree_if_missing(repo_root / "python_scripts", run_dir / "python_scripts")

        # Config resources referenced as `config/...` (e.g. HMM profiles).
        _copytree_if_missing(repo_root / "config", run_dir / "config")

        # Default adapter file referenced as `tests/adapters_anchored.fasta`.
        _copyfile_if_missing(
            repo_root / "tests" / "adapters_anchored.fasta",
            run_dir / "tests" / "adapters_anchored.fasta",
        )

        # Runtime-managed assets referenced as `runtime/classifiers/...` or `runtime/adapters/...`.
        # Prefer symlinks (so uploads remain visible), fall back to copying if needed.
        runtime_dir = run_dir / "runtime"
        runtime_dir.mkdir(parents=True, exist_ok=True)

        def _symlink_or_copy_dir(src: Path, dest: Path) -> None:
            if dest.exists() or dest.is_symlink():
                return
            if not src.exists():
                return
            try:
                rel_target = os.path.relpath(src.resolve(), start=dest.parent.resolve())
                dest.symlink_to(rel_target, target_is_directory=True)
            except Exception:
                shutil.copytree(src, dest)

        _symlink_or_copy_dir(settings.classifier_root, runtime_dir / "classifiers")
        _symlink_or_copy_dir(settings.adapter_root, runtime_dir / "adapters")

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
            self._stage_run_support_files(run_dir)
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
        submission = submission.copy(update={"runtime": runtime})

        # Record run metadata early
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
            "config_overrides": submission.config_overrides,
            "input_dir": submission.input_dir,
            "samples_csv": submission.samples_csv,
            "notes": submission.notes,
            "pid": None,
            "log_path": str(run_dir / "logs" / "snakemake.log"),
        }
        self._set_metadata(run_id, record)

        config_path = self._write_rendered_config(
            run_dir, submission.workflow, submission.config_overrides
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
        """Wait for a process and update run metadata."""
        return_code = process.wait()
        metadata = self._get_metadata(run_id)
        if not metadata:
            return
        record = self._metadata_to_dict(metadata)
        record["return_code"] = int(return_code)
        record["ended_at"] = datetime.utcnow().isoformat()
        record["status"] = "completed" if return_code == 0 else "failed"
        record["message"] = (
            f"Exited with code {return_code}" if return_code != 0 else "Run completed successfully"
        )
        self._set_metadata(run_id, record)

    def get_status(self, run_id: str) -> RunStatus:
        metadata = self._get_metadata(run_id)
        if not metadata:
            raise StateManagementError(
                f"Unknown run_id {run_id}",
                run_id=run_id,
                suggestion="Check run ID or list all runs",
            )
        record = self._metadata_to_dict(metadata)
        # Refresh status if process ended since last check.
        if record.get("status") == "running" and record.get("pid"):
            pid = int(record["pid"])
            if not self._pid_alive(pid):
                record["status"] = "completed" if record.get("return_code", 1) == 0 else "failed"
                record["ended_at"] = record.get("ended_at") or datetime.utcnow().isoformat()
                record["message"] = record.get("message") or "Process ended"
                self._set_metadata(run_id, record)
        return RunStatus(**record)

    def cancel_run(self, run_id: str) -> RunStatus:
        metadata = self._get_metadata(run_id)
        if not metadata:
            raise StateManagementError(
                f"Unknown run_id {run_id}",
                run_id=run_id,
                suggestion="Check run ID or list all runs",
            )
        record = self._metadata_to_dict(metadata)
        pid = record.get("pid")
        if pid:
            try:
                if hasattr(os, "killpg"):
                    os.killpg(int(pid), signal.SIGTERM)  # POSIX
                else:
                    os.kill(int(pid), signal.SIGTERM)  # Windows fallback
                record["message"] = f"Terminated PID {pid}"
            except ProcessLookupError:
                record["message"] = f"PID {pid} was not running"
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

    def _pid_alive(self, pid: int) -> bool:
        try:
            os.kill(pid, 0)
        except OSError:
            return False
        return True

    def tail_logs(self, run_id: str, lines: Optional[int] = None) -> Dict[str, Any]:
        metadata = self._get_metadata(run_id)
        if not metadata:
            raise StateManagementError(
                f"Unknown run_id {run_id}",
                run_id=run_id,
                suggestion="Check run ID or list all runs",
            )
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

    def package_artifacts(self, run_id: str) -> ArtifactResponse:
        metadata = self._get_metadata(run_id)
        if not metadata:
            raise StateManagementError(
                f"Unknown run_id {run_id}",
                run_id=run_id,
                suggestion="Check run ID or list all runs",
            )
        record = self._metadata_to_dict(metadata)
        run_dir = Path(record["run_dir"])
        archive_path = settings.artifact_root / f"{run_id}.tar.gz"

        try:
            archive_path.parent.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            raise FileProcessingError(
                f"Failed to create artifact directory: {archive_path.parent}",
                filepath=str(archive_path.parent),
                suggestion="Check permissions and disk space",
            ) from exc

        try:
            subprocess.run(
                ["tar", "-czf", str(archive_path), "-C", str(run_dir), "."],
                check=True,
            )
        except subprocess.CalledProcessError as exc:
            raise FileProcessingError(
                f"Failed to create archive: {archive_path}",
                filepath=str(archive_path),
                suggestion=f"Check disk space and run directory contents: {run_dir}",
            ) from exc

        record["artifact_path"] = str(archive_path)
        self._set_metadata(run_id, record)
        return ArtifactResponse(run_id=run_id, archive_path=str(archive_path))

    def delete_run(self, run_id: str) -> Dict[str, Any]:
        metadata = self._get_metadata(run_id)
        if not metadata:
            raise StateManagementError(
                f"Unknown run_id {run_id}",
                run_id=run_id,
                suggestion="Check run ID or list all runs",
            )
        record = self._metadata_to_dict(metadata)

        # Best-effort cancel before deleting (even if status is stale).
        try:
            self.cancel_run(run_id)
        except Exception as exc:
            logger.debug("Best-effort cancel failed for %s: %s", run_id, exc)

        removed: List[str] = []

        def _safe_remove(target: Path, allowed_root: Path) -> None:
            nonlocal removed
            try:
                target = target.resolve()
            except FileNotFoundError:
                return
            if allowed_root not in target.parents and target != allowed_root:
                raise FileProcessingError(
                    f"Refusing to delete path outside allowed root: {target}",
                    filepath=str(target),
                    suggestion="Verify file paths are within expected directories",
                )
            if target.exists():
                try:
                    if target.is_dir():
                        shutil.rmtree(target)
                    else:
                        target.unlink()
                    removed.append(str(target))
                except OSError as exc:
                    raise FileProcessingError(
                        f"Failed to delete {target}: {exc}",
                        filepath=str(target),
                        suggestion="Check file permissions",
                    ) from exc

        run_dir = Path(record["run_dir"])
        try:
            _safe_remove(run_dir, settings.run_root)
        except Exception as exc:
            # If the configured run_dir is nested under run_root, also clean .snakemake within it if present
            logger.debug("Best-effort run_dir removal failed for %s: %s", run_dir, exc)

        # Remove embedded .snakemake if under run_root/run_id
        snakemake_dir = run_dir / ".snakemake"
        try:
            _safe_remove(snakemake_dir, settings.run_root)
        except Exception as exc:
            logger.debug("Best-effort snakemake dir removal failed for %s: %s", snakemake_dir, exc)

        archive_path = record.get("artifact_path")
        if archive_path:
            try:
                _safe_remove(Path(archive_path), settings.artifact_root)
            except Exception as exc:
                logger.debug("Best-effort archive removal failed for %s: %s", archive_path, exc)

        # Clean up any stray Snakemake cache under repo_root if no runs remain.
        try:
            if not self._all_runs():
                repo_snakemake = settings.repo_root / ".snakemake"
                _safe_remove(repo_snakemake, settings.repo_root)
        except Exception as exc:
            logger.debug("Best-effort repo snakemake removal failed: %s", exc)

        self._delete_run_state(run_id)
        return {"run_id": run_id, "removed_paths": removed}
