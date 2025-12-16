import json
import os
import shlex
import shutil
import signal
import subprocess
import sys
import threading
from datetime import datetime
from pathlib import Path, PurePosixPath
from typing import Any, Dict, List, Optional

import yaml

from .config import settings
from .schemas import (
    ArtifactResponse,
    RenderConfigResponse,
    RunStatus,
    RunSubmissionRequest,
    RuntimeType,
    WorkflowType,
)


def _deep_merge(base: Dict[str, Any], overrides: Dict[str, Any]) -> Dict[str, Any]:
    """Recursively merge overrides into base."""
    merged = dict(base)
    for key, value in overrides.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            merged[key] = _deep_merge(base[key], value)
        else:
            merged[key] = value
    return merged


class RunStateStore:
    """JSON-backed run metadata store."""

    def __init__(self, state_file: Path):
        self.state_file = state_file
        state_file.parent.mkdir(parents=True, exist_ok=True)
        if not state_file.exists():
            self._write({})
        self._data = self._read()

    def _read(self) -> Dict[str, Any]:
        with self.state_file.open() as fh:
            return json.load(fh)

    def _write(self, data: Dict[str, Any]) -> None:
        self.state_file.parent.mkdir(parents=True, exist_ok=True)
        with self.state_file.open("w") as fh:
            json.dump(data, fh, indent=2, default=str)

    def set(self, run_id: str, payload: Dict[str, Any]) -> None:
        self._data[run_id] = payload
        self._write(self._data)

    def get(self, run_id: str) -> Optional[Dict[str, Any]]:
        return self._data.get(run_id)

    def all(self) -> Dict[str, Any]:
        return self._data

    def delete(self, run_id: str) -> None:
        if run_id in self._data:
            del self._data[run_id]
            self._write(self._data)


class SnakemakeRunner:
    """Runtime-aware Snakemake command builder and launcher."""

    def __init__(self):
        self.repo_root = settings.repo_root
        self.container_repo = PurePosixPath("/workspace")

    def _snakefile_for(self, workflow: WorkflowType) -> Path:
        snakefile = settings.snakefiles.get(workflow.value)
        if snakefile is None:
            raise ValueError(f"No snakefile configured for workflow {workflow}")
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

    def _container_path(self, path: Path) -> str:
        """Map a host path under the repo into the container workspace."""
        try:
            rel = path.resolve().relative_to(self.repo_root.resolve())
        except ValueError as exc:
            raise ValueError(
                f"Path {path} must live under the repository root {self.repo_root} for container runtimes. "
                "Move it under the repo or provide an explicit bind mapping."
            ) from exc
        return str(self.container_repo / rel.as_posix())

    def _normalize_host_path(self, path: Path) -> str:
        """Return a POSIX-ish string for docker binds even on Windows paths."""
        return path.resolve().as_posix()

    def build_command(self, submission: RunSubmissionRequest, config_path: Path, run_dir: Path) -> List[str]:
        snakefile = self._snakefile_for(submission.workflow)
        cores = submission.cores or settings.default_cores
        workdir_host = run_dir.resolve()

        base_snake_args = [
            "snakemake",
            "-s",
            str(snakefile),
            "--configfile",
            str(config_path),
            "--cores",
            str(cores),
            "--printshellcmds",
            "--directory",
            str(workdir_host),
        ]

        if submission.runtime == RuntimeType.conda:
            return [sys.executable, "-m"] + base_snake_args

        container_image = submission.container_image or settings.container_image
        repo_bind = f"{self._normalize_host_path(self.repo_root)}:{self.container_repo}"
        binds = [repo_bind] + self._resolve_binds(run_dir, submission)

        snakefile_in_container = self._container_path(snakefile)
        config_in_container = self._container_path(config_path)
        workdir_container = self._container_path(workdir_host)

        base_snake_args = [
            "snakemake",
            "-s",
            snakefile_in_container,
            "--configfile",
            config_in_container,
            "--cores",
            str(cores),
            "--printshellcmds",
            "--directory",
            workdir_container,
        ]

        if submission.runtime == RuntimeType.docker:
            image_ref = container_image
            if image_ref.startswith("docker://"):
                image_ref = image_ref[len("docker://") :]
            cmd: List[str] = ["docker", "run", "--rm"]
            for bind in binds:
                cmd += ["-v", bind]
            cmd += ["-w", workdir_container, image_ref] + base_snake_args
            return cmd

        cache_dir = Path(submission.cache_dir or settings.singularity_cache)
        cache_dir.mkdir(parents=True, exist_ok=True)
        cmd = [
            "apptainer",
            "exec",
            "-B",
            ",".join(binds),
            "--pwd",
            workdir_container,
            "--cache-dir",
            str(cache_dir),
            container_image,
        ] + base_snake_args
        return cmd

    def launch(self, command: List[str], log_path: Path) -> subprocess.Popen:
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
        self.store = RunStateStore(settings.state_file)
        self.runner = SnakemakeRunner()
        settings.run_root.mkdir(parents=True, exist_ok=True)
        settings.artifact_root.mkdir(parents=True, exist_ok=True)
        settings.staging_root.mkdir(parents=True, exist_ok=True)
        settings.singularity_cache.mkdir(parents=True, exist_ok=True)

    def _load_default_config(self, workflow: WorkflowType) -> Dict[str, Any]:
        path = settings.default_configs.get(workflow.value)
        if not path:
            raise ValueError(f"No default config registered for workflow {workflow}")
        with Path(path).open() as fh:
            return yaml.safe_load(fh)

    def render_config(self, workflow: WorkflowType, overrides: Dict[str, Any]) -> RenderConfigResponse:
        base_config = self._load_default_config(workflow)
        merged = _deep_merge(base_config, overrides or {})
        rendered = yaml.safe_dump(merged, sort_keys=False)
        return RenderConfigResponse(workflow=workflow, merged=rendered)

    def default_config_sections(self, workflow: WorkflowType) -> Dict[str, Any]:
        return self._load_default_config(workflow)

    def _write_rendered_config(self, run_dir: Path, workflow: WorkflowType, overrides: Dict[str, Any]) -> Path:
        rendered = self.render_config(workflow, overrides)
        config_path = run_dir / "rendered_config.yaml"
        config_path.write_text(rendered.merged)
        return config_path

    def submit_run(self, submission: RunSubmissionRequest) -> RunStatus:
        timestamp = datetime.utcnow().strftime("%Y%m%d%H%M%S")
        run_id = f"{submission.workflow.value}-{timestamp}"
        run_dir = settings.run_root / run_id
        run_dir.mkdir(parents=True, exist_ok=True)
        try:
            runtime = submission.runtime or RuntimeType(settings.default_runtime)
        except ValueError as exc:
            raise ValueError(f"Unsupported runtime '{settings.default_runtime}'") from exc
        submission = submission.copy(update={"runtime": runtime})

        # Record run metadata early
        record = {
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
        self.store.set(run_id, record)

        config_path = self._write_rendered_config(
            run_dir, submission.workflow, submission.config_overrides
        )

        command = self.runner.build_command(submission, config_path, run_dir)
        record["command"] = " ".join(shlex.quote(token) for token in command)

        if submission.dry_run:
            record["status"] = "staged"
            record["message"] = f"Dry run. Command: {record['command']}"
            self.store.set(run_id, record)
            return RunStatus(**record)

        log_path = Path(record["log_path"])
        process = self.runner.launch(command, log_path)
        record["pid"] = process.pid
        record["status"] = "running"
        record["started_at"] = datetime.utcnow().isoformat()
        record["message"] = f"Running with PID {process.pid}"
        self.store.set(run_id, record)

        monitor_thread = threading.Thread(
            target=self._watch_process, args=(run_id, process), daemon=True
        )
        monitor_thread.start()
        return RunStatus(**record)

    def _watch_process(self, run_id: str, process: subprocess.Popen) -> None:
        """Wait for a process and update run metadata."""
        return_code = process.wait()
        record = self.store.get(run_id)
        if not record:
            return
        record["return_code"] = return_code
        record["ended_at"] = datetime.utcnow().isoformat()
        record["status"] = "completed" if return_code == 0 else "failed"
        record["message"] = (
            f"Exited with code {return_code}" if return_code != 0 else "Run completed successfully"
        )
        self.store.set(run_id, record)

    def get_status(self, run_id: str) -> RunStatus:
        record = self.store.get(run_id)
        if not record:
            raise ValueError(f"Unknown run_id {run_id}")
        # Refresh status if the process ended since the last check.
        if record.get("status") == "running" and record.get("pid"):
            pid = int(record["pid"])
            if not self._pid_alive(pid):
                record["status"] = "completed" if record.get("return_code", 1) == 0 else "failed"
                record["ended_at"] = record.get("ended_at") or datetime.utcnow().isoformat()
                record["message"] = record.get("message") or "Process ended"
                self.store.set(run_id, record)
        return RunStatus(**record)

    def cancel_run(self, run_id: str) -> RunStatus:
        record = self.store.get(run_id)
        if not record:
            raise ValueError(f"Unknown run_id {run_id}")
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
                record["message"] = f"Could not terminate PID {pid}: {exc}"
        record["status"] = "cancelled"
        record["ended_at"] = datetime.utcnow().isoformat()
        self.store.set(run_id, record)
        return RunStatus(**record)

    def _pid_alive(self, pid: int) -> bool:
        try:
            os.kill(pid, 0)
        except OSError:
            return False
        return True

    def tail_logs(self, run_id: str, lines: Optional[int] = None) -> Dict[str, Any]:
        record = self.store.get(run_id)
        if not record:
            raise ValueError(f"Unknown run_id {run_id}")
        run_dir = Path(record["run_dir"])
        log_dir = run_dir / "logs"
        target_log = Path(record.get("log_path") or (log_dir / "snakemake.log"))
        if not target_log.exists():
            return {"run_id": run_id, "tail": ["log not yet available"]}

        num_lines = lines or settings.log_tail_lines
        content = target_log.read_text().splitlines()[-num_lines:]
        return {"run_id": run_id, "tail": content}

    def package_artifacts(self, run_id: str) -> ArtifactResponse:
        record = self.store.get(run_id)
        if not record:
            raise ValueError(f"Unknown run_id {run_id}")
        run_dir = Path(record["run_dir"])
        archive_path = settings.artifact_root / f"{run_id}.tar.gz"
        archive_path.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            f"tar -czf {archive_path} -C {run_dir} .",
            shell=True,
            check=False,
        )
        record["artifact_path"] = str(archive_path)
        self.store.set(run_id, record)
        return ArtifactResponse(run_id=run_id, archive_path=str(archive_path))

    def delete_run(self, run_id: str) -> Dict[str, Any]:
        record = self.store.get(run_id)
        if not record:
            raise ValueError(f"Unknown run_id {run_id}")

        # Best-effort cancel before deleting (even if status is stale).
        try:
            self.cancel_run(run_id)
        except Exception:
            pass

        removed: List[str] = []

        def _safe_remove(target: Path, allowed_root: Path) -> None:
            nonlocal removed
            try:
                target = target.resolve()
            except FileNotFoundError:
                return
            if allowed_root not in target.parents and target != allowed_root:
                raise ValueError(f"Refusing to delete path outside allowed root: {target}")
            if target.exists():
                if target.is_dir():
                    shutil.rmtree(target)
                else:
                    target.unlink()
                removed.append(str(target))

        run_dir = Path(record["run_dir"])
        try:
            _safe_remove(run_dir, settings.run_root)
        except Exception:
            # If the configured run_dir is nested under run_root, also clean .snakemake within it if present
            pass

        # Remove embedded .snakemake if under run_root/run_id
        snakemake_dir = run_dir / ".snakemake"
        try:
            _safe_remove(snakemake_dir, settings.run_root)
        except Exception:
            pass

        archive_path = record.get("artifact_path")
        if archive_path:
            try:
                _safe_remove(Path(archive_path), settings.artifact_root)
            except Exception:
                pass

        # Clean up any stray Snakemake cache under repo_root if no runs remain.
        try:
            if not self.store.all():
                repo_snakemake = settings.repo_root / ".snakemake"
                _safe_remove(repo_snakemake, settings.repo_root)
        except Exception:
            pass

        self.store.delete(run_id)
        return {"run_id": run_id, "removed_paths": removed}
