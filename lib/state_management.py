"""State management for pipeline runs."""

import fcntl
import json
from dataclasses import asdict, dataclass
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional, cast

from lib.exceptions import StateManagementError


class RunStatus(str, Enum):
    """Lifecycle states for a workflow execution."""

    QUEUED = "queued"
    STAGING = "staging"
    RUNNING = "running"
    CANCELLED = "cancelled"
    FAILED = "failed"
    COMPLETED = "completed"


@dataclass
class RunMetadata:
    """Metadata for an individual run."""

    run_id: str
    workflow: str
    runtime: str
    status: str
    submitted_at: str
    run_dir: str
    run_name: Optional[str] = None
    sample_source: Optional[str] = None
    config_overrides: Optional[Dict[str, Any]] = None
    input_dir: Optional[str] = None
    samples_csv: Optional[str] = None
    notes: Optional[str] = None
    scheduler_job_id: Optional[str] = None
    pid: Optional[int] = None
    log_path: Optional[str] = None
    started_at: Optional[str] = None
    ended_at: Optional[str] = None
    return_code: Optional[int] = None
    message: Optional[str] = None
    command: Optional[str] = None
    artifact_path: Optional[str] = None
    keep_outputs: Optional[bool] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "RunMetadata":
        """Create from dictionary."""
        return cls(**data)


class RunStateStore:
    """JSON-backed run metadata store."""

    def __init__(self, state_file: Path):
        """
        Initialize run state store.

        Args:
            state_file: Path to the state JSON file
        """
        self.state_file = state_file
        state_file.parent.mkdir(parents=True, exist_ok=True)
        if not state_file.exists():
            self._write({})
        self._data = self._read()

    def _read(self) -> Dict[str, Any]:
        """Read state from disk with file locking."""
        try:
            with self.state_file.open("r") as fh:
                fcntl.flock(fh, fcntl.LOCK_SH)
                try:
                    return cast(Dict[str, Any], json.load(fh))
                finally:
                    fcntl.flock(fh, fcntl.LOCK_UN)
        except json.JSONDecodeError as exc:
            raise StateManagementError(
                f"Invalid JSON in state file: {self.state_file}",
                run_id=None,
                suggestion="Check state file integrity or delete it to start fresh",
            ) from exc
        except FileNotFoundError:
            return {}

    def _write(self, data: Dict[str, Any]) -> None:
        """Write state to disk with exclusive file locking."""
        self.state_file.parent.mkdir(parents=True, exist_ok=True)
        with self.state_file.open("w") as fh:
            fcntl.flock(fh, fcntl.LOCK_EX)
            try:
                json.dump(data, fh, indent=2, default=str)
            finally:
                fcntl.flock(fh, fcntl.LOCK_UN)

    def create_run(self, run_id: str, metadata: RunMetadata) -> None:
        """
        Create a new run record.

        Args:
            run_id: Unique run identifier
            metadata: Run metadata to store

        Raises:
            StateManagementError: If run_id already exists
        """
        if run_id in self._data:
            raise StateManagementError(
                f"Run ID {run_id} already exists",
                run_id=run_id,
                suggestion="Use a different run ID or delete existing run",
            )
        self._data[run_id] = metadata.to_dict()
        self._write(self._data)

    def update_run(self, run_id: str, updates: Dict[str, Any]) -> None:
        """
        Update an existing run record.

        Args:
            run_id: Run identifier
            updates: Dictionary of fields to update

        Raises:
            StateManagementError: If run_id doesn't exist
        """
        if run_id not in self._data:
            raise StateManagementError(
                f"Unknown run_id {run_id}",
                run_id=run_id,
                suggestion="Check the run ID or list all runs",
            )
        self._data[run_id].update(updates)
        self._write(self._data)

    def get_run(self, run_id: str) -> Optional[RunMetadata]:
        """
        Get a run record.

        Args:
            run_id: Run identifier

        Returns:
            RunMetadata if found, None otherwise
        """
        data = self._data.get(run_id)
        if data:
            return RunMetadata.from_dict(data)
        return None

    def list_runs(self) -> Dict[str, RunMetadata]:
        """
        List all runs.

        Returns:
            Dictionary mapping run_id to RunMetadata
        """
        runs = {}
        for run_id, data in self._data.items():
            runs[run_id] = RunMetadata.from_dict(data)
        return runs

    def delete_run(self, run_id: str) -> None:
        """
        Delete a run record.

        Args:
            run_id: Run identifier

        Raises:
            StateManagementError: If run_id doesn't exist
        """
        if run_id not in self._data:
            raise StateManagementError(
                f"Unknown run_id {run_id}", run_id=run_id, suggestion="Check the run ID"
            )
        del self._data[run_id]
        self._write(self._data)

    def clear_all(self) -> None:
        """Clear all run records."""
        self._data = {}
        self._write(self._data)

    def count_runs(self, status: Optional[RunStatus] = None) -> int:
        """
        Count runs by status.

        Args:
            status: Optional status filter. If None, counts all runs.

        Returns:
            Number of runs
        """
        if status is None:
            return len(self._data)
        return sum(1 for run in self._data.values() if run.get("status") == status.value)

    def find_runs_by_status(self, status: RunStatus) -> List[RunMetadata]:
        """
        Find all runs with a given status.

        Args:
            status: Status to filter by

        Returns:
            List of RunMetadata objects
        """
        return [
            RunMetadata.from_dict(data)
            for data in self._data.values()
            if data.get("status") == status.value
        ]

    def get_recent_runs(self, limit: int = 10) -> List[RunMetadata]:
        """
        Get most recent runs.

        Args:
            limit: Maximum number of runs to return

        Returns:
            List of RunMetadata objects, sorted by submitted_at (newest first)
        """
        runs = [
            (data.get("submitted_at", ""), RunMetadata.from_dict(data))
            for data in self._data.values()
        ]
        # Sort by submitted_at, reverse to get newest first
        runs.sort(key=lambda x: x[0], reverse=True)
        return [metadata for _, metadata in runs[:limit]]
