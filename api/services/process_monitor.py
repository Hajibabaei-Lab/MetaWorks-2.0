"""Process lifecycle monitoring — launch, watch, PID management."""

import logging
import os
import signal
from datetime import datetime
from typing import Any, Dict

logger = logging.getLogger(__name__)


class ProcessMonitor:
    """Manages subprocess lifecycle for pipeline runs."""

    @staticmethod
    def is_alive(pid: int) -> bool:
        try:
            os.kill(pid, 0)
        except OSError:
            return False
        return True

    @staticmethod
    def terminate(pid: int) -> str:
        """Terminate a process by PID. Returns status message."""
        try:
            if hasattr(os, "killpg"):
                os.killpg(int(pid), signal.SIGTERM)
            else:
                os.kill(int(pid), signal.SIGTERM)
            return f"Terminated PID {pid}"
        except ProcessLookupError:
            return f"PID {pid} was not running"

    def resolve_status(
        self, record: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Refresh run record if process ended since last check."""
        if record.get("status") == "running" and record.get("pid"):
            pid = int(record["pid"])
            if not self.is_alive(pid):
                return_code = record.get("return_code", 1)
                record["status"] = "completed" if return_code == 0 else "failed"
                record["ended_at"] = record.get("ended_at") or datetime.utcnow().isoformat()
                record["message"] = record.get("message") or "Process ended"
        return record

    def finalize_completed(
        self, run_id: str, return_code: int
    ) -> Dict[str, Any]:
        """Build a finalization record for a completed process."""
        return {
            "return_code": int(return_code),
            "ended_at": datetime.utcnow().isoformat(),
            "status": "completed" if return_code == 0 else "failed",
            "message": (
                f"Exited with code {return_code}"
                if return_code != 0
                else "Run completed successfully"
            ),
        }
