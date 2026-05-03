"""Run file staging, cleanup, artifact packaging, and deletion."""

import logging
import os
import shutil
import subprocess
from pathlib import Path
from typing import Any, Dict, List

from lib.exceptions import FileProcessingError

from ..config import settings

logger = logging.getLogger(__name__)


class RunFileService:
    """Manages run directory lifecycle: staging, cleanup, artifacts, deletion."""

    def stage_support_files(self, run_dir: Path) -> None:
        """
        Stage repo-relative resources into the run directory.

        The Snakemake runner uses ``--directory <run_dir>``, so many workflow
        paths like ``workflow/scripts/...``, ``config/...``, or ``tests/...``
        must exist relative to *run_dir*.
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

        _copytree_if_missing(repo_root / "workflow" / "scripts", run_dir / "workflow" / "scripts")
        _copytree_if_missing(repo_root / "config", run_dir / "config")

        _copyfile_if_missing(
            repo_root / "tests" / "adapters_anchored.fasta",
            run_dir / "tests" / "adapters_anchored.fasta",
        )

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

    def package_artifacts(self, run_dir: Path, run_id: str) -> Path:
        """Create a tar.gz archive of the run directory. Returns archive path."""
        archive_path = settings.artifact_root / f"{run_id}.tar.gz"

        try:
            archive_path.parent.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            raise FileProcessingError(
                f"Failed to create artifact directory: {archive_path.parent}",
                filepath=str(archive_path.parent),
                suggestion="Check permissions and disk space",
            ) from exc

        if not run_dir.exists():
            raise FileProcessingError(
                f"Run outputs are no longer available for {run_id}",
                filepath=str(run_dir),
                suggestion="Run outputs may have been cleaned.",
            )

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

        return archive_path

    def cleanup_run_files(self, record: Dict[str, Any]) -> List[str]:
        """Best-effort cleanup of run files under configured roots."""
        removed: List[str] = []

        def _safe_remove(target: Path, allowed_root: Path) -> None:
            nonlocal removed
            try:
                target = target.resolve()
            except FileNotFoundError:
                return
            if allowed_root not in target.parents and target != allowed_root:
                return
            if not target.exists():
                return
            if target.is_dir():
                shutil.rmtree(target)
            else:
                target.unlink()
            removed.append(str(target))

        run_dir_raw = record.get("run_dir")
        if run_dir_raw:
            try:
                _safe_remove(Path(run_dir_raw), settings.run_root)
            except Exception as exc:
                logger.debug("Best-effort run_dir cleanup failed for %s: %s", run_dir_raw, exc)
            try:
                _safe_remove(Path(run_dir_raw) / ".snakemake" / "log", settings.run_root)
            except Exception as exc:
                logger.debug(
                    "Best-effort snakemake log cleanup failed for %s: %s",
                    run_dir_raw,
                    exc,
                )

        log_path = record.get("log_path")
        if log_path:
            try:
                _safe_remove(Path(log_path), settings.run_root)
            except Exception as exc:
                logger.debug("Best-effort log cleanup failed for %s: %s", log_path, exc)

        return removed

    def delete_run_files(self, record: Dict[str, Any]) -> List[str]:
        """Delete all files associated with a run. Returns list of removed paths."""
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
        log_path = record.get("log_path")

        if log_path:
            try:
                _safe_remove(Path(log_path), settings.run_root)
            except Exception as exc:
                logger.debug("Best-effort log removal failed for %s: %s", log_path, exc)

        try:
            _safe_remove(run_dir, settings.run_root)
        except Exception as exc:
            logger.debug("Best-effort run_dir removal failed for %s: %s", run_dir, exc)

        snakemake_dir = run_dir / ".snakemake"
        try:
            _safe_remove(snakemake_dir, settings.run_root)
        except Exception as exc:
            logger.debug("Best-effort snakemake dir removal failed for %s: %s", snakemake_dir, exc)

        snakemake_log_dir = run_dir / ".snakemake" / "log"
        try:
            _safe_remove(snakemake_log_dir, settings.run_root)
        except Exception as exc:
            logger.debug(
                "Best-effort snakemake log dir removal failed for %s: %s",
                snakemake_log_dir,
                exc,
            )

        archive_path = record.get("artifact_path")
        if archive_path:
            try:
                _safe_remove(Path(archive_path), settings.artifact_root)
            except Exception as exc:
                logger.debug("Best-effort archive removal failed for %s: %s", archive_path, exc)

        return removed

    def cleanup_repo_snakemake(self) -> None:
        """Remove stray .snakemake cache under repo_root."""
        repo_snakemake = settings.repo_root / ".snakemake"
        if repo_snakemake.exists():
            try:
                shutil.rmtree(repo_snakemake)
            except Exception as exc:
                logger.debug("Best-effort repo snakemake removal failed: %s", exc)

    def cleanup_after_download(
        self, archive_path: str, record: Dict[str, Any], retention_policy: str
    ) -> List[str]:
        """Remove run files and temporary archive after a successful artifact download."""
        removed: List[str] = []
        if retention_policy == "until_download":
            removed = self.cleanup_run_files(record)

        try:
            archive = Path(archive_path).resolve()
            if settings.artifact_root.resolve() in archive.parents and archive.exists():
                archive.unlink()
                removed.append(str(archive))
        except Exception as exc:
            logger.debug("Best-effort archive cleanup after download failed: %s", exc)

        return removed
