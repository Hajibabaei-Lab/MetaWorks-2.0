"""Singularity/Apptainer runtime handler for executing Snakemake."""

import subprocess
from pathlib import Path, PurePosixPath
from typing import List, Optional

from lib.exceptions import RuntimeError as MWRuntimeError

from .base import RuntimeHandler


class SingularityRuntime(RuntimeHandler):
    """Handler for executing Snakemake in a Singularity/Apptainer container."""

    def __init__(
        self,
        repo_root: Path,
        container_repo: PurePosixPath = PurePosixPath("/workspace"),
        cache_dir: Optional[Path] = None,
    ):
        """
        Initialize Singularity runtime handler.

        Args:
            repo_root: Path to the repository root
            container_repo: Path inside the container corresponding to repo_root
            cache_dir: Directory for Singularity cache (optional)
        """
        super().__init__(repo_root)
        self.container_repo = container_repo
        self.cache_dir = cache_dir

    def get_type(self) -> str:
        """Return the runtime type identifier."""
        return "singularity"

    def validate(self) -> bool:
        """
        Validate that Singularity/Apptainer is available.

        Returns:
            True if Singularity/Apptainer is available, False otherwise
        """
        try:
            # Try apptainer first (newer name)
            result = subprocess.run(
                ["apptainer", "--version"], capture_output=True, text=True, check=False
            )
            if result.returncode == 0:
                return True
            # Fall back to singularity
            result = subprocess.run(
                ["singularity", "--version"], capture_output=True, text=True, check=False
            )
            return result.returncode == 0
        except (ImportError, FileNotFoundError):
            return False

    def build_command(
        self,
        snakefile: Path,
        config_path: Path,
        workdir: Path,
        cores: int,
        image: Optional[str] = None,
        bind_paths: Optional[List[str]] = None,
        cache_dir: Optional[Path] = None,
        **kwargs,
    ) -> List[str]:
        """
        Build the command to execute Snakemake with Singularity/Apptainer.

        Args:
            snakefile: Path to the Snakefile
            config_path: Path to the configuration file
            workdir: Working directory for the run
            cores: Number of cores to use
            image: Singularity image to use (optional, uses default if not provided)
            bind_paths: Additional bind paths to mount (optional)
            cache_dir: Singularity cache directory (optional)
            **kwargs: Additional parameters (ignored)

        Returns:
            List of command arguments

        Raises:
            MWRuntimeError: If paths are not under repo_root
        """
        # Determine executable name (apptainer or singularity)
        executable = self._get_executable()

        # Normalize image reference
        image_ref = image or "metaworks.sif"

        # Build bind mounts
        binds = self._build_binds(bind_paths or [])

        # Map paths to container
        snakefile_in_container = self._container_path(snakefile)
        config_in_container = self._container_path(config_path)
        workdir_container = self._container_path(workdir)

        # Build Snakemake arguments
        snake_args = [
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

        # Build Singularity command
        cmd = [executable, "exec", "-B", ",".join(binds)]
        cmd += ["--pwd", workdir_container]

        # Add cache directory if specified
        cache_path = cache_dir or self.cache_dir
        if cache_path:
            cache_path = Path(cache_path)
            cache_path.mkdir(parents=True, exist_ok=True)
            cmd += ["--cache-dir", str(cache_path)]

        cmd += [image_ref] + snake_args

        return cmd

    def _get_executable(self) -> str:
        """
        Determine which executable to use (apptainer or singularity).

        Returns:
            Executable name
        """
        try:
            result = subprocess.run(
                ["apptainer", "--version"], capture_output=True, text=True, check=False
            )
            if result.returncode == 0:
                return "apptainer"
        except (ImportError, FileNotFoundError):
            pass
        return "singularity"

    def _build_binds(self, additional_binds: List[str]) -> List[str]:
        """
        Build list of Singularity bind mounts.

        Args:
            additional_binds: Additional bind paths to include

        Returns:
            List of bind mount strings
        """
        binds = []
        # Add repo bind
        repo_bind = f"{self._normalize_path(self.repo_root)}:{self.container_repo}"
        binds.append(repo_bind)
        # Add additional binds
        for bind in additional_binds:
            if bind and bind not in binds:
                binds.append(bind)
        return binds

    def _container_path(self, path: Path) -> str:
        """
        Map a host path under the repo into the container workspace.

        Args:
            path: Host path to map

        Returns:
            Path inside the container

        Raises:
            MWRuntimeError: If path is not under repo_root
        """
        try:
            rel = path.resolve().relative_to(self.repo_root.resolve())
        except ValueError as exc:
            raise MWRuntimeError(
                f"Path {path} must live under the repository root {self.repo_root} "
                f"for Singularity runtime.",
                suggestion="Move it under the repo or provide an explicit bind mapping.",
            ) from exc
        return str(self.container_repo / rel.as_posix())
