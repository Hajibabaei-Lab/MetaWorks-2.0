"""Docker runtime handler for executing Snakemake."""

from pathlib import Path, PurePosixPath
from typing import List, Optional

from lib.exceptions import RuntimeError as MWRuntimeError

from .base import RuntimeHandler


class DockerRuntime(RuntimeHandler):
    """Handler for executing Snakemake in a Docker container."""

    def __init__(
        self, repo_root: Path, container_repo: PurePosixPath = PurePosixPath("/workspace")
    ):
        """
        Initialize Docker runtime handler.

        Args:
            repo_root: Path to the repository root
            container_repo: Path inside the container corresponding to repo_root
        """
        super().__init__(repo_root)
        self.container_repo = container_repo

    def get_type(self) -> str:
        """Return the runtime type identifier."""
        return "docker"

    def validate(self) -> bool:
        """
        Validate that Docker is available.

        Returns:
            True if Docker is available, False otherwise
        """
        try:
            import subprocess

            result = subprocess.run(
                ["docker", "--version"], capture_output=True, text=True, check=False
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
        **kwargs,
    ) -> List[str]:
        """
        Build the command to execute Snakemake with Docker.

        Args:
            snakefile: Path to the Snakefile
            config_path: Path to the configuration file
            workdir: Working directory for the run
            cores: Number of cores to use
            image: Docker image to use (optional, uses default if not provided)
            bind_paths: Additional bind paths to mount (optional)
            **kwargs: Additional parameters (ignored)

        Returns:
            List of command arguments

        Raises:
            MWRuntimeError: If paths are not under repo_root
        """
        # Normalize image reference
        image_ref = image or "metaworks:latest"
        if image_ref.startswith("docker://"):
            image_ref = image_ref[len("docker://") :]

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

        # Build Docker command
        cmd = ["docker", "run", "--rm"]
        for bind in binds:
            cmd += ["-v", bind]
        cmd += ["-w", workdir_container, image_ref] + snake_args

        return cmd

    def _build_binds(self, additional_binds: List[str]) -> List[str]:
        """
        Build list of Docker bind mounts.

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
                f"for Docker runtime.",
                suggestion="Move it under the repo or provide an explicit bind mapping.",
            ) from exc
        return str(self.container_repo / rel.as_posix())
