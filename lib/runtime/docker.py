"""Docker runtime handler for executing Snakemake."""

import os
from pathlib import Path, PurePosixPath
from typing import List, Optional

from .base import RuntimeHandler


class DockerRuntime(RuntimeHandler):
    """Handler for executing Snakemake in a Docker container."""

    def __init__(
        self, repo_root: Path, container_repo: PurePosixPath = PurePosixPath("/workspace")
    ):
        super().__init__(repo_root, container_repo)

    def get_type(self) -> str:
        return "docker"

    def validate(self) -> bool:
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
        image_ref = image or "metaworks:latest"
        if image_ref.startswith("docker://"):
            image_ref = image_ref[len("docker://"):]

        binds = self._build_binds(bind_paths or [])

        snakefile_in_container = self._container_path(snakefile)
        config_in_container = self._container_path(config_path)
        workdir_container = self._container_path(workdir)

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

        cmd = ["docker", "run", "--rm"]
        if hasattr(os, "getuid") and hasattr(os, "getgid"):
            cmd += ["--user", f"{os.getuid()}:{os.getgid()}"]
        for bind in binds:
            cmd += ["-v", bind]
        cmd += [
            "-e",
            f"HOME={workdir_container}",
            "-e",
            f"XDG_CACHE_HOME={workdir_container}/.cache",
        ]
        cmd += ["-w", workdir_container, image_ref] + snake_args

        return cmd
