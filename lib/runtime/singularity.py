"""Singularity/Apptainer runtime handler for executing Snakemake."""

import subprocess
from pathlib import Path, PurePosixPath
from typing import List, Optional

from .base import RuntimeHandler


class SingularityRuntime(RuntimeHandler):
    """Handler for executing Snakemake in a Singularity/Apptainer container."""

    def __init__(
        self,
        repo_root: Path,
        container_repo: PurePosixPath = PurePosixPath("/workspace"),
        cache_dir: Optional[Path] = None,
    ):
        super().__init__(repo_root, container_repo)
        self.cache_dir = cache_dir

    def get_type(self) -> str:
        return "apptainer"

    def validate(self) -> bool:
        try:
            result = subprocess.run(
                ["apptainer", "--version"], capture_output=True, text=True, check=False
            )
            if result.returncode == 0:
                return True
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
        executable = self._get_executable()

        image_ref = image or "metaworks.sif"

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

        cmd = [executable, "exec", "-B", ",".join(binds)]
        cmd += ["--pwd", workdir_container]

        cache_path = cache_dir or self.cache_dir
        if cache_path:
            cache_path = Path(cache_path)
            cache_path.mkdir(parents=True, exist_ok=True)
            cmd += ["--cache-dir", str(cache_path)]

        cmd += [image_ref] + snake_args

        return cmd

    def _get_executable(self) -> str:
        try:
            result = subprocess.run(
                ["apptainer", "--version"], capture_output=True, text=True, check=False
            )
            if result.returncode == 0:
                return "apptainer"
        except (ImportError, FileNotFoundError):
            pass
        return "singularity"
