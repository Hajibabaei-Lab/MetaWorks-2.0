"""Conda runtime handler for executing Snakemake."""

import sys
from pathlib import Path
from typing import List

from .base import RuntimeHandler


class CondaRuntime(RuntimeHandler):
    """Handler for executing Snakemake in a Conda environment."""

    def get_type(self) -> str:
        return "conda"

    def validate(self) -> bool:
        try:
            import conda  # noqa: F401

            return True
        except ImportError:
            return False

    def build_command(
        self, snakefile: Path, config_path: Path, workdir: Path, cores: int, **kwargs
    ) -> List[str]:
        return [
            sys.executable,
            "-m",
            "snakemake",
            "-s",
            self._normalize_path(snakefile),
            "--configfile",
            self._normalize_path(config_path),
            "--cores",
            str(cores),
            "--printshellcmds",
            "--directory",
            self._normalize_path(workdir),
        ]
