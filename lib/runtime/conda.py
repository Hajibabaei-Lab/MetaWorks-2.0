"""Conda runtime handler for executing Snakemake."""

import sys
from pathlib import Path
from typing import List

from .base import RuntimeHandler


class CondaRuntime(RuntimeHandler):
    """Handler for executing Snakemake in a Conda environment."""

    def get_type(self) -> str:
        """Return the runtime type identifier."""
        return "conda"

    def validate(self) -> bool:
        """
        Validate that Conda is available.

        Returns:
            True if conda is available, False otherwise
        """
        try:
            import conda  # noqa: F401

            return True
        except ImportError:
            return False

    def build_command(
        self, snakefile: Path, config_path: Path, workdir: Path, cores: int, **kwargs
    ) -> List[str]:
        """
        Build the command to execute Snakemake with Conda.

        Args:
            snakefile: Path to the Snakefile
            config_path: Path to the configuration file
            workdir: Working directory for the run
            cores: Number of cores to use
            **kwargs: Additional parameters (ignored for conda)

        Returns:
            List of command arguments
        """
        base_args = [
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
        return base_args
