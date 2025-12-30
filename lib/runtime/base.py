"""Base runtime handler interface."""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import List


class RuntimeHandler(ABC):
    """Abstract base class for runtime handlers."""

    def __init__(self, repo_root: Path):
        """
        Initialize runtime handler.

        Args:
            repo_root: Path to the repository root
        """
        self.repo_root = repo_root

    @abstractmethod
    def build_command(
        self, snakefile: Path, config_path: Path, workdir: Path, cores: int, **kwargs
    ) -> List[str]:
        """
        Build the command to execute Snakemake in this runtime.

        Args:
            snakefile: Path to the Snakefile
            config_path: Path to the configuration file
            workdir: Working directory for the run
            cores: Number of cores to use
            **kwargs: Additional runtime-specific parameters

        Returns:
            List of command arguments
        """
        pass

    @abstractmethod
    def get_type(self) -> str:
        """Return the runtime type identifier."""
        pass

    @abstractmethod
    def validate(self) -> bool:
        """Validate that the runtime is available and properly configured."""
        pass

    def _normalize_path(self, path: Path) -> str:
        """
        Normalize a path for use in commands.

        Args:
            path: Path to normalize

        Returns:
            Normalized path string
        """
        return path.resolve().as_posix()
