"""Base runtime handler interface."""

from abc import ABC, abstractmethod
from pathlib import Path, PurePosixPath
from typing import List

from lib.exceptions import MetaWorksRuntimeError as MWRuntimeError


class RuntimeHandler(ABC):
    """Abstract base class for runtime handlers."""

    def __init__(self, repo_root: Path, container_repo: PurePosixPath = None):
        self.repo_root = repo_root
        self.container_repo = container_repo

    @abstractmethod
    def build_command(
        self, snakefile: Path, config_path: Path, workdir: Path, cores: int, **kwargs
    ) -> List[str]:
        pass

    @abstractmethod
    def get_type(self) -> str:
        pass

    @abstractmethod
    def validate(self) -> bool:
        pass

    def _normalize_path(self, path: Path) -> str:
        return path.resolve().as_posix()

    def _build_binds(self, additional_binds: List[str]) -> List[str]:
        binds = []
        if self.container_repo is not None:
            repo_bind = f"{self._normalize_path(self.repo_root)}:{self.container_repo}"
            binds.append(repo_bind)
        for bind in additional_binds:
            if bind and bind not in binds:
                binds.append(bind)
        return binds

    def _container_path(self, path: Path) -> str:
        try:
            rel = path.resolve().relative_to(self.repo_root.resolve())
        except ValueError as exc:
            raise MWRuntimeError(
                f"Path {path} must live under the repository root {self.repo_root} "
                f"for container runtime.",
                suggestion="Move it under the repo or provide an explicit bind mapping.",
            ) from exc
        return str(self.container_repo / rel.as_posix())
