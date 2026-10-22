"""Wrappers for external command discovery and execution."""

from __future__ import annotations

import shutil

from .exceptions import DependencyError


def require_command(executable: str) -> str:
    resolved = shutil.which(executable)
    if resolved is None:
        raise DependencyError(f"Required executable not found on PATH: '{executable}'.")
    return resolved
