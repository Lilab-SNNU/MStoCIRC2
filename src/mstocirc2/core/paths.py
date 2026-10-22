"""Path helpers used across pipeline stages."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path


def ensure_directory(path: str | Path) -> Path:
    directory = Path(path).expanduser()
    directory.mkdir(parents=True, exist_ok=True)
    return directory


def make_default_module_output_dir(
    module_name: str,
    parent: str | Path | None = None,
) -> Path:
    root = Path(parent).expanduser() if parent else Path.cwd()
    stamp = datetime.now().strftime("%y-%m-%d")
    counter = 1
    while True:
        candidate = root / f"MStoCIRC2_{module_name}_{stamp}.{counter}"
        try:
            candidate.mkdir(parents=True, exist_ok=False)
            return candidate
        except FileExistsError:
            counter += 1
