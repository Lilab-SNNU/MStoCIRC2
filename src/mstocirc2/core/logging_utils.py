"""Central logging configuration."""

from __future__ import annotations

import logging
import shutil
import sys
from contextlib import contextmanager
from contextvars import ContextVar
from datetime import datetime
from pathlib import Path
from typing import Iterator

_MODULE_TAG: ContextVar[str | None] = ContextVar("mstocirc2_module_tag", default=None)

_LOGGER_MODULE_TAGS = {
    "ORF": "ORF",
    "Search": "Search",
    "Eval": "Eval",
    "DEA": "DEA",
    "circrna_dea": "DEA",
}

_PACKAGE_MODULE_PREFIXES = (
    ("mstocirc2.ORF", "ORF"),
    ("mstocirc2.Search", "Search"),
    ("mstocirc2.Eval", "Eval"),
    ("mstocirc2.DEA", "DEA"),
)


def _resolve_module_tag(logger_name: str) -> str | None:
    tag = _MODULE_TAG.get()
    if tag:
        return tag
    if logger_name in _LOGGER_MODULE_TAGS:
        return _LOGGER_MODULE_TAGS[logger_name]
    for prefix, mapped in _PACKAGE_MODULE_PREFIXES:
        if logger_name == prefix or logger_name.startswith(prefix + "."):
            return mapped
    return None


@contextmanager
def module_log_context(module_name: str) -> Iterator[None]:
    token = _MODULE_TAG.set(module_name)
    try:
        yield
    finally:
        _MODULE_TAG.reset(token)


class TerminalProgressBar:
    """Minimal terminal progress bar that does not flood the log history."""

    def __init__(self, total: int, prefix: str, width: int = 28) -> None:
        self.total = max(int(total), 0)
        self.prefix = prefix
        self.width = max(int(width), 10)
        self._enabled = self.total > 0 and sys.stderr.isatty()
        self._closed = False

    def update(self, current: int) -> None:
        if not self._enabled or self._closed:
            return
        current = max(0, min(int(current), self.total))
        columns = shutil.get_terminal_size((100, 20)).columns
        bar_width = min(self.width, max(10, columns - len(self.prefix) - 28))
        filled = int(bar_width * current / self.total) if self.total else bar_width
        bar = "#" * filled + "-" * (bar_width - filled)
        percent = 100.0 * current / self.total if self.total else 100.0
        sys.stderr.write(f"\r{self.prefix} [{bar}] {current}/{self.total} ({percent:5.1f}%)")
        sys.stderr.flush()

    def close(self) -> None:
        if not self._enabled or self._closed:
            return
        self._closed = True
        sys.stderr.write("\n")
        sys.stderr.flush()


class CompactCliFormatter(logging.Formatter):
    """Compact terminal formatter for end users."""

    def format(self, record: logging.LogRecord) -> str:
        message = record.getMessage()
        tag = _resolve_module_tag(record.name)
        if record.levelno >= logging.ERROR:
            prefix = "[ERROR]"
        elif record.levelno >= logging.WARNING:
            prefix = "[WARN]"
        elif record.levelno == logging.DEBUG:
            prefix = "[DEBUG]"
        else:
            prefix = f"[{tag}]" if tag else None
            if prefix:
                return f"{prefix} {message}"
            return message
        if tag and prefix != "[DEBUG]":
            return f"[{tag}] {prefix} {message}"
        if prefix == "[DEBUG]":
            return f"{prefix} {record.name}: {message}"
        return f"{prefix} {message}"


def configure_logging(verbose: int = 0, quiet: bool = False) -> None:
    """Configure a single global logging policy for the CLI."""
    if quiet:
        level = logging.ERROR
    elif verbose >= 2:
        level = logging.DEBUG
    elif verbose == 1:
        level = logging.INFO
    else:
        level = logging.INFO

    handler = logging.StreamHandler()
    handler.setFormatter(CompactCliFormatter())
    logging.basicConfig(level=level, handlers=[handler], force=True)


@contextmanager
def module_file_log(output_dir: str | Path) -> Iterator[Path]:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    log_path = out / f"log_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.txt"
    handler = logging.FileHandler(log_path, encoding="utf-8")
    handler.setFormatter(CompactCliFormatter())
    handler.setLevel(logging.NOTSET)
    root = logging.getLogger()
    root.addHandler(handler)
    try:
        yield log_path
    finally:
        root.removeHandler(handler)
        handler.close()
