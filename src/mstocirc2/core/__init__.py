"""Shared infrastructure for the MStoCIRC2 CLI."""

from .exceptions import CLIUsageError, DependencyError, MStoCIRC2Error
from .logging_utils import TerminalProgressBar, configure_logging, module_file_log, module_log_context
from .paths import (
    ensure_directory,
    make_default_module_output_dir,
)
from .preflight import check_duplicate_cli_options, preflight_cli_args
from .process import require_command

__all__ = [
    "CLIUsageError",
    "DependencyError",
    "MStoCIRC2Error",
    "TerminalProgressBar",
    "configure_logging",
    "module_file_log",
    "module_log_context",
    "ensure_directory",
    "make_default_module_output_dir",
    "check_duplicate_cli_options",
    "preflight_cli_args",
    "require_command",
]
