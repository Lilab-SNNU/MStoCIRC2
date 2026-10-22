from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Any, Iterable, Sequence

from .exceptions import CLIUsageError

_NONE_TOKENS = {"", "none", "None", "NONE"}
_PYTHON_NAMES = {
    "python",
    "python3",
    "python.exe",
    "python3.exe",
    "pythonw.exe",
}


def _is_none(value: Any) -> bool:
    if value is None:
        return True
    if isinstance(value, str) and value.strip() in _NONE_TOKENS:
        return True
    return False


def _as_paths(value: Any) -> list[Path]:
    if _is_none(value):
        return []
    if isinstance(value, (list, tuple)):
        paths: list[Path] = []
        for item in value:
            paths.extend(_as_paths(item))
        return paths
    return [Path(str(value)).expanduser()]


def _option_label(dest: str, option_strings: Sequence[str] | None = None) -> str:
    if option_strings:
        return " / ".join(option_strings)
    return f"--{dest.replace('_', '-')}"


def _looks_like_python_interpreter(path: Path) -> bool:
    name = path.name.lower()
    if name in _PYTHON_NAMES:
        return True
    if name.startswith("python") and path.is_file():
        return True
    return False


def _read_text_prefix(path: Path, nbytes: int = 8192) -> str:
    with path.open("rb") as handle:
        chunk = handle.read(nbytes)
    if not chunk:
        return ""
    if b"\x00" in chunk:
        raise CLIUsageError(f"File appears to be binary, not text: '{path}'")
    try:
        return chunk.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise CLIUsageError(
            f"File is not valid UTF-8 text: '{path}'. {exc}"
        ) from exc


def require_existing_file(path_value: Any, label: str) -> Path:
    path = Path(str(path_value)).expanduser()
    if not path.exists():
        raise CLIUsageError(f"{label} not found: '{path}'")
    if not path.is_file():
        raise CLIUsageError(f"{label} is not a file: '{path}'")
    return path.resolve()


def require_existing_path(path_value: Any, label: str) -> Path:
    path = Path(str(path_value)).expanduser()
    if not path.exists():
        raise CLIUsageError(f"{label} not found: '{path}'")
    return path.resolve()


def require_existing_dir(path_value: Any, label: str) -> Path:
    path = Path(str(path_value)).expanduser()
    if not path.exists():
        raise CLIUsageError(f"{label} not found: '{path}'")
    if not path.is_dir():
        raise CLIUsageError(f"{label} is not a directory: '{path}'")
    return path.resolve()


def require_executable(path_value: Any, label: str) -> Path:
    path = require_existing_file(path_value, label)
    if os.name != "nt" and not os.access(path, os.X_OK):
        raise CLIUsageError(f"{label} is not executable: '{path}'")
    return path


def require_readable_text_file(path_value: Any, label: str) -> Path:
    path = require_existing_file(path_value, label)
    _read_text_prefix(path)
    return path


def require_fasta_file(path_value: Any, label: str) -> Path:
    path = require_existing_file(path_value, label)
    if _looks_like_python_interpreter(path):
        raise CLIUsageError(
            f"{label} looks like a Python interpreter: '{path}'. "
            "Use '-dp/--deepcip-python' for DeepCIP Python or '-py/--python-bin' for FragPipe Python. "
            f"{label} must point to a FASTA file."
        )
    text = _read_text_prefix(path)
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if not stripped.startswith(">"):
            raise CLIUsageError(
                f"{label} does not look like a FASTA file "
                f"(first non-empty line should start with '>'): '{path}'"
            )
        return path.resolve()
    raise CLIUsageError(f"{label} is empty: '{path}'")


def _subparser_for_command(root_parser: argparse.ArgumentParser, command: str) -> argparse.ArgumentParser | None:
    for action in root_parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            return action.choices.get(command)
    return None


def check_duplicate_cli_options(argv: Sequence[str], root_parser: argparse.ArgumentParser) -> None:
    command = None
    for token in argv:
        if not token.startswith("-"):
            command = token
            break
    if not command:
        return
    subparser = _subparser_for_command(root_parser, command)
    if subparser is None:
        return

    option_to_action: dict[str, argparse.Action] = {}
    for action in subparser._actions:
        if not action.option_strings:
            continue
        if isinstance(action, (argparse._HelpAction,)):
            continue
        for opt in action.option_strings:
            option_to_action[opt] = action

    counts: dict[int, int] = {}
    labels: dict[int, str] = {}
    nargs_map: dict[int, Any] = {}
    action_types: dict[int, type] = {}
    seen_command = False
    for token in argv:
        if not seen_command:
            if token == command:
                seen_command = True
            continue
        if token == "--":
            break
        if not token.startswith("-"):
            continue
        key = token.split("=", 1)[0]
        action = option_to_action.get(key)
        if action is None:
            continue
        aid = id(action)
        counts[aid] = counts.get(aid, 0) + 1
        labels[aid] = _option_label(action.dest, action.option_strings)
        nargs_map[aid] = action.nargs
        action_types[aid] = type(action)

    for aid, count in counts.items():
        if count <= 1:
            continue
        nargs = nargs_map[aid]
        action_type = action_types[aid]
        if action_type in (
            argparse._StoreTrueAction,
            argparse._StoreFalseAction,
            argparse._CountAction,
            argparse._AppendAction,
            argparse._AppendConstAction,
        ):
            continue
        if nargs in ("+", "*"):
            raise CLIUsageError(
                f"Option {labels[aid]} was provided {count} times. "
                f"List multiple values after a single {labels[aid]}, for example: "
                f"{labels[aid].split(' / ')[0]} a.fasta b.fasta. "
                f"Repeating the flag overwrites earlier values."
            )
        raise CLIUsageError(
            f"Option {labels[aid]} was provided {count} times. Provide it only once "
            f"(later values silently overwrite earlier ones)."
        )


def _validate_optional_file(value: Any, label: str, *, fasta: bool = False) -> None:
    if _is_none(value):
        return
    if fasta:
        require_fasta_file(value, label)
    else:
        require_readable_text_file(value, label)


def _validate_optional_dir(value: Any, label: str) -> None:
    if _is_none(value):
        return
    require_existing_dir(value, label)


def _validate_deepcip_path(value: Any, label: str) -> None:
    if _is_none(value):
        return
    path = require_existing_dir(value, label)
    script = path / "DeepCIP.py"
    if not script.is_file():
        raise CLIUsageError(
            f"{label} does not contain DeepCIP.py: '{path}'. "
            "Point -dc/--deepcip-path to the DeepCIP installation directory."
        )


def _validate_optional_executable(value: Any, label: str) -> None:
    if _is_none(value):
        return
    require_executable(value, label)


def _validate_fasta_list(value: Any, label: str) -> None:
    paths = _as_paths(value)
    if not paths:
        raise CLIUsageError(f"{label} is required.")
    for path in paths:
        require_fasta_file(path, label)


def preflight_cli_args(args: argparse.Namespace) -> None:
    command = getattr(args, "command", None)
    if not command:
        return

    if command in {"orf", "nonquant", "quant"}:
        _validate_fasta_list(getattr(args, "canonical_protein", None), "-cp/--canonical-protein")
        _validate_optional_file(getattr(args, "circ_seq_file", None), "-cs/--circ-seq", fasta=True)
        _validate_optional_file(getattr(args, "circ_info", None), "-ci/--circ-info")
        _validate_optional_file(getattr(args, "dna_seq", None), "-ds/--dna-seq", fasta=True)
        _validate_optional_file(getattr(args, "dna_gff", None), "-dg/--dna-gff")
        _validate_optional_file(getattr(args, "find_circ", None), "-fc/--find-circ")
        _validate_optional_file(getattr(args, "circexplorer", None), "-ce/--circexplorer")
        _validate_optional_file(getattr(args, "CIRI", None), "-ciri/--ciri")

    if command == "search":
        _validate_optional_file(getattr(args, "circ_orf", None), "-co/--circ-orf", fasta=True)
        require_fasta_file(getattr(args, "linear_protein"), "-cp/--canonical-protein")
        require_readable_text_file(getattr(args, "manifest"), "-mf/--manifest")
        require_executable(getattr(args, "fragpipe_bin"), "-fb/--fragpipe-bin")
        _validate_optional_file(getattr(args, "workflow", None), "-wf/--workflow")
        _validate_optional_dir(getattr(args, "tools_folder", None), "-td/--tools-dir")
        _validate_optional_executable(getattr(args, "diann_bin", None), "-db/--diann-bin")
        _validate_optional_executable(getattr(args, "python_bin", None), "-py/--python-bin")

    if command in {"nonquant", "quant"}:
        require_readable_text_file(getattr(args, "manifest"), "-mf/--manifest")
        require_executable(getattr(args, "fragpipe_bin"), "-fb/--fragpipe-bin")
        _validate_optional_file(getattr(args, "workflow", None), "-wf/--workflow")
        _validate_optional_dir(getattr(args, "tools_folder", None), "-td/--tools-dir")
        _validate_optional_executable(getattr(args, "diann_bin", None), "-db/--diann-bin")
        _validate_optional_executable(getattr(args, "python_bin", None), "-py/--python-bin")
        _validate_deepcip_path(getattr(args, "deepcip_path", None), "-dc/--deepcip-path")
        _validate_optional_executable(getattr(args, "deepcip_python", None), "-dp/--deepcip-python")

    if command == "eval":
        require_fasta_file(getattr(args, "canonical_protein"), "-cp/--canonical-protein")
        require_existing_path(getattr(args, "path_ms_input"), "-ms/--ms-input")
        ms_path = Path(str(getattr(args, "path_ms_input"))).expanduser()
        if ms_path.is_file():
            require_readable_text_file(ms_path, "-ms/--ms-input")
        _validate_optional_dir(getattr(args, "circ_file_input", None), "-fi/--circ-file-input")
        _validate_optional_file(getattr(args, "circ_seq", None), "-cs/--circ-seq", fasta=True)
        _validate_optional_file(getattr(args, "circ_info", None), "-ci/--circ-info")
        _validate_optional_file(getattr(args, "circ_orf", None), "-co/--circ-orf", fasta=True)
        _validate_optional_file(getattr(args, "circ_mapping", None), "-cm/--circ-mapping")
        _validate_deepcip_path(getattr(args, "deepcip_path", None), "-dc/--deepcip-path")
        _validate_optional_executable(getattr(args, "deepcip_python", None), "-dp/--deepcip-python")

    if command == "dea":
        require_readable_text_file(getattr(args, "peptide_matrix"), "-pm/--peptide-matrix")
        require_readable_text_file(getattr(args, "design"), "-de/--design")
        _validate_optional_file(getattr(args, "circrna_reference", None), "-cr/--circrna-reference")
