from __future__ import annotations

import argparse
import os
import re
from pathlib import Path
from typing import Any, Sequence

from .exceptions import CLIUsageError

_NONE_TOKENS = {"", "none", "None", "NONE"}
_PYTHON_NAMES = {
    "python",
    "python3",
    "python.exe",
    "python3.exe",
    "pythonw.exe",
}
_START_CODON_RE = re.compile(r"^[ACGT]{3}$")


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


def _path_looks_like_python(value: str) -> bool:
    name = Path(value).name.lower()
    if name in _PYTHON_NAMES:
        return True
    return name.startswith("python")


def _read_text_prefix(path: Path, nbytes: int = 8192) -> str:
    with path.open("rb") as handle:
        chunk = handle.read(nbytes)
    if not chunk:
        return ""
    if b"\x00" in chunk:
        raise CLIUsageError(f"File appears to be binary, not text: '{path}'")
    try:
        return chunk.decode("utf-8")
    except UnicodeDecodeError:
        try:
            return chunk.decode("latin-1")
        except UnicodeDecodeError as exc:
            raise CLIUsageError(
                f"File cannot be decoded as text: '{path}'. {exc}"
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


def _collect_option_values(argv: Sequence[str], option_keys: set[str]) -> list[str]:
    values: list[str] = []
    i = 0
    while i < len(argv):
        token = argv[i]
        key = token.split("=", 1)[0]
        if key in option_keys:
            if "=" in token:
                values.append(token.split("=", 1)[1])
                i += 1
                continue
            i += 1
            while i < len(argv) and not argv[i].startswith("-"):
                values.append(argv[i])
                i += 1
            continue
        i += 1
    return values


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
    dests: dict[int, str] = {}
    nargs_map: dict[int, Any] = {}
    action_types: dict[int, type] = {}
    option_keys_by_action: dict[int, set[str]] = {}
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
        dests[aid] = action.dest
        nargs_map[aid] = action.nargs
        action_types[aid] = type(action)
        option_keys_by_action.setdefault(aid, set()).update(action.option_strings)

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
        label = labels[aid]
        hint = ""
        if dests[aid] in {"canonical_protein", "linear_protein"}:
            values = _collect_option_values(argv, option_keys_by_action[aid])
            if any(_path_looks_like_python(v) for v in values):
                hint = (
                    " One of the values looks like a Python interpreter; "
                    "use '-dp/--deepcip-python' for DeepCIP Python or '-py/--python-bin' for FragPipe Python."
                )
        if nargs in ("+", "*"):
            raise CLIUsageError(
                f"Option {label} was provided {count} times. "
                f"List multiple values after a single {label}, for example: "
                f"{label.split(' / ')[0]} a.fasta b.fasta. "
                f"Repeating the flag overwrites earlier values.{hint}"
            )
        raise CLIUsageError(
            f"Option {label} was provided {count} times. Provide it only once "
            f"(later values silently overwrite earlier ones).{hint}"
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


def _validate_orf_input_mode(args: argparse.Namespace) -> None:
    if not _is_none(getattr(args, "circ_seq_file", None)):
        return
    if _is_none(getattr(args, "dna_seq", None)) or _is_none(getattr(args, "dna_gff", None)):
        raise CLIUsageError(
            "Genome-based ORF prediction requires '-ds/--dna-seq' and '-dg/--dna-gff' "
            "when '-cs/--circ-seq' is not provided."
        )
    source_attrs = ("circ_info", "find_circ", "circexplorer", "CIRI")
    if all(_is_none(getattr(args, name, None)) for name in source_attrs):
        raise CLIUsageError(
            "Provide at least one circRNA source via '-ci/--circ-info', '-fc/--find-circ', "
            "'-ce/--circexplorer', or '-ciri/--ciri'."
        )


def _validate_start_codon(value: Any) -> None:
    if value is None:
        return
    raw = str(value)
    if re.search(r"\s", raw):
        raise CLIUsageError(
            "Invalid '-sc/--start-codon' value. Use comma-separated DNA codons, "
            "for example: ATG,TTG,CTG,GTG."
        )
    codons = [c.strip().upper() for c in raw.split(",") if c.strip()]
    if not codons:
        raise CLIUsageError("At least one valid start codon must be provided via '-sc/--start-codon'.")
    for codon in codons:
        if not _START_CODON_RE.fullmatch(codon):
            raise CLIUsageError(
                f"Invalid start codon '{codon}'. Start codons must be 3-base DNA codons using A/C/G/T."
            )


def _validate_dea_strategy(value: Any) -> None:
    if _is_none(value):
        return
    from ..DEA.workflow_strategy import STRATEGY_REGISTRY

    name = str(value)
    if name not in STRATEGY_REGISTRY:
        valid = ", ".join(sorted(STRATEGY_REGISTRY.keys()))
        raise CLIUsageError(f"Unknown '-st/--strategy' value '{name}'. Valid values: {valid}")


def _validate_dea_organism(value: Any) -> None:
    if _is_none(value):
        return
    from ..DEA.kegg_registry import normalize_kegg_organism

    try:
        normalize_kegg_organism(str(value))
    except ValueError as exc:
        raise CLIUsageError(str(exc)) from exc


def _validate_eval_input_combo(args: argparse.Namespace) -> None:
    has_fi = not _is_none(getattr(args, "circ_file_input", None))
    has_co = not _is_none(getattr(args, "circ_orf", None))
    if not has_fi and not has_co:
        raise CLIUsageError(
            "Eval requires either '-fi/--circ-file-input' or '-co/--circ-orf'."
        )


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
        _validate_orf_input_mode(args)
        _validate_start_codon(getattr(args, "start_codon", None))
        min_orf_len = getattr(args, "min_orf_len", None)
        if min_orf_len is not None and int(min_orf_len) < 1:
            raise CLIUsageError("'-ml/--min-orf-len' must be >= 1.")
        flank_aa_len = getattr(args, "flank_aa_len", None)
        if flank_aa_len is not None and int(flank_aa_len) < 0:
            raise CLIUsageError("'-fl/--flank-aa-len' must be >= 0.")

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

    if command == "quant":
        _validate_dea_strategy(getattr(args, "strategy", None))
        _validate_dea_organism(getattr(args, "organism", None))

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
        _validate_eval_input_combo(args)

    if command == "dea":
        require_readable_text_file(getattr(args, "peptide_matrix"), "-pm/--peptide-matrix")
        require_readable_text_file(getattr(args, "design"), "-de/--design")
        _validate_optional_file(getattr(args, "circrna_reference", None), "-cr/--circrna-reference")
        _validate_dea_strategy(getattr(args, "strategy", None))
        _validate_dea_organism(getattr(args, "organism", None))
