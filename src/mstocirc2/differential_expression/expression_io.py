"""Expression matrix I/O – parsing peptide matrices (DIA‑NN, FragPipe), circRNA reference, design file, and sample‑name matching."""
from __future__ import annotations
import re
import logging
from pathlib import Path
from typing import List, Tuple, Dict, Set
import pandas as pd

log = logging.getLogger(__name__)
_TABLE_ENCODINGS = ("utf-8", "utf-8-sig", "utf-16", "utf-16-le", "utf-16-be")


def _read_table(path: str) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    last_error: Exception | None = None
    for encoding in _TABLE_ENCODINGS:
        try:
            df = pd.read_csv(p, sep="\t", low_memory=False, encoding=encoding)
            if df.shape[1] == 1:
                df = pd.read_csv(p, sep=",", low_memory=False, encoding=encoding)
            return df
        except UnicodeError as exc:
            last_error = exc
            continue
        except Exception:
            raise

    raise UnicodeDecodeError(
        getattr(last_error, "encoding", "utf-8"),
        getattr(last_error, "object", b""),
        getattr(last_error, "start", 0),
        getattr(last_error, "end", 1),
        (
            f"Unable to decode table '{path}'. Supported text encodings include UTF-8, "
            "UTF-8 with BOM, and UTF-16."
        ),
    ) from last_error


def detect_peptide_matrix_format(df: pd.DataFrame) -> str:
    cols = set(df.columns)
    if "Protein.Group" in cols and "Stripped.Sequence" in cols:
        return "diann"
    if ("Peptide Sequence" in cols and "Protein" in cols
            and any(c.endswith(" Intensity") and not c.endswith(" MaxLFQ Intensity")
                    for c in df.columns)):
        return "fragpipe"
    return "unknown"


def normalize_fragpipe_matrix(
    df: pd.DataFrame, design_samples: List[str] | None = None,
) -> pd.DataFrame:
    intensity_pattern = re.compile(r"^(.+) Intensity$")
    maxlfq_suffix = " MaxLFQ Intensity"
    sample_cols: Dict[str, str] = {}
    for c in df.columns:
        if c.endswith(maxlfq_suffix):
            continue
        m = intensity_pattern.match(c)
        if m:
            sample_name = m.group(1).strip()
            if sample_name.endswith(" MaxLFQ"):
                continue
            sample_cols[sample_name] = c
    if not sample_cols:
        raise ValueError(
            "-pm FragPipe peptide matrix has no '<sample> Intensity' columns."
        )
    out_sample_to_col: Dict[str, str]
    if design_samples:
        matched = match_samples_to_columns(list(design_samples), list(sample_cols.keys()))
        out_sample_to_col = {
            design_s: sample_cols[matrix_s] for design_s, matrix_s in matched.items()
        }
    else:
        out_sample_to_col = dict(sample_cols)
    out_cols = {
        "Stripped.Sequence": df["Peptide Sequence"],
        "Protein.Group": df["Protein"].astype(str).str.strip(),
    }
    if "Gene" in df.columns:
        out_cols["Genes"] = df["Gene"]
    elif "Genes" in df.columns:
        out_cols["Genes"] = df["Genes"]
    for sample_name, col_name in out_sample_to_col.items():
        out_cols[sample_name] = pd.to_numeric(df[col_name], errors="coerce")
    out = pd.DataFrame(out_cols)
    log.info(
        f"FragPipe matrix normalized: {len(out)} rows, {len(out_sample_to_col)} samples kept"
    )
    return out


def load_peptide_matrix(
    path: str, design_samples: List[str] | None = None,
) -> Tuple[pd.DataFrame, str]:
    df = _read_table(path)
    fmt = detect_peptide_matrix_format(df)
    if fmt == "diann":
        log.info("Peptide matrix format: DIA-NN pr_matrix")
        return df, fmt
    if fmt == "fragpipe":
        log.info("Peptide matrix format: FragPipe combined_peptide")
        return normalize_fragpipe_matrix(df, design_samples), fmt
    raise ValueError(
        "-pm peptide matrix format is not recognized. Expected a DIA-NN pr_matrix "
        "(columns including Protein.Group and Stripped.Sequence) or a FragPipe "
        f"combined_peptide table (Peptide Sequence, Protein, and '<sample> Intensity'). "
        f"Got columns: {list(df.columns)[:10]}..."
    )


_COUNT_SUFFIX_RE = re.compile(r"\(\d+(?:\|\d+)?\)\s*$")


def parse_circrna_reference(path: str) -> Tuple[pd.DataFrame, Set[str], Set[str]]:
    df = _read_table(path)
    missing = [c for c in ("circ_ORF_ID", "peptide_seq") if c not in df.columns]
    if missing:
        raise ValueError(f"CircRNA_Reference missing column(s): {missing}")
    circ_orf_ids: Set[str] = set(df["circ_ORF_ID"].dropna().astype(str).str.strip())
    circ_orf_ids.discard("")
    circ_ref_peptides: Set[str] = set()
    raw_total = 0
    for entry in df["peptide_seq"].dropna().astype(str):
        for piece in entry.split(";"):
            piece = piece.strip()
            if not piece:
                continue
            raw_total += 1
            piece = _COUNT_SUFFIX_RE.sub("", piece).strip()
            if piece:
                circ_ref_peptides.add(piece.upper())
    if not circ_ref_peptides:
        raise ValueError("CircRNA_Reference has no valid peptide_seq entries")
    return df, circ_ref_peptides, circ_orf_ids


def _design_condition_ok(series: pd.Series) -> bool:
    vals = series.dropna().astype(str).str.strip()
    vals = vals[(vals != "") & (vals.str.lower() != "nan")]
    return vals.nunique() >= 1


def _read_design_table(path: str) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    last_unicode: Exception | None = None
    candidates: list[pd.DataFrame] = []
    for encoding in _TABLE_ENCODINGS:
        try:
            for sep in ("\t", ",", r"\s+"):
                try:
                    if sep == r"\s+":
                        df = pd.read_csv(
                            p,
                            sep=sep,
                            engine="python",
                            encoding=encoding,
                        )
                    else:
                        df = pd.read_csv(p, sep=sep, low_memory=False, encoding=encoding)
                except Exception:
                    continue
                if df.shape[1] < 2:
                    continue
                df.columns = [str(c).strip() for c in df.columns]
                candidates.append(df)
            if candidates:
                break
        except UnicodeError as exc:
            last_unicode = exc
            candidates = []
            continue

    if not candidates and last_unicode is not None:
        raise UnicodeDecodeError(
            getattr(last_unicode, "encoding", "utf-8"),
            getattr(last_unicode, "object", b""),
            getattr(last_unicode, "start", 0),
            getattr(last_unicode, "end", 1),
            (
                f"Unable to decode design file '{path}'. Supported text encodings include UTF-8, "
                "UTF-8 with BOM, and UTF-16."
            ),
        ) from last_unicode

    for df in candidates:
        cols_lower = {str(c).strip().lower(): c for c in df.columns}
        if "sample" in cols_lower and "condition" in cols_lower:
            out = df.rename(
                columns={
                    cols_lower["sample"]: "sample",
                    cols_lower["condition"]: "condition",
                }
            )
            if _design_condition_ok(out["condition"]):
                return out[["sample", "condition"]].copy()

    if candidates:
        df = candidates[0]
        cols_lower = {str(c).strip().lower(): c for c in df.columns}
        if "sample" in cols_lower and "condition" in cols_lower:
            out = df.rename(
                columns={
                    cols_lower["sample"]: "sample",
                    cols_lower["condition"]: "condition",
                }
            )
            return out[["sample", "condition"]].copy()
        raise ValueError(
            "expected columns named 'sample' and 'condition', "
            f"but found columns: {list(df.columns)[:12]}."
        )

    raise ValueError(
        "could not parse design file as tab-, comma-, or whitespace-delimited "
        "table with at least two columns."
    )


def parse_design(path: str) -> pd.DataFrame:
    df = _read_design_table(path)
    missing = [c for c in ("sample", "condition") if c not in df.columns]
    if missing:
        raise ValueError(
            "expected columns named 'sample' and 'condition', "
            f"but missing {missing} (found columns: {list(df.columns)[:12]}). "
            "Use a two-column design table (tab, comma, or whitespace separated), "
            "for example design.from_annotation.tsv from `mstocirc2 quant`."
        )
    df = df[["sample", "condition"]].copy()
    df["sample"] = df["sample"].astype(str).str.strip()
    df["condition"] = df["condition"].astype(str).str.strip()
    empty_cond = (df["condition"] == "") | (df["condition"].str.lower() == "nan")
    if bool(empty_cond.all()):
        raise ValueError(
            "the 'condition' column is empty for all rows. "
            "Check that each line separates sample and condition with a tab, comma, "
            "or spaces (for example: sample_name<TAB>group_label)."
        )
    if bool(empty_cond.any()):
        bad = df.loc[empty_cond, "sample"].astype(str).tolist()[:5]
        raise ValueError(
            f"the 'condition' column is empty for {int(empty_cond.sum())} row(s), "
            f"including sample(s): {bad}."
        )
    df = df.loc[~empty_cond].copy()
    conds = df["condition"].unique().tolist()
    if len(conds) != 2:
        raise ValueError(
            f"exactly 2 conditions are required, got {len(conds)}: {conds}"
        )
    return df


_STRIP_SUFFIXES = {".d", ".raw", ".mzml", ".mzxml", ".wiff", ".tsv", ".csv"}


def _strip_known_suffixes(name: str) -> str:
    base = str(name)
    while True:
        stem, ext = Path(base).stem, Path(base).suffix
        if ext.lower() in _STRIP_SUFFIXES:
            base = stem
            continue
        break
    return base


def _path_without_known_suffix(path_like: str) -> str | None:
    raw = str(path_like)
    suffix = Path(raw).suffix
    if not suffix or suffix.lower() not in _STRIP_SUFFIXES:
        return None
    if len(raw) >= len(suffix) and raw[-len(suffix):].lower() == suffix.lower():
        return raw[: -len(suffix)]
    return None


def match_samples_to_columns(samples: List[str], columns: List[str]) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    unmatched: List[str] = []
    used_cols: Set[str] = set()
    col_list = [str(c) for c in columns]

    for sample in samples:
        s = str(sample).strip()
        hit = None
        if s in col_list and s not in used_cols:
            hit = s
        if hit is None:
            for col in col_list:
                if col in used_cols:
                    continue
                if Path(col).name == s:
                    hit = col
                    break
        if hit is None:
            for col in col_list:
                if col in used_cols:
                    continue
                no_suffix = _path_without_known_suffix(col)
                if no_suffix is not None and no_suffix == s:
                    hit = col
                    break
        if hit is None:
            for col in col_list:
                if col in used_cols:
                    continue
                if _strip_known_suffixes(Path(col).name) == s:
                    hit = col
                    break
        if hit is None:
            unmatched.append(s)
        else:
            mapping[s] = hit
            used_cols.add(hit)

    if unmatched:
        preview = col_list[:8]
        raise ValueError(
            f"-de sample name(s) {unmatched} could not be matched to any column in "
            f"-pm after exact matching, path-stripped matching, suffix-stripped "
            f"matching, and path+suffix-stripped matching. Example -pm columns: {preview}."
        )
    return mapping
