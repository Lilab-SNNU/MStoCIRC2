"""Expression matrix I/O – parsing peptide matrices (DIA‑NN, FragPipe), circRNA reference, design file, and sample‑name matching."""
from __future__ import annotations
import re
import logging
from pathlib import Path
from typing import List, Tuple, Dict, Set
import pandas as pd

from .homology_refining import normalize_protein_id

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
    missing = [c for c in ("search_ORF_IDs", "peptide_seq") if c not in df.columns]
    if missing:
        raise ValueError(f"CircRNA_Reference missing column(s): {missing}")
    search_orf_ids: Set[str] = set()
    for entry in df["search_ORF_IDs"].dropna().astype(str):
        for piece in entry.split(";"):
            piece = normalize_protein_id(piece.strip())
            if piece:
                search_orf_ids.add(piece)
    if not search_orf_ids:
        raise ValueError("CircRNA_Reference has no valid search_ORF_IDs entries")
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
    return df, circ_ref_peptides, search_orf_ids


def _design_condition_ok(series: pd.Series) -> bool:
    vals = series.dropna().astype(str).str.strip()
    vals = vals[(vals != "") & (vals.str.lower() != "nan")]
    return vals.nunique() >= 1


def _as_design_frame(df: pd.DataFrame) -> pd.DataFrame | None:
    if df.shape[1] < 2:
        return None
    cols_lower = {str(c).strip().lower(): c for c in df.columns}
    if "sample" in cols_lower and "condition" in cols_lower:
        out = df.rename(
            columns={
                cols_lower["sample"]: "sample",
                cols_lower["condition"]: "condition",
            }
        )[["sample", "condition"]].copy()
    else:
        out = df.iloc[:, :2].copy()
        out.columns = ["sample", "condition"]
        if len(out) > 0:
            first_sample = str(out.iloc[0, 0]).strip().lower()
            first_cond = str(out.iloc[0, 1]).strip().lower()
            if first_sample == "sample" and first_cond == "condition":
                out = out.iloc[1:].reset_index(drop=True)
    return out


def _read_design_table(path: str) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    last_unicode: Exception | None = None
    named_ok: list[pd.DataFrame] = []
    positional_ok: list[pd.DataFrame] = []
    for encoding in _TABLE_ENCODINGS:
        try:
            for sep in ("\t", ",", r"\s+"):
                for use_header in (True, False):
                    try:
                        kwargs: dict = {
                            "sep": sep,
                            "encoding": encoding,
                            "header": 0 if use_header else None,
                        }
                        if sep == r"\s+":
                            kwargs["engine"] = "python"
                        else:
                            kwargs["low_memory"] = False
                        df = pd.read_csv(p, **kwargs)
                    except Exception:
                        continue
                    if df.shape[1] < 2:
                        continue
                    if use_header:
                        df.columns = [str(c).strip() for c in df.columns]
                    coerced = _as_design_frame(df)
                    if coerced is None or coerced.empty:
                        continue
                    cols_lower = {str(c).strip().lower(): c for c in df.columns}
                    if (
                        use_header
                        and "sample" in cols_lower
                        and "condition" in cols_lower
                        and _design_condition_ok(coerced["condition"])
                    ):
                        named_ok.append(coerced)
                    elif (not use_header) and _design_condition_ok(coerced["condition"]):
                        positional_ok.append(coerced)
            if named_ok or positional_ok:
                break
        except UnicodeError as exc:
            last_unicode = exc
            named_ok = []
            positional_ok = []
            continue

    if not named_ok and not positional_ok and last_unicode is not None:
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

    if named_ok:
        return named_ok[0]
    if positional_ok:
        return positional_ok[0]
    raise ValueError(
        "could not parse design file as a two-column tab-, comma-, or whitespace-delimited table."
    )


def parse_design(path: str) -> pd.DataFrame:
    df = _read_design_table(path)
    df = df[["sample", "condition"]].copy()
    df["sample"] = df["sample"].astype(str).str.strip()
    df["condition"] = df["condition"].astype(str).str.strip()
    empty_cond = (df["condition"] == "") | (df["condition"].str.lower() == "nan")
    if bool(empty_cond.all()):
        raise ValueError("the condition values are empty for all rows.")
    if bool(empty_cond.any()):
        bad = df.loc[empty_cond, "sample"].astype(str).tolist()[:5]
        raise ValueError(
            f"the condition value is empty for {int(empty_cond.sum())} row(s), "
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
_UNCALIBRATED_TOKEN = "_uncalibrated"


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


def _strip_uncalibrated_token(name: str) -> str | None:
    raw = str(name)
    marker = _UNCALIBRATED_TOKEN
    if raw.lower().endswith(marker):
        return raw[: -len(marker)]
    return None


def _de_sample_forms(sample: str) -> List[str]:
    s = str(sample).strip()
    forms: List[str] = []
    for value in (
        s,
        Path(s).name,
        _path_without_known_suffix(s),
        _strip_known_suffixes(Path(s).name),
    ):
        if value is None:
            continue
        text = str(value).strip()
        if text and text not in forms:
            forms.append(text)
    return forms


def _find_pm_column(de_form: str, col_list: List[str], used_cols: Set[str]) -> str | None:
    if de_form in col_list and de_form not in used_cols:
        return de_form
    for col in col_list:
        if col in used_cols:
            continue
        if Path(col).name == de_form:
            return col
    for col in col_list:
        if col in used_cols:
            continue
        no_suffix = _path_without_known_suffix(col)
        if no_suffix is not None and no_suffix == de_form:
            return col
    for col in col_list:
        if col in used_cols:
            continue
        if _strip_known_suffixes(Path(col).name) == de_form:
            return col
    for col in col_list:
        if col in used_cols:
            continue
        stem = _strip_known_suffixes(Path(col).name)
        cleaned = _strip_uncalibrated_token(stem)
        if cleaned is not None and cleaned == de_form:
            return col
    return None


def match_samples_to_columns(samples: List[str], columns: List[str]) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    unmatched: List[str] = []
    used_cols: Set[str] = set()
    col_list = [str(c) for c in columns]

    for sample in samples:
        s = str(sample).strip()
        hit = None
        for de_form in _de_sample_forms(s):
            hit = _find_pm_column(de_form, col_list, used_cols)
            if hit is not None:
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
            f"-pm after exact matching, path/suffix normalization, and "
            f"_uncalibrated stripping. Example -pm columns: {preview}."
        )
    return mapping
