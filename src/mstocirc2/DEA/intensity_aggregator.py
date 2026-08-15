"""Peptide‑to‑protein intensity aggregation via directLFQ or median summarisation."""
from __future__ import annotations
import logging
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Set, Tuple
import numpy as np
import pandas as pd
from .homology_refining import (
    STRIPPED_SEQ_CANDIDATES, PROTEIN_ASSIGN_CANDIDATES,
    _first_present, is_circrna_search_id,
)

log = logging.getLogger(__name__)


def _build_long(peptide_df: pd.DataFrame, sample_columns: Dict[str, str]) -> pd.DataFrame:
    seq_col = _first_present(peptide_df, STRIPPED_SEQ_CANDIDATES)
    prot_col = _first_present(peptide_df, PROTEIN_ASSIGN_CANDIDATES)
    reverse_map = {v: k for k, v in sample_columns.items()}
    sample_cols_matrix = list(sample_columns.values())
    df = peptide_df[[seq_col, prot_col] + sample_cols_matrix].copy()
    df = df.rename(columns={seq_col: "peptide", prot_col: "protein_id"})
    df = df.rename(columns=reverse_map)
    df["protein_id"] = df["protein_id"].astype(str).str.strip()
    df = df[df["protein_id"] != ""]
    design_samples = list(sample_columns.keys())
    for s in design_samples:
        df[s] = pd.to_numeric(df[s], errors="coerce")
    return df[["protein_id", "peptide"] + design_samples]


def _median_rollup(long_df: pd.DataFrame, design_samples: List[str]) -> pd.DataFrame:
    logdf = long_df.copy()
    for s in design_samples:
        vals = logdf[s].to_numpy(dtype=float, copy=True)
        with np.errstate(divide="ignore", invalid="ignore"):
            vals = np.where(vals > 0, np.log2(vals), np.nan)
        logdf[s] = vals
    agg = logdf.groupby("protein_id")[design_samples].median()
    return np.power(2.0, agg)


def _run_directlfq(long_df: pd.DataFrame, design_samples: List[str]) -> pd.DataFrame:
    try:
        import directlfq  # noqa: F401
    except ImportError:
        raise ImportError("'directlfq' package required. Install: pip install directlfq")
    agg_df = long_df.groupby(["protein_id", "peptide"])[design_samples].sum().reset_index()
    wide = agg_df.set_index(["protein_id", "peptide"])[design_samples]
    wide.index.names = ["protein", "ion"]
    wide = wide.replace(0, np.nan)
    wide = np.log2(wide)
    log.info(
        "directLFQ: prepared wide matrix rows=%d, cols=%d",
        len(wide),
        len(wide.columns),
    )
    with tempfile.TemporaryDirectory(prefix="mstocirc2_directlfq_") as tmpdir:
        tmpdir_path = Path(tmpdir)
        input_path = tmpdir_path / "input.pkl"
        output_path = tmpdir_path / "output.pkl"
        wide.to_pickle(input_path)
        proc = subprocess.run(
            [
                sys.executable,
                "-m",
                "mstocirc2.DEA.directlfq_bridge",
                str(input_path),
                str(output_path),
            ],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            check=False,
        )
        log.info(
            "directLFQ bridge exited with code=%d, stdout_chars=%d, stderr_chars=%d",
            proc.returncode,
            len(proc.stdout or ""),
            len(proc.stderr or ""),
        )
        if proc.returncode != 0:
            raise RuntimeError(
                "directLFQ roll-up failed. "
                + (proc.stderr.strip() or proc.stdout.strip() or "directLFQ subprocess failed.")
            )
        if not output_path.exists():
            raise RuntimeError("directLFQ subprocess finished without writing an output matrix.")
        protein_df = pd.read_pickle(output_path)

    protein_df = protein_df.set_index("protein")
    protein_df.index.name = "protein_id"
    return protein_df[design_samples]


def filter_proteins_by_group_detection(
    protein_matrix: pd.DataFrame,
    design: pd.DataFrame,
) -> Tuple[pd.DataFrame, Dict[str, int]]:
    matched = [
        s for s in design["sample"].astype(str) if s in protein_matrix.columns
    ]
    if not matched:
        raise ValueError("No design samples match protein matrix columns for detection filtering.")
    keep_mask = pd.Series(False, index=protein_matrix.index)
    for _, grp in design.groupby("condition", sort=False):
        samples = [s for s in grp["sample"].astype(str) if s in protein_matrix.columns]
        if not samples:
            continue
        n = len(samples)
        vals = protein_matrix[samples].apply(pd.to_numeric, errors="coerce")
        detected = vals.notna() & (vals > 0)
        keep_mask = keep_mask | (detected.sum(axis=1) > (n / 2.0))
    filtered = protein_matrix.loc[keep_mask].copy()
    summary = {
        "input_protein_count": int(len(protein_matrix)),
        "proteins_retained": int(len(filtered)),
        "proteins_removed_low_detection": int(len(protein_matrix) - len(filtered)),
    }
    return filtered, summary


def rollup_peptides_to_proteins(
    peptide_df: pd.DataFrame,
    sample_columns: Dict[str, str],
    expression_matrix_type: str,
    search_orf_ids: Set[str],
    pre_filter_counts: pd.Series | None = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    design_samples = list(sample_columns.keys())
    long_df = _build_long(peptide_df, sample_columns)
    log.info(
        "Roll-up builder: long_df rows=%d, proteins=%d, peptides=%d",
        len(long_df),
        long_df["protein_id"].nunique(),
        long_df["peptide"].nunique(),
    )
    retained_counts = long_df.groupby("protein_id")["peptide"].nunique().rename(
        "retained_peptide_count")
    if expression_matrix_type.lower().startswith("directlfq"):
        matrix = _run_directlfq(long_df, design_samples)
    else:
        matrix = _median_rollup(long_df, design_samples)
    if pre_filter_counts is not None:
        input_counts = pre_filter_counts.rename("input_peptide_count")
        summary = pd.concat([input_counts, retained_counts], axis=1)
        summary["input_peptide_count"] = summary["input_peptide_count"].fillna(
            summary["retained_peptide_count"]).astype(int)
    else:
        summary = retained_counts.to_frame()
        summary["input_peptide_count"] = summary["retained_peptide_count"]
    summary["removed_peptides"] = (
        summary["input_peptide_count"] - summary["retained_peptide_count"])
    summary["is_circrna"] = summary.index.to_series().apply(
        lambda pid: is_circrna_search_id(pid, search_orf_ids))
    summary = summary.reset_index()
    summary = summary[["protein_id", "input_peptide_count",
                       "retained_peptide_count", "removed_peptides", "is_circrna"]]
    matrix = matrix.loc[matrix.notna().any(axis=1)]
    summary = summary[summary["protein_id"].isin(matrix.index)].reset_index(drop=True)
    return matrix, summary
