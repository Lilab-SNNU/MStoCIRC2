"""circRNA‑mRNA homologous peptide filter – removes peptides that appear in both circRNA and linear proteomes."""
from __future__ import annotations
import re
from typing import List, Tuple, Set, Dict
import pandas as pd

STRIPPED_SEQ_CANDIDATES = ["Stripped.Sequence", "Sequence", "peptide_seq",
                            "PeptideSeq", "Peptide"]
PROTEIN_ASSIGN_CANDIDATES = ["Protein.Group", "Protein.Ids",
                              "Leading razor protein", "Proteins", "Protein IDs"]

_PE_SUFFIX_RE = re.compile(r"\s*PE=\d+\s*$", re.IGNORECASE)


def _first_present(df: pd.DataFrame, candidates: List[str]) -> str:
    for c in candidates:
        if c in df.columns:
            return c
    raise ValueError(f"Missing columns: {candidates}. Got: {list(df.columns)[:15]}")


def normalize_protein_id(pid: str) -> str:
    return _PE_SUFFIX_RE.sub("", str(pid).strip()).strip()


def is_circrna_search_id(pid: str, search_orf_ids: Set[str]) -> bool:
    if not pid or not search_orf_ids:
        return False
    return normalize_protein_id(pid) in search_orf_ids


def filter_peptides(
    peptide_df: pd.DataFrame,
    circ_ref_peptides: Set[str],
    search_orf_ids: Set[str],
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, int]]:
    seq_col = _first_present(peptide_df, STRIPPED_SEQ_CANDIDATES)
    prot_col = _first_present(peptide_df, PROTEIN_ASSIGN_CANDIDATES)
    input_n = len(peptide_df)
    report_rows: List[Dict[str, str]] = []
    keep = []
    circ_kept = 0
    for seq, prot in zip(peptide_df[seq_col].astype(str),
                         peptide_df[prot_col].astype(str)):
        seq_upper = seq.strip().upper()
        pid = prot.strip()
        is_circ = is_circrna_search_id(pid, search_orf_ids)
        if is_circ:
            if seq_upper in circ_ref_peptides:
                keep.append(True)
                circ_kept += 1
            else:
                keep.append(False)
                report_rows.append({
                    "peptide_seq": seq_upper,
                    "protein_group": pid,
                    "removed_reason": "mrna_homology",
                })
        else:
            keep.append(True)
    out = peptide_df.loc[keep].copy().reset_index(drop=True)
    summary = {
        "input_peptide_count": int(input_n),
        "peptides_removed_mrna_homology": int(sum(1 for k in keep if not k)),
        "peptides_retained": int(len(out)),
        "circrna_peptides_retained": int(circ_kept),
    }
    report = pd.DataFrame(report_rows, columns=[
        "peptide_seq", "protein_group", "removed_reason"])
    return out, report, summary
