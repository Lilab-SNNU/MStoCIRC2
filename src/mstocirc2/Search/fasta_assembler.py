"""FASTA assembly – reading, decoy generation (Met‑retention), header tagging, and merging."""
from __future__ import annotations
import logging
from pathlib import Path
from typing import List, Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .header_tagger import tag_header

log = logging.getLogger(__name__)

_DEFAULT_CONTAM_FASTA = Path(__file__).resolve().parent / "contam.fasta"


def _record_title(record: SeqRecord) -> str:
    rid = (record.id or "").strip()
    desc = (record.description or "").strip()
    if not desc or desc == rid:
        return rid
    if desc.startswith(rid + " "):
        return desc
    return f"{rid} {desc}".strip()


def _generate_decoy(record: SeqRecord) -> SeqRecord:
    seq_str = str(record.seq)
    if seq_str.startswith("M"):
        decoy_seq = "M" + seq_str[1:][::-1]
    else:
        decoy_seq = seq_str[::-1]

    title = _record_title(record)
    if title.lower().startswith("rev_"):
        decoy_title = title
    else:
        decoy_title = f"rev_{title}"
    parts = decoy_title.split(None, 1)
    decoy_id = parts[0]
    decoy_desc = parts[1] if len(parts) > 1 else ""
    return SeqRecord(Seq(decoy_seq), id=decoy_id, name=decoy_id, description=decoy_desc)


def _read_records(fasta_path: str) -> List[SeqRecord]:
    if not fasta_path:
        return []
    path = Path(fasta_path)
    if not path.exists():
        log.warning(f"FASTA file not found: {fasta_path}")
        return []
    if path.is_dir():
        log.warning(f"FASTA path is a directory, skipping: {fasta_path}")
        return []
    return list(SeqIO.parse(str(path), "fasta"))


def _append_target_decoy(
    merged: List[SeqRecord],
    records: List[SeqRecord],
    pe_tag: str,
) -> None:
    for rec in records:
        rec.description = tag_header(rec.description, pe_tag)
        merged.append(rec)
        decoy = _generate_decoy(rec)
        decoy.description = tag_header(decoy.description, pe_tag)
        merged.append(decoy)


def assemble_database(
    circ_orf_path: str,
    linear_protein_path: str,
    output_dir: str,
    contam_path: Optional[str] = None,
) -> Path:
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "decoys-contam-database.fasta"

    merged: List[SeqRecord] = []

    circ_records = _read_records(circ_orf_path)
    _append_target_decoy(merged, circ_records, "PE=2")
    log.info(f"circRNA: {len(circ_records)} targets + {len(circ_records)} decoys")

    linear_records = _read_records(linear_protein_path)
    _append_target_decoy(merged, linear_records, "PE=1")
    log.info(f"linear proteins: {len(linear_records)} targets + {len(linear_records)} decoys")

    contam_fasta = Path(contam_path) if contam_path else _DEFAULT_CONTAM_FASTA
    if not contam_fasta.exists():
        raise RuntimeError(f"Contaminant FASTA not found: {contam_fasta}")
    contam_records = _read_records(str(contam_fasta))
    if not contam_records:
        raise RuntimeError(f"Contaminant FASTA is empty or unreadable: {contam_fasta}")
    _append_target_decoy(merged, contam_records, "PE=1")
    log.info(f"contaminants: {len(contam_records)} targets + {len(contam_records)} decoys ({contam_fasta})")

    if not merged:
        raise RuntimeError("No sequences to write. At least one input FASTA must be provided.")

    SeqIO.write(merged, str(out_path), "fasta")
    log.info(f"Unified database written → {out_path}  ({len(merged)} total sequences)")
    return out_path.resolve()
