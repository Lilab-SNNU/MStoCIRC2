import re
from pathlib import Path
from typing import Dict, List

from ..orf_prediction.bsj_truncator import (
    HEADER_RE,
    orf_coords_to_zero_based,
    orf_effective_coding_nt_len,
)
from .mass_spec_parser import extract_circ_info_id


def parse_bsj_location(bsj_location_str: str) -> List[float]:
    if not bsj_location_str or str(bsj_location_str).strip() == '0':
        return []
    res = []
    for x in str(bsj_location_str).split('|'):
        x = x.strip().split()[0]
        if x:
            try:
                res.append(float(x))
            except ValueError:
                pass
    return res


def is_peptide_crossing_bsj(pep_start_0based: int, pep_end_0based: int, bsj_list: List[float]) -> bool:
    if not bsj_list:
        return False
    for b_val in bsj_list:
        b_int = int(b_val)
        b_dec = round(b_val - b_int, 1)
        b_0based = b_int - 1
        if b_dec == 0.3:
            if pep_start_0based <= b_0based and pep_end_0based >= b_0based + 2:
                return True
        else:
            if pep_start_0based <= b_0based < pep_end_0based:
                return True
    return False


def is_standard_bsj_string(bsj: str | None) -> bool:
    if not bsj or str(bsj).strip() == '0':
        return False
    for token in str(bsj).split('|'):
        token = token.strip()
        if token and re.match(r'^[\d.|]+$', token) and '.' in token:
            return True
    return False


def extract_bsj_token_from_header(header: str) -> str | None:
    header_token = header.split()[0]
    match = HEADER_RE.match(header_token)
    if match:
        bps = (match.group('bps') or '').strip()
        if bps and re.match(r'^[\d.|]+$', bps) and '.' in bps:
            return bps
    parts = header.split(':')
    if len(parts) < 2:
        return None
    candidates: List[str] = []
    if len(parts) >= 3:
        candidates.extend([parts[2], parts[-1]])
    else:
        candidates.append(parts[1])
    for candidate in candidates:
        if candidate and re.match(r'^[\d.|]+$', candidate) and '.' in candidate:
            return candidate
    return None


def compute_bsj_from_coords(header_str: str, dic_circ_len: Dict[str, int]) -> str | None:
    match = re.search(r'\((-?\d+),(-?\d+)\)', header_str)
    if not match:
        return None
    orf_start_1based = int(match.group(1))
    orf_end_1based = int(match.group(2))
    orf_start_0, orf_nt_len = orf_coords_to_zero_based(orf_start_1based, orf_end_1based)
    if orf_nt_len <= 0:
        return '0'
    base_id = header_str.split('(')[0]
    circ_id_local = base_id.split('-ORF')[0] if '-ORF' in base_id else base_id
    circ_id_stripped = extract_circ_info_id(circ_id_local)
    circ_len = dic_circ_len.get(circ_id_local, 0) or dic_circ_len.get(circ_id_stripped, 0)
    if circ_len <= 0:
        return None
    if orf_start_0 + orf_nt_len <= circ_len:
        return '0'
    effective_nt_len = orf_effective_coding_nt_len(header_str, orf_start_1based, orf_end_1based)
    breakpoints: List[str] = []
    j1_dist = circ_len - orf_start_0
    current_j_dist = j1_dist
    while current_j_dist < effective_nt_len:
        aa_idx = (current_j_dist // 3) + 1
        rem = current_j_dist % 3
        if aa_idx <= 0:
            aa_idx = 1
        bp_str = f"{aa_idx - 1}.3" if rem == 0 else f"{aa_idx}.{rem}"
        breakpoints.append(bp_str)
        current_j_dist += circ_len
    return '|'.join(breakpoints) if breakpoints else '0'


def _load_circ_len_map(circ_seq: str | None) -> Dict[str, int]:
    if not circ_seq:
        return {}
    circ_seq_path = Path(circ_seq)
    if not circ_seq_path.exists():
        return {}
    lengths: Dict[str, int] = {}
    current_name: str | None = None
    with circ_seq_path.open('r', encoding='utf-8', errors='ignore') as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith('>'):
                current_name = line[1:].split(':')[0]
                lengths[current_name] = 0
            elif current_name is not None:
                lengths[current_name] = lengths.get(current_name, 0) + len(line)
    return lengths


def _resolve_bsj_for_header(
    header: str,
    dic_circ_len: Dict[str, int],
    force_compute_all: bool,
) -> str:
    if force_compute_all:
        computed = compute_bsj_from_coords(header, dic_circ_len)
        return computed if computed is not None else '0'
    token = extract_bsj_token_from_header(header)
    if token:
        return token
    computed = compute_bsj_from_coords(header, dic_circ_len)
    if computed is not None:
        return computed
    return '0'


def build_source_bsj_location(
    circ_orf: str,
    circ_mapping: str,
    circ_seq: str | None = None,
) -> Dict[str, str]:
    dic_circ_len = _load_circ_len_map(circ_seq)
    header_samples: List[str] = []
    circ_mapping_path = Path(circ_mapping) if circ_mapping else None
    if circ_mapping_path and circ_mapping_path.exists():
        with circ_mapping_path.open('r', encoding='utf-8', errors='ignore') as handle:
            handle.readline()
            for line in handle:
                cols = line.strip('\n').split('\t')
                if len(cols) >= 2:
                    header_samples.append(cols[1].split()[0])
    else:
        circ_orf_path = Path(circ_orf)
        if circ_orf_path.exists():
            with circ_orf_path.open('r', encoding='utf-8', errors='ignore') as handle:
                for line in handle:
                    if line.startswith('>'):
                        header_samples.append(line.strip()[1:].split()[0])

    standard_count = sum(
        1 for header in header_samples if is_standard_bsj_string(extract_bsj_token_from_header(header))
    )
    force_compute_all = bool(header_samples) and standard_count == 0

    dic_bsj: Dict[str, str] = {}
    if circ_mapping_path and circ_mapping_path.exists():
        with circ_mapping_path.open('r', encoding='utf-8', errors='ignore') as handle:
            handle.readline()
            for line in handle:
                cols = line.strip('\n').split('\t')
                if len(cols) < 2:
                    continue
                source_header = cols[1].split()[0]
                source_base = source_header.split('(')[0]
                src_bsj = _resolve_bsj_for_header(source_header, dic_circ_len, force_compute_all)
                if src_bsj is not None:
                    dic_bsj[source_base] = src_bsj
                elif source_base not in dic_bsj:
                    out_base = cols[0].split('(')[0]
                    dic_bsj[source_base] = dic_bsj.get(out_base, '0')
                out_base = cols[0].split('(')[0]
                if out_base not in dic_bsj and source_base in dic_bsj:
                    dic_bsj[out_base] = dic_bsj[source_base]
    else:
        circ_orf_path = Path(circ_orf)
        if circ_orf_path.exists():
            with circ_orf_path.open('r', encoding='utf-8', errors='ignore') as handle:
                for line in handle:
                    if line.startswith('>'):
                        header = line.strip()[1:].split()[0]
                        base_id = header.split('(')[0]
                        dic_bsj[base_id] = _resolve_bsj_for_header(header, dic_circ_len, force_compute_all)
    return dic_bsj
