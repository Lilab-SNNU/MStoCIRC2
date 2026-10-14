from ..orf_prediction.bsj_truncator import extract_gene_from_orf_header, is_meaningful_gene
from .evidence_scorer import (
    build_source_bsj_location,
    is_peptide_crossing_bsj,
    parse_bsj_location,
)

def build_identity_mapping_from_orf(
    circ_orf: str,
) -> tuple[dict[str, str], dict[str, str], dict[str, tuple[int, int]]]:
    dic_output_to_source: dict[str, str] = {}
    dic_source_aa: dict[str, str] = {}
    dic_source_orf_coords: dict[str, tuple[int, int]] = {}
    if not Path(circ_orf).exists():
        return dic_output_to_source, dic_source_aa, dic_source_orf_coords
    cur_id = ""
    cur_seq = ""

    def _commit(header: str, seq: str) -> None:
        if not header or not seq:
            return
        dic_output_to_source[header] = header
        dic_source_aa[header] = seq
        match = re.search(r"\((-?\d+),(-?\d+)\)", header)
        if match:
            coords = (int(match.group(1)), int(match.group(2)))
            dic_source_orf_coords[header] = coords

    with Path(circ_orf).open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line.startswith(">"):
                _commit(cur_id, cur_seq)
                cur_id = line[1:].strip()
                cur_seq = ""
            elif cur_id:
                cur_seq += line
        _commit(cur_id, cur_seq)
    return dic_output_to_source, dic_source_aa, dic_source_orf_coords


def resolve_full_orf_sequence(source_header: str, dic_source_aa: dict[str, str]) -> str:
    return dic_source_aa.get(source_header, "")


def resolve_bsj_truncated_orf_sequence(
    source_header: str,
    dic_orf_seqs_fallback: dict[str, str],
    dic_source_search_orfs: dict[str, list[str]],
) -> str:
    seq = dic_orf_seqs_fallback.get(source_header, "")
    if not seq:
        seq = dic_orf_seqs_fallback.get(source_header.split("(")[0], "")
    if not seq:
        for out_id in dic_source_search_orfs.get(source_header, []):
            seq = dic_orf_seqs_fallback.get(out_id, "") or dic_orf_seqs_fallback.get(out_id.split("(")[0], "")
            if seq:
                break
    return seq


def compute_peptide_coverage_on_corf(corf: str, peptide_entries_str: str) -> float:
    corf_clean = corf.strip("*").upper()
    corf_clean_il = corf_clean.replace("I", "L")
    corf_s = len(corf_clean)
    if corf_s == 0:
        return 0.0
    coverage_mask = [False] * corf_s
    entries = [x.strip().split("###") for x in peptide_entries_str.strip(";").split(";") if x.strip()]
    for parts in entries:
        if len(parts) < 2:
            continue
        peptide = parts[1].upper()
        peptide_il = peptide.replace("I", "L")
        start = 0
        while True:
            idx = corf_clean_il.find(peptide_il, start)
            if idx == -1:
                break
            end = idx + len(peptide)
            for k in range(idx, min(end, corf_s)):
                coverage_mask[k] = True
            start = idx + 1
    return sum(coverage_mask) / corf_s

        return "none", "none", "none", "none", "none", "none", {}, 0, "", []
    peptide_rows: list[tuple[str, str, int, int, str, str]] = []
    for peptide in sorted(all_unique_peps):
        total_count = sum(dic_pro_sore[peptide].values())
        ms_hit_count = len(dic_pro_ms_counts.get(peptide, set()))
        ms_evidence = ";".join(
            f"{mgf_id}:{dic_pro_sore[peptide][mgf_id]}"
            for mgf_id in sorted(dic_pro_sore[peptide].keys())
        )
        if ms_evidence:
            ms_evidence += ";"
        bsj_spanning = "yes" if peptide in stats_cross["unique_pep"] else "no"
        peptide_rows.append(
            (peptide.upper(), circ_orf_id, ms_hit_count, total_count, ms_evidence, bsj_spanning)
        )
        peptide_rows,
def write_circ_peptide_tsv(peptide_rows: list[tuple[str, str, int, int, str, str]], file_path_out: str) -> None:
    header = "\t".join(
        [
            "circ_peptide",
            "circ_ORF_ID",
            "MS_count",
            "spectral_count",
            "MS_sample_evidence",
            "BSJ_Spanning",
        ]
    )
    with Path(file_path_out).open("w", encoding="utf-8") as handle:
        handle.write(header + "\n")
        for peptide, orf_id, ms_count, spectral_count, ms_evidence, bsj_spanning in peptide_rows:
            handle.write(
                f"{peptide}\t{orf_id}\t{ms_count}\t{spectral_count}\t{ms_evidence}\t{bsj_spanning}\n"
            )


def calculate_and_write_scores(
    results_buffer: list[dict[str, Any]],
    file_path_out: str,
    has_info: bool,
    include_gene: bool = True,
) -> None:
    header_cols = ["circ_ORF_ID", "search_ORF_IDs", "circRNA_ID"]
    if include_gene:
        header_cols.append("gene")
            "BSJ_peptide_coverage",
            "cORF_BSJ_position",
            "BSJ_peptide_support",
        score_m6a = 0.15 * (1 - math.exp(-0.5 * stats["m"]))
        score_cov = 0.10 * stats["cov"]
            parts.pop(4 if include_gene else 3)
    include_gene: bool = True,
    cross_unique_pep = int((data["detailed_stats"].get("cross") or {}).get("unique_pep", 0))
    bsj_peptide_support = "True" if cross_unique_pep > 0 else "False"
    row_parts = [
        circ_orf_id,
        data["search_orf_ids"],
        data["circ_info_id"],
    ]
    if include_gene:
        row_parts.append(data["gene"])
    row_parts.extend(
        [
            data["circ_strand"],
            ires_str,
            circ_m6a_count,
            circ_orf_seq,
            str(orf_len),
            data["marge_peptide_bsj_multiple"],
            f"{data['bsj_peptide_coverage']:.4f}",
            data["bsj_location"],
            str(data["raw_sum"]),
            str(data["number_all"]),
            str(data["unique_peptide_count"]),
            data["peptide_seq_str"],
            bsj_peptide_support,
            data["mgf_score_all"],
        ]
    row_circ_predict = "\t".join(row_parts)
        "cov": float(data["bsj_peptide_coverage"]),
        idx_seq_cov = header_parts.index("BSJ_peptide_coverage")
        idx_bsj_pos = header_parts.index("cORF_BSJ_position")
        return
    try:
        idx_circ_gene = header_parts.index("gene")
    except ValueError:
        idx_circ_gene = None
        required_indices = [idx_corf_id, idx_corf_seq, idx_seq_cov, idx_bsj_pos]
        if idx_circ_gene is not None:
            required_indices.append(idx_circ_gene)
        if len(parts) <= max(required_indices):
        circ_gene = parts[idx_circ_gene] if idx_circ_gene is not None else "NA"
        if idx_circ_gene is not None:
            ax.text(1, -y_max, f"Gene_name:{circ_gene}", fontsize=12)
        ax.text(1, -y_max / (0.66 if idx_circ_gene is not None else 1), f"Sequence cover:{float_circ_seq_coverage:.2f}%", fontsize=12)
        return
    try:
        idx_circ_gene = header_parts.index("gene")
    except ValueError:
        idx_circ_gene = None
        required_indices = [idx_full_orf, idx_corf_seq, idx_corf_map]
        if idx_circ_gene is not None:
            required_indices.append(idx_circ_gene)
        if len(parts) <= max(required_indices):
        circ_gene = parts[idx_circ_gene] if idx_circ_gene is not None else ""
    if circ_mapping and Path(circ_mapping).exists():
        dic_output_to_source, dic_source_aa, dic_source_orf_coords = load_mapping_context(circ_mapping)
    else:
        dic_output_to_source, dic_source_aa, dic_source_orf_coords = build_identity_mapping_from_orf(circ_orf)
    dic_circ_orf_bsj_location = build_source_bsj_location(circ_orf, circ_mapping or "", circ_seq)
    all_peptide_rows: list[tuple[str, str, int, int, str, str]] = []
        full_seq = resolve_full_orf_sequence(source_header, dic_source_aa)
        if not full_seq:
        truncated_seq = resolve_bsj_truncated_orf_sequence(
            source_header,
            dic_orf_seqs_fallback,
            dic_source_search_orfs,
        )
            peptide_rows,
        ) = circ_orf_tsv_peptide_map(source_header, full_seq, dic_id_peptide, dic_circ_orf_bsj_location)
        if marge_peptide_bsj_multiple == "none":
        bsj_peptide_coverage = 0.0
        if truncated_seq and source_header in dic_id_peptide:
            bsj_peptide_coverage = compute_peptide_coverage_on_corf(
                truncated_seq,
                dic_id_peptide[source_header],
            )
        all_peptide_rows.extend(peptide_rows)

        gene = extract_gene_from_orf_header(source_header)
            "bsj_peptide_coverage": bsj_peptide_coverage,
        candidate_orfs.append((source_header, full_seq))
            "mass-spec parsing and canonical-protein filtering, so `circ_translation_evidence.tsv` was not generated."
    include_gene = any(is_meaningful_gene(data["gene"]) for data in precomputed.values())
        result = process_orf(
            orf_id,
            orf_seq,
            dic_ires_results,
            precomputed,
            dic_m6a_results,
            run_ires,
            include_gene=include_gene,
        )
    circ_translation_evidence_path = output_dir / "circ_translation_evidence.tsv"
    circ_peptide_path = output_dir / "circ_peptide.tsv"
    calculate_and_write_scores(results_buffer, str(circ_translation_evidence_path), has_info, include_gene)
    write_circ_peptide_tsv(all_peptide_rows, str(circ_peptide_path))
    if not circ_translation_evidence_path.exists():
            f"Translation evaluation did not generate the expected result file: '{circ_translation_evidence_path}'."
        )
    if not circ_peptide_path.exists():
        raise MStoCIRC2Error(
            f"Translation evaluation did not generate the expected result file: '{circ_peptide_path}'."
    draw_msgf_second(str(circ_translation_evidence_path), 10, str(output_dir), dic_id_peptide)
    draw_circ_corf(str(output_dir), str(circ_translation_evidence_path), dic_output_to_source, dic_circ_id_seq)
    log.info("Evaluation stage completed: %s", circ_translation_evidence_path)