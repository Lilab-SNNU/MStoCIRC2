

_MEANINGLESS_GENE_VALUES = frozenset({"", "none", "unknown", "na", "nan", "null"})


def is_meaningful_gene(gene: str) -> bool:
    return (gene or "").strip().lower() not in _MEANINGLESS_GENE_VALUES


def strip_gene_from_orf_header(header: str) -> str:
    match = HEADER_RE.match(header.split()[0])
    if not match:
        return header
    tag = "ORF" + match.group("orfnum")
    suffix = match.group("suffix") or ""
    if suffix:
        tag += "_" + suffix
    rebuilt = f"{match.group('base')}-{tag}({match.group('start')},{match.group('end')})"
    bps = match.group("bps") or ""
    if bps:
        rebuilt += ":" + bps
    return rebuilt


def extract_gene_from_orf_header(header: str) -> str:
    parts = header.split(":")
    if len(parts) >= 3:
        gene = parts[1].strip()
        if is_meaningful_gene(gene):
            return gene
    return "unknown"


def records_have_meaningful_gene(records: List[Dict[str, Any]]) -> bool:
    for record in records:
        for key in ("header", "source_header"):
            header = record.get(key, "")
            if header and is_meaningful_gene(extract_gene_from_orf_header(header)):
                return True
    return False


def normalize_orf_record_genes(records: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    if records_have_meaningful_gene(records):
        return records
    normalized: List[Dict[str, Any]] = []
    for record in records:
        updated = dict(record)
        for key in ("header", "source_header"):
            if key in updated and updated[key]:
                updated[key] = strip_gene_from_orf_header(updated[key])
        normalized.append(updated)
    return normalized
