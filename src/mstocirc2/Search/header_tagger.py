"""Regex-based PE tag processing: strip existing PE=\\d+, normalise spaces, append new PE tag."""
from __future__ import annotations
import re

_PE_PATTERN = re.compile(r"\bPE=\d+")
_SPACES = re.compile(r" +")


def tag_header(header: str, pe_tag: str) -> str:
    header = _PE_PATTERN.sub("", header)
    header = _SPACES.sub(" ", header).strip()
    if not header:
        return pe_tag
    return f"{header} {pe_tag}"
