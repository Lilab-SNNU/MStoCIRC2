"""Workflow strategy dataclass and registry mapping --strategy values to the four pipeline stages."""
from __future__ import annotations
from dataclasses import dataclass, asdict
from typing import Dict


@dataclass(frozen=True)
class WorkflowStrategy:
    expression_matrix_type: str
    imputation: str
    dea_method: str

    def as_dict(self) -> Dict[str, str]:
        return asdict(self)


STRATEGY_REGISTRY: Dict[str, WorkflowStrategy] = {
    "fragpipe-dda":    WorkflowStrategy("directLFQ",                "SeqKNN",  "DEqMS"),
    "maxquant-dda":    WorkflowStrategy("directLFQ",                "Impseq",  "DEqMS"),
    "diann-dia":       WorkflowStrategy("directLFQ",                "MinDet",  "limma"),
    "spectronaut-dia": WorkflowStrategy("directLFQ",                "Impseq",  "ROTS"),
    "fragpipe-tmt":    WorkflowStrategy("TMT-Integrator abundance", "SeqKNN",  "limma"),
    "maxquant-tmt":    WorkflowStrategy("reporter intensity",       "bpca",    "proDA"),
    "generic":         WorkflowStrategy("directLFQ",                "MinProb", "limma"),
}


def resolve_strategy(name: str) -> WorkflowStrategy:
    if name not in STRATEGY_REGISTRY:
        valid = ", ".join(sorted(STRATEGY_REGISTRY.keys()))
        raise ValueError(f"Unknown --strategy '{name}'. Valid values: {valid}")
    return STRATEGY_REGISTRY[name]
