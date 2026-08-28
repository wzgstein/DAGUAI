from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Mapping


class ContractError(ValueError):
    """Raised when an audit object is internally inconsistent."""


def _nonempty_unique(values: Iterable[str], name: str) -> tuple[str, ...]:
    cleaned = tuple(str(v).strip() for v in values if str(v).strip())
    if not cleaned:
        raise ContractError(f"{name} must not be empty")
    if len(set(cleaned)) != len(cleaned):
        raise ContractError(f"{name} contains duplicates")
    return cleaned


@dataclass(frozen=True)
class SupportContract:
    universe: str
    targets: tuple[str, ...]
    filter_name: str
    required_blocks: tuple[str, ...]
    coverage_denominator: int | None = None

    def __post_init__(self) -> None:
        if self.universe not in {"common", "model_native", "union"}:
            raise ContractError(f"unsupported support universe: {self.universe}")
        object.__setattr__(self, "targets", _nonempty_unique(self.targets, "targets"))
        object.__setattr__(
            self,
            "required_blocks",
            _nonempty_unique(self.required_blocks, "required_blocks"),
        )
        if not self.filter_name.strip():
            raise ContractError("filter_name must not be empty")
        if self.coverage_denominator is not None:
            if self.coverage_denominator < len(self.targets):
                raise ContractError("coverage_denominator cannot be smaller than target count")


@dataclass(frozen=True)
class SplitContract:
    name: str
    unit: str
    groups: Mapping[str, str]
    folds: Mapping[str, int]
    partner_access: str
    fold_pure_preprocessing: bool

    def __post_init__(self) -> None:
        if self.unit not in {"target", "relation_component", "cell_line", "study", "donor"}:
            raise ContractError(f"unsupported split unit: {self.unit}")
        if self.partner_access not in {"allowed", "forbidden", "not_applicable"}:
            raise ContractError(f"unsupported partner_access: {self.partner_access}")
        if set(self.groups) != set(self.folds):
            raise ContractError("groups and folds must describe exactly the same targets")
        if not self.folds or min(self.folds.values()) < 0:
            raise ContractError("split must contain non-negative fold ids")
        if self.unit == "relation_component":
            group_to_folds: dict[str, set[int]] = {}
            for target, group in self.groups.items():
                group_to_folds.setdefault(group, set()).add(int(self.folds[target]))
            broken = [group for group, ids in group_to_folds.items() if len(ids) != 1]
            if broken:
                raise ContractError(
                    "relation components cross outer folds: "
                    + ", ".join(sorted(broken)[:5])
                )
        if not self.fold_pure_preprocessing:
            raise ContractError("confirmatory splits require fold-pure preprocessing")


@dataclass(frozen=True)
class RepresentationContract:
    name: str
    role: str
    dimension: int
    locked: bool = True

    def __post_init__(self) -> None:
        if self.role not in {"semantic", "static", "value", "context", "null"}:
            raise ContractError(f"unsupported representation role: {self.role}")
        if self.dimension <= 0:
            raise ContractError("representation dimension must be positive")


@dataclass(frozen=True)
class InferenceContract:
    cluster: str
    bootstrap_repetitions: int = 50_000
    signflip_repetitions: int = 50_000
    seed: int = 20_260_828

    def __post_init__(self) -> None:
        if not self.cluster.strip():
            raise ContractError("cluster must not be empty")
        if self.bootstrap_repetitions < 1_000:
            raise ContractError("bootstrap_repetitions must be at least 1000")
        if self.signflip_repetitions < 1_000:
            raise ContractError("signflip_repetitions must be at least 1000")


@dataclass(frozen=True)
class AuditContract:
    experiment_id: str
    stage: str
    truth_boundary: str
    support: SupportContract
    split: SplitContract
    representations: tuple[RepresentationContract, ...]
    inference: InferenceContract
    primary_contrast: tuple[str, str]
    metadata: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.stage not in {"development", "sensitivity", "confirmation"}:
            raise ContractError(f"unsupported stage: {self.stage}")
        if len(self.truth_boundary.strip()) < 20:
            raise ContractError("truth_boundary is too short to be informative")
        names = _nonempty_unique((r.name for r in self.representations), "representations")
        before, after = self.primary_contrast
        if before not in names or after not in names:
            raise ContractError("primary contrast references an undeclared representation")
        if self.stage == "confirmation" and not all(r.locked for r in self.representations):
            raise ContractError("all confirmation representations must be locked")
        if set(self.support.targets) != set(self.split.folds):
            raise ContractError("support targets and split targets differ")
