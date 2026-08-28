from __future__ import annotations

from collections import defaultdict, deque
from collections.abc import Iterable, Mapping, Sequence

import numpy as np
from sklearn.model_selection import GroupKFold


class SplitError(ValueError):
    """Raised when relation groups cannot support the requested split."""


def build_bipartite_components(
    memberships: Iterable[tuple[str, str]],
) -> dict[str, str]:
    """Return a gene-to-component mapping while retaining overlapping relations."""

    adjacency: dict[str, set[str]] = defaultdict(set)
    genes: set[str] = set()
    for gene_raw, relation_raw in memberships:
        gene = str(gene_raw).strip()
        relation = str(relation_raw).strip()
        if not gene or not relation:
            continue
        gene_node = f"g:{gene}"
        relation_node = f"r:{relation}"
        adjacency[gene_node].add(relation_node)
        adjacency[relation_node].add(gene_node)
        genes.add(gene)
    if not genes:
        raise SplitError("no valid gene-relation memberships")

    import hashlib

    seen: set[str] = set()
    result: dict[str, str] = {}
    for start in sorted(adjacency):
        if start in seen:
            continue
        queue: deque[str] = deque([start])
        seen.add(start)
        nodes: list[str] = []
        while queue:
            node = queue.popleft()
            nodes.append(node)
            for next_node in sorted(adjacency[node]):
                if next_node not in seen:
                    seen.add(next_node)
                    queue.append(next_node)
        digest = hashlib.sha256("\n".join(sorted(nodes)).encode()).hexdigest()[:16]
        component_id = f"component-{digest}"
        for node in nodes:
            if node.startswith("g:"):
                result[node[2:]] = component_id
    return result


def grouped_kfold_assignments(
    targets: Sequence[str],
    groups: Mapping[str, str],
    n_splits: int,
) -> dict[str, int]:
    """Assign deterministic relation-component-held-out folds."""

    ordered = sorted(dict.fromkeys(str(target) for target in targets))
    if not ordered:
        raise SplitError("targets must not be empty")
    missing = [target for target in ordered if target not in groups]
    if missing:
        raise SplitError(f"missing groups for targets: {missing[:5]}")
    unique_groups = {groups[target] for target in ordered}
    if len(unique_groups) < n_splits:
        raise SplitError(f"need at least {n_splits} groups, found {len(unique_groups)}")

    x = np.zeros((len(ordered), 1), dtype=float)
    y = np.zeros(len(ordered), dtype=float)
    group_array = np.asarray([groups[target] for target in ordered], dtype=object)
    fold_map: dict[str, int] = {}
    splitter = GroupKFold(n_splits=n_splits)
    for fold, (_, test_indices) in enumerate(splitter.split(x, y, groups=group_array)):
        for index in test_indices:
            fold_map[ordered[int(index)]] = fold

    for group in unique_groups:
        fold_ids = {fold_map[target] for target in ordered if groups[target] == group}
        if len(fold_ids) != 1:
            raise AssertionError(f"group {group} crosses folds")
    return fold_map


def assert_no_group_leakage(
    train_targets: Sequence[str],
    test_targets: Sequence[str],
    groups: Mapping[str, str],
) -> None:
    train_groups = {groups[target] for target in train_targets}
    test_groups = {groups[target] for target in test_targets}
    overlap = train_groups.intersection(test_groups)
    if overlap:
        raise SplitError(f"group leakage detected: {sorted(overlap)[:5]}")
