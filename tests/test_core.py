import pandas as pd
import pytest

from relation_audit.contracts import ContractError, SplitContract
from relation_audit.inference import paired_cluster_inference
from relation_audit.splits import build_bipartite_components, grouped_kfold_assignments


def test_overlapping_relations_form_one_component():
    components = build_bipartite_components(
        [("A", "c1"), ("B", "c1"), ("B", "c2"), ("C", "c2"), ("D", "c3")]
    )
    assert components["A"] == components["B"] == components["C"]
    assert components["D"] != components["A"]


def test_relation_component_cannot_cross_folds():
    with pytest.raises(ContractError, match="cross outer folds"):
        SplitContract(
            "bad",
            "relation_component",
            {"A": "g1", "B": "g1"},
            {"A": 0, "B": 1},
            "forbidden",
            True,
        )


def test_grouped_assignments_preserve_components():
    memberships = []
    for index in range(6):
        memberships.extend(
            [(f"G{index}a", f"C{index}"), (f"G{index}b", f"C{index}")]
        )
    components = build_bipartite_components(memberships)
    folds = grouped_kfold_assignments(sorted(components), components, n_splits=3)
    for group in set(components.values()):
        assert len({folds[target] for target, value in components.items() if value == group}) == 1


def test_cluster_inference_is_deterministic_and_positive():
    frame = pd.DataFrame(
        {
            "cluster": ["a", "a", "b", "b", "c", "c", "d", "d"],
            "difference": [0.2, 0.1, 0.3, 0.4, 0.15, 0.2, 0.25, 0.3],
        }
    )
    first = paired_cluster_inference(
        frame,
        "difference",
        "cluster",
        bootstrap_repetitions=2000,
        signflip_repetitions=2000,
        seed=7,
    )
    second = paired_cluster_inference(
        frame,
        "difference",
        "cluster",
        bootstrap_repetitions=2000,
        signflip_repetitions=2000,
        seed=7,
    )
    assert first == second
    assert first.mean_difference > 0
    assert first.ci_low > 0
