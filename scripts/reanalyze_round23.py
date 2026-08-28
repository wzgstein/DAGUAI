#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from relation_audit.inference import paired_cluster_inference

SEED = 20260828
KEY_REPRESENTATIONS = [
    "GenePT_NCBI_UniProt",
    "GO_C",
    "SG_Static",
    "SG_H00",
    "SG_D12",
    "SG_H03",
    "SG_H06",
    "SG_H12",
    "SG_CLSM06",
    "SG_MEANM12",
    "SG_DECM12",
    "SG_TARGETM12",
    "SG_Random",
]


def load_raw(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    required = {
        "cell_line",
        "perturbation_gene",
        "complex_id",
        "representation",
        "method",
        "cosine_gain",
    }
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")
    raw = frame.loc[frame["method"].eq("raw_knn")].copy()
    raw["complex_id"] = raw["complex_id"].astype(str)
    duplicated = raw.duplicated(
        ["cell_line", "perturbation_gene", "representation"],
        keep=False,
    )
    if duplicated.any():
        raise ValueError(f"{path}: duplicate raw_knn target/representation rows")
    return raw


def pivot_raw(frame: pd.DataFrame) -> pd.DataFrame:
    metadata = (
        frame[["cell_line", "perturbation_gene", "complex_id"]]
        .drop_duplicates()
        .set_index(["cell_line", "perturbation_gene"])
    )
    values = frame.pivot(
        index=["cell_line", "perturbation_gene"],
        columns="representation",
        values="cosine_gain",
    )
    return metadata.join(values, how="inner").reset_index()


def infer_table(
    wide: pd.DataFrame,
    representations: list[str],
    *,
    reference: str | None,
    bootstrap_reps: int,
    signflip_reps: int,
    seed: int,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for offset, representation in enumerate(representations):
        if representation not in wide:
            continue
        work = wide[["complex_id", representation] + ([reference] if reference else [])].copy()
        work["audit_value"] = (
            work[representation] - work[reference]
            if reference
            else work[representation]
        )
        result = paired_cluster_inference(
            work,
            value_col="audit_value",
            cluster_col="complex_id",
            bootstrap_repetitions=bootstrap_reps,
            signflip_repetitions=signflip_reps,
            seed=seed + offset,
        )
        rows.append(
            {
                "representation": representation,
                "reference": reference or "zero",
                "n_rows": result.n_rows,
                "n_complexes": result.n_clusters,
                "equal_complex_mean_difference": result.mean_difference,
                "ci_2.5": result.ci_low,
                "ci_97.5": result.ci_high,
                "positive_complex_fraction": result.positive_cluster_fraction,
                "sign_flip_p_two_sided": result.signflip_p_two_sided,
            }
        )
    return pd.DataFrame(rows).sort_values(
        "equal_complex_mean_difference",
        ascending=False,
        ignore_index=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mask", type=Path, required=True)
    parser.add_argument("--zero", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--bootstrap-reps", type=int, default=50_000)
    parser.add_argument("--signflip-reps", type=int, default=50_000)
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    mask = pivot_raw(load_raw(args.mask))
    zero = pivot_raw(load_raw(args.zero))

    key = infer_table(
        mask,
        KEY_REPRESENTATIONS,
        reference=None,
        bootstrap_reps=args.bootstrap_reps,
        signflip_reps=args.signflip_reps,
        seed=args.seed,
    )
    mask_context_reps = [
        "SG_H00",
        "SG_H03",
        "SG_D03",
        "SG_H06",
        "SG_D06",
        "SG_H09",
        "SG_D09",
        "SG_H12",
        "SG_D12",
        "SG_CLSM03",
        "SG_CLSM06",
        "SG_CLSM09",
        "SG_CLSM12",
        "SG_MEANM12",
        "SG_DECM12",
        "SG_TARGETM12",
    ]
    zero_context_reps = [
        "SG_H00",
        "SG_H03",
        "SG_D03",
        "SG_H06",
        "SG_D06",
        "SG_H09",
        "SG_D09",
        "SG_H12",
        "SG_D12",
        "SG_CLSZ03",
        "SG_CLSZ06",
        "SG_CLSZ09",
        "SG_CLSZ12",
        "SG_MEANZ12",
        "SG_DECZ12",
        "SG_TARGETZ12",
    ]
    mask_context = infer_table(
        mask,
        mask_context_reps,
        reference="SG_Static",
        bootstrap_reps=args.bootstrap_reps,
        signflip_reps=args.signflip_reps,
        seed=args.seed + 10_000,
    )
    zero_context = infer_table(
        zero,
        zero_context_reps,
        reference="SG_Static",
        bootstrap_reps=args.bootstrap_reps,
        signflip_reps=args.signflip_reps,
        seed=args.seed + 20_000,
    )
    semantic = infer_table(
        mask,
        ["GenePT_NCBI_UniProt", "GO_C", "SG_GeneNorm", "SG_Value", "SG_Random"],
        reference="SG_Static",
        bootstrap_reps=args.bootstrap_reps,
        signflip_reps=args.signflip_reps,
        seed=args.seed + 30_000,
    )

    outputs = {
        "round23_key_gains_complex_bootstrap.csv": key,
        "round23_mask_contextual_minus_static.csv": mask_context,
        "round23_zero_contextual_minus_static.csv": zero_context,
        "round23_semantic_minus_static.csv": semantic,
    }
    for name, frame in outputs.items():
        frame.to_csv(args.out / name, index=False)

    verdict = {
        "stage": "post_hoc_audit",
        "n_target_line_rows": int(mask.shape[0]),
        "n_complexes": int(mask["complex_id"].nunique()),
        "static_accessibility_positive": bool(
            key.loc[key["representation"].eq("SG_Static"), "ci_2.5"].iloc[0] > 0
        ),
        "any_mask_context_ci_lower_above_zero": bool((mask_context["ci_2.5"] > 0).any()),
        "any_zero_context_ci_lower_above_zero": bool((zero_context["ci_2.5"] > 0).any()),
        "contextual_promotion": "fail",
        "truth_boundary": (
            "Round 23 perturbation-consensus pseudo-state and target-value query; "
            "post-hoc assigned-complex audit, not preregistered confirmation."
        ),
    }
    (args.out / "VERDICT.json").write_text(json.dumps(verdict, indent=2) + "\n")
    print(json.dumps(verdict, indent=2))


if __name__ == "__main__":
    main()
