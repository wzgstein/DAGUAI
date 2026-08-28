# Relation Availability and Utilization Audit for Single-Cell Foundation Models

> 中文摘要：本分支把原项目重构为一个**支持感知、干预落地、可归因**的审计框架。核心问题不是“哪个模型分数最高”，而是：关系信号是否存在、静态 token 是否可访问、contextual 计算是否带来增量、语义先验之后是否仍有剩余增量，以及这些结论能否在独立数据与污染校正后成立。

## Current decision

The project continues as an audit, not as a new-tokenizer project.

- **Interventional relation availability:** supported on the current development data.
- **Narrow-support static accessibility:** supported for scGPT on the Round 23 intersection.
- **Fixed contextual-over-static increment:** not supported at the current gate.
- **Tokenizer/RMT training:** frozen until the preregistered confirmation gate passes.

The governing comparison is nested rather than substitutive:

```text
S       = static gene-token geometry
S + V   = static geometry plus baseline-expression/value encoding
S + V+C = static/value terms plus contextual transformation
Q       = external semantic priors (for example GenePT / GO-C)
```

The decisive question is whether `C` improves held-out perturbation-response prediction after paying for `S`, `V`, and `Q` on the same support and split.

## Repository map

```text
docs/PROJECT_CHARTER.md                 strategic scope and claim boundaries
docs/ROUND23_INDEPENDENT_REANALYSIS.md  post-hoc complex-cluster audit
docs/ROUND24_PREREGISTRATION.md         locked confirmation design
config/round24_confirmation.yaml        machine-readable primary choices
schema/experiment_manifest.schema.json  provenance/support/split contract
src/relation_audit/                     reusable contracts, splits, metrics, inference
scripts/reanalyze_round23.py            reproducible complex-level reanalysis
tests/                                  leakage and inference unit tests
results/round23_reanalysis/             frozen audit tables and verdict
```

## Reproduce the Round 23 audit

```bash
python -m pip install -e '.[dev]'
python scripts/reanalyze_round23.py \
  --mask /path/to/mask-artifact/ALL_GENE_LEVEL_RESULTS.csv \
  --zero /path/to/zero-artifact/ALL_GENE_LEVEL_RESULTS.csv \
  --out /tmp/round23_reanalysis \
  --bootstrap-reps 50000 \
  --signflip-reps 50000
pytest
```

The script treats the assigned CORUM complex as the outer resampling unit. This is an audit of the executed Round 23 object; it does not repair overlapping-complex membership retrospectively.

## Non-claims

This branch does **not** currently claim that:

- scGPT contextual states are universally uninformative;
- a perturbation-consensus pseudo-state is equivalent to an individual control cell;
- CORUM partner averaging is a mathematical ceiling;
- within-relation retrieval measures unseen-complex generalization;
- a per-target oracle is a deployable representation selector;
- a new tokenizer, RMT objective, or foundation model is warranted.

## Decision rule

Model training is reopened only if a representation and readout locked on development data show a positive, practically material contextual increment on both independent confirmation cell lines, survive relation-component inference and common-support analysis, and remain positive after semantic priors are included.

See [`docs/ROUND24_PREREGISTRATION.md`](docs/ROUND24_PREREGISTRATION.md) for the full gate.
