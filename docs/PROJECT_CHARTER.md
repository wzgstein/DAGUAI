# Project Charter

**Status date:** 2026-08-28  
**Project state:** confirmation-stage audit; architecture training frozen.

## 1. Object of study

The project studies whether single-cell foundation models contain **target-specific relation information that is usable for predicting genetic-perturbation responses**, after separating five different objects:

1. **Availability (`A`)** — a relation source such as CORUM identifies partners whose perturbation responses transfer better than matched random partners.
2. **Static accessibility (`S`)** — a model's static gene-token geometry retrieves useful neighbors above dimension- and covariance-matched controls.
3. **Value contribution (`V`)** — baseline expression/value encoding contributes beyond the static token.
4. **Contextual increment (`C`)** — transformer/decoder computation contributes beyond `S+V`.
5. **Semantic-residual transfer (`I`, `T`)** — `C` remains useful after external semantic priors and survives independent data, support, split, and contamination checks.

A representation is not promoted merely because it has a positive raw score. Promotion requires a positive **nested increment** on the same targets, folds, nuisance payment, and statistical units.

## 2. Governing contradiction

The current evidence shows positive relation structure and positive static scGPT accessibility on a narrow support, while every fixed Round 23 contextual readout is no better than the same model's static table under complex-cluster reanalysis.

The project therefore faces a specific contradiction:

> Relation information is present and partly encoded in token geometry, but the contextual computation intended to exploit cell state has not demonstrated incremental use of that information.

The next experiment must decide whether this is a real property of the model or an artifact of the previous state construction, support intersection, outcome-based filtering, and non-fold-pure preprocessing.

## 3. Decisive local joint

The decisive comparison is:

```text
Q + S + V     versus     Q + S + V + C
```

under all of the following conditions:

- actual control-cell state rather than a median of perturbations;
- all-finite targets as the primary population;
- fold-pure response preprocessing;
- common target support as the primary estimand;
- relation-component-aware generalization and inference;
- representations and hyperparameters locked before confirmation;
- independent HepG2 and Jurkat confirmation data.

## 4. Workstreams

### Track A — relation availability and within-relation retrieval

Purpose: establish that a relation source identifies response-transfer partners and that static geometry can access them.

- Partner responses may appear in the training pool by design.
- The estimand is **within-relation interpolation**, not unseen-relation generalization.
- Inference is clustered by full relation-graph connected component.
- The CORUM partner mean is called a **partner reference**, never a ceiling.

### Track B — learned generalization and nested attribution

Purpose: test whether representation features predict responses for relation components excluded from training.

- Outer split is relation-component-held-out.
- Explicit response averaging over test-component partners is prohibited.
- Models are nested: `S`, `S+V`, `S+V+C`, `Q`, `Q+S+V`, `Q+S+V+C`.
- The primary promotion estimand is the paired held-out improvement of the final nested model.

### Track C — support and contamination audit

Purpose: quantify how much of a result is created by target selection, vocabulary coverage, knowledge-source overlap, or possible pretraining exposure.

- common support is primary;
- model-native support is secondary and reported with coverage;
- every result is bound to a target-universe hash;
- checkpoint and dataset dates are registered;
- non-CORUM relation sensitivity is required before a broad biological claim.

## 5. Frozen routes

The following routes remain frozen until the Round 24 promotion gate passes:

- new gene tokenizer;
- RMT-based training objective;
- exposure or contextual-utilization objective retraining;
- large adapter or end-to-end fine-tuning;
- architecture search over layers/readouts on confirmation data.

A small cross-fitted linear readout is an evaluation instrument, not a reopened training route.

## 6. Claim discipline

Every claim must carry:

- data version and hash;
- checkpoint and representation identity;
- support-universe hash and coverage;
- target filter;
- split unit and grouping graph;
- nuisance-payment protocol;
- primary statistical unit;
- raw paired estimate and uncertainty;
- whether the analysis was development, sensitivity, or confirmation.

No claim may be upgraded from exploratory to confirmatory after inspecting confirmation labels.

## 7. Stop condition

If the locked contextual increment fails on the independent confirmation pair, the project closes as a negative utilization audit. Static accessibility and relation availability may remain positive claims, but tokenizer/RMT development stays closed.
