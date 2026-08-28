# Round 24 Preregistration: Fold-Pure Control-State Confirmation

**Status:** design locked in this document before confirmation labels are inspected.  
**Date:** 2026-08-28  
**Development data:** Replogle 2022 K562 and RPE1.  
**Independent confirmation data:** Nadig 2025 HepG2 and Jurkat.  
**Primary model:** whole-human scGPT checkpoint used in Round 23.  
**Primary relation source:** CORUM, with full overlapping-membership graph retained.

## 1. Confirmatory question

Does a fixed scGPT contextual transformation improve genetic-perturbation response prediction beyond the same model's static gene token, baseline-value encoding, and external semantic priors?

The confirmatory nested comparison is:

```text
M0 = Q + S + V
M1 = Q + S + V + C
```

where:

- `Q` is the locked semantic stack (GenePT and GO-C, block-normalized);
- `S` is the locked static scGPT gene embedding;
- `V = H00 - S` is the target's baseline-expression/value contribution in a control-cell context;
- `C = D12 - H00` is the locked final decoder contextual transformation;
- all blocks are fitted only through a cross-fitted linear readout; representation extraction remains frozen.

`D12` is selected as the primary contextual state because it is a trained final generative readout, not because it won an oracle search in confirmation data. `H03` is a single prespecified secondary sensitivity representation; it cannot independently promote the project.

## 2. Data populations

### 2.1 Development correction

Re-run the protocol on Replogle K562 and RPE1 to verify that the implementation removes known Round 23 leakage:

- use only control cells to construct the cell-line state;
- remove the `energy_test_p_value < 0.001` primary filter;
- do not select response features using held-out target responses;
- use the full overlapping relation graph for group construction;
- use dimension/covariance-matched null representations.

Development data may be used to set regularization grids, block scaling, neighbor count, and numerical tolerances. Once set, they are written to the configuration hash and cannot be changed for confirmation.

### 2.2 Confirmation pair

HepG2 and Jurkat are analyzed independently. No confirmation response labels may be used to choose target filters, layers, readouts, block weights, relation source, nuisance rank, regularization grid, or success thresholds.

The two cell lines are separate replications. A pooled result cannot rescue a sign reversal in either line.

### 2.3 Primary target population

A target enters the primary common-support universe when all conditions hold:

1. a finite perturbation-response vector can be estimated;
2. the target is a single-gene perturbation with at least the minimum cell count;
3. it maps unambiguously to the locked human gene-symbol table;
4. it is in the scGPT vocabulary and static embedding table;
5. its token is present in the locked control-state sequence policy;
6. it is covered by both semantic blocks;
7. it belongs to the CORUM relation graph and has at least two eligible partners in the development-defined support policy.

The target list and SHA-256 hash are emitted before model scoring.

### 2.4 Secondary populations

- model-native support, reported together with numerator/denominator coverage;
- significant-only targets, using the original `p < 0.001` rule, as sensitivity only;
- high- and low-control-expression strata;
- relation-degree strata;
- a non-CORUM physical/pathway relation source selected and frozen after a label-free coverage audit.

No secondary population can override a failed primary gate.

## 3. Response construction and nuisance payment

### 3.1 Control reference

For each cell line, use only annotated control cells. The primary control profile is the gene-wise median of normalized control expression. Control-cell bootstrap states are used as a sensitivity analysis of representation stability.

### 3.2 Perturbation response

For each target, estimate a pseudobulk perturbation response relative to controls using the dataset's prespecified normalized expression layer. The exact transformation is frozen in the manifest. All-finite targets are retained; response significance is not an eligibility criterion.

### 3.3 Fold-pure response space

The response gene universe is defined by identifier intersection and data-quality rules that do not inspect perturbation effects. Within every outer fold:

- nuisance PCA is fitted on training perturbations only;
- all scalers and response transformations are fitted on training rows only;
- held-out target rows do not influence variance selection, rank selection, or scaling;
- the paid response is the residual after projection on the training-fitted nuisance subspace.

Primary nuisance rank is 200, inherited from the development protocol. Locked sensitivity ranks are 50, 100, and 400 where mathematically feasible.

## 4. Representation extraction

### 4.1 Static block (`S`)

Use the checkpoint's gene-token embedding after the same checkpoint-key conversion documented in Round 23. Static vectors are L2-normalized only inside the locked block-normalization pipeline.

### 4.2 Value block (`V`)

Encode the target token at its control-state expression value and compute:

```text
V_g = H00_g(control state) - S_g
```

This separates baseline value/expression encoding from transformer context.

### 4.3 Context block (`C`)

Compute the primary contextual delta:

```text
C_g = D12_g(control state) - H00_g(control state)
```

The target is not set to zero in the primary analysis. Native `mask_value=-1` and zero-value queries are sensitivity analyses only. They are not interpreted as biological knockouts.

### 4.4 Semantic block (`Q`)

GenePT and GO-C are each standardized on training targets, reduced only with training-fitted transforms, block-normalized, and concatenated. The semantic stack is fixed before confirmation.

### 4.5 Block scaling

Every block receives unit expected squared norm on training data before concatenation. No confirmation-label-based weighting is permitted.

## 5. Evaluation tracks

### 5.1 Track A: availability and within-relation retrieval

This track asks whether known partners carry transferable response information.

- Outer unit: held-out target gene.
- Relation partners may remain in the training set by definition.
- Partner reference: equal-weight mean of available training-partner residual responses.
- Null: degree-, baseline-expression-, perturbation-impact-, and coverage-matched random genes; at least 100 locked repetitions.
- Inference: aggregate target-level paired gains to full CORUM bipartite connected components, then bootstrap/sign-flip components.

Claims are restricted to within-relation interpolation.

### 5.2 Track B: relation-component-held-out nested attribution

This is the promotion track.

- Build a bipartite graph of genes and CORUM complexes, retaining all memberships.
- Connected components are indivisible outer groups.
- No training row may share a relation component with a test row.
- Explicit partner-response averaging is prohibited.
- Fit ridge readouts for the nested feature sets below using training data only:

```text
B0 = S
B1 = S + V
B2 = S + V + C
B3 = Q
B4 = Q + S + V
B5 = Q + S + V + C   # primary model
```

The regularization grid is fixed on development data. Inner tuning uses group-aware splits nested within the outer training data.

### 5.3 Cross-cell-line transfer

As a secondary confirmatory analysis, fit the locked readout on HepG2 and evaluate on Jurkat, and vice versa, restricted to the common target and response-gene universe. No target from the test line is used for fitting transformations.

## 6. Primary outcomes

### 6.1 Primary metric

The primary metric is held-out cosine gain over the matched-random response baseline in the nuisance-paid response space.

For nested attribution, the primary estimand is:

```text
Delta_context = cosine_gain(B5) - cosine_gain(B4)
```

computed as a paired target difference and aggregated equally over relation components.

### 6.2 Key secondary metrics

- out-of-sample `R²` / normalized MSE increment of `B5` over `B4`;
- `B2 - B1`, isolating context before semantic priors;
- partner-reference availability over matched random;
- static accessibility over label-permuted and covariance-matched controls;
- positive-component fraction;
- cross-line rank correlation of component effects;
- coverage-performance frontier.

### 6.3 Multiplicity

Only `Delta_context` for `D12-H00`, common support, rank 200, relation-component-held-out, and cosine gain is primary. All other layers, ranks, query values, supports, and metrics are secondary and labeled as such. No per-target or per-complex oracle may promote the project.

## 7. Statistical inference

For each confirmation cell line separately:

1. compute paired target-level model differences;
2. average within the full relation connected component;
3. report the equal-component mean;
4. form a 95% percentile cluster bootstrap interval with 50,000 deterministic resamples;
5. report a two-sided component sign-flip p-value with 50,000 draws;
6. report the number of targets and components and the effective support hash.

A two-line equal-weight summary is descriptive. It does not replace line-specific replication.

## 8. Promotion gate

Tokenizer/adapter/objective training is reopened only when **all** conditions hold:

1. the Track A partner reference beats matched random with a positive 95% component-cluster interval in both HepG2 and Jurkat;
2. static scGPT accessibility beats its dimension/covariance-matched null with a positive interval in both lines;
3. primary `Delta_context = B5-B4` is positive in both lines;
4. the 95% component-cluster interval for `Delta_context` excludes zero in each line;
5. the equal-line point estimate is at least `0.01` cosine-gain units;
6. the sign of `Delta_context` is unchanged on the prespecified rank and control-state bootstrap sensitivities;
7. no result depends on significant-only filtering or model-native-only support;
8. the contamination/data-date registry contains no known direct confirmation-data exposure by the checkpoint.

Failure of any item keeps architecture training frozen. There is no “near miss” promotion based on an oracle layer or pooled p-value.

## 9. Falsification and interpretation

### Gate passes

The supported statement is narrow:

> A locked contextual scGPT transformation contains target-specific perturbation-response information beyond static/value/semantic features under the tested support and data conditions.

The next allowed step is a small adapter or objective-level intervention, not a new tokenizer by default.

### Gate fails

The project closes as a negative contextual-utilization audit while retaining any independently supported availability and static-accessibility claims. The failure does not imply that all scFM contextual states are universally useless; it bounds the tested model, state construction, response task, and support.

## 10. Required artifacts

Every run must emit source and environment provenance; data/checkpoint/asset hashes; control and target-universe hashes; complete relation memberships and component assignments; split and leakage assertions; representation block statistics; per-target predictions and paired differences; cluster-level summaries; a machine-readable promotion verdict; and a truth-boundary statement.
