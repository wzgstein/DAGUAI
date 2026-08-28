from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class PairedClusterResult:
    n_rows: int
    n_clusters: int
    mean_difference: float
    ci_low: float
    ci_high: float
    positive_cluster_fraction: float
    signflip_p_two_sided: float


def _cluster_means(frame: pd.DataFrame, value_col: str, cluster_col: str) -> np.ndarray:
    required = {value_col, cluster_col}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"missing columns: {sorted(missing)}")
    work = frame[[value_col, cluster_col]].dropna()
    if work.empty:
        raise ValueError("no finite clustered observations")
    values = (
        work.groupby(cluster_col, sort=True, observed=True)[value_col]
        .mean()
        .to_numpy(dtype=float)
    )
    if not np.isfinite(values).all():
        raise ValueError("non-finite cluster mean")
    return values


def paired_cluster_inference(
    frame: pd.DataFrame,
    value_col: str,
    cluster_col: str,
    *,
    bootstrap_repetitions: int = 50_000,
    signflip_repetitions: int = 50_000,
    seed: int = 20_260_828,
    chunk_size: int = 5_000,
) -> PairedClusterResult:
    """Equal-cluster bootstrap and sign-flip inference for a paired difference."""

    if bootstrap_repetitions < 1_000 or signflip_repetitions < 1_000:
        raise ValueError("at least 1000 repetitions are required")
    if chunk_size <= 0:
        raise ValueError("chunk_size must be positive")

    values = _cluster_means(frame, value_col, cluster_col)
    n_clusters = int(values.size)
    if n_clusters < 2:
        raise ValueError("at least two clusters are required")
    rng = np.random.default_rng(seed)

    bootstrap = np.empty(bootstrap_repetitions, dtype=float)
    offset = 0
    while offset < bootstrap_repetitions:
        take = min(chunk_size, bootstrap_repetitions - offset)
        indices = rng.integers(0, n_clusters, size=(take, n_clusters))
        bootstrap[offset : offset + take] = values[indices].mean(axis=1)
        offset += take

    observed = float(values.mean())
    sign_statistics = np.empty(signflip_repetitions, dtype=float)
    offset = 0
    while offset < signflip_repetitions:
        take = min(chunk_size, signflip_repetitions - offset)
        signs = rng.choice(np.array([-1.0, 1.0]), size=(take, n_clusters))
        sign_statistics[offset : offset + take] = (signs * values).mean(axis=1)
        offset += take
    p_value = float(
        (np.count_nonzero(np.abs(sign_statistics) >= abs(observed)) + 1)
        / (signflip_repetitions + 1)
    )

    clean_rows = frame[[value_col, cluster_col]].dropna()
    return PairedClusterResult(
        n_rows=int(clean_rows.shape[0]),
        n_clusters=n_clusters,
        mean_difference=observed,
        ci_low=float(np.quantile(bootstrap, 0.025)),
        ci_high=float(np.quantile(bootstrap, 0.975)),
        positive_cluster_fraction=float(np.mean(values > 0)),
        signflip_p_two_sided=p_value,
    )
