from __future__ import annotations

import numpy as np


def cosine_similarity_rows(observed: np.ndarray, predicted: np.ndarray) -> np.ndarray:
    observed = np.asarray(observed, dtype=float)
    predicted = np.asarray(predicted, dtype=float)
    if observed.shape != predicted.shape or observed.ndim != 2:
        raise ValueError("observed and predicted must be equally shaped 2D arrays")
    numerator = np.einsum("ij,ij->i", observed, predicted)
    denominator = np.linalg.norm(observed, axis=1) * np.linalg.norm(predicted, axis=1)
    output = np.full(observed.shape[0], np.nan, dtype=float)
    valid = np.isfinite(numerator) & np.isfinite(denominator) & (denominator > 0)
    output[valid] = numerator[valid] / denominator[valid]
    return output


def normalized_mse_rows(observed: np.ndarray, predicted: np.ndarray) -> np.ndarray:
    observed = np.asarray(observed, dtype=float)
    predicted = np.asarray(predicted, dtype=float)
    if observed.shape != predicted.shape or observed.ndim != 2:
        raise ValueError("observed and predicted must be equally shaped 2D arrays")
    numerator = np.mean((observed - predicted) ** 2, axis=1)
    denominator = np.mean(observed**2, axis=1)
    output = np.full(observed.shape[0], np.nan, dtype=float)
    valid = np.isfinite(numerator) & np.isfinite(denominator) & (denominator > 0)
    output[valid] = numerator[valid] / denominator[valid]
    return output
