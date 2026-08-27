#!/usr/bin/env python3
"""Apply data-interface-only compatibility fixes to the audited Round 21 source.

The patch does not alter cell selection, gene selection, model weights, layer
selection, response residualization, random controls, or evaluation metrics.
It only:
1. converts categorical control labels through pandas StringDtype;
2. bypasses anndata/scipy backed sparse slicing with an explicit float32
   CSR/CSC subset reader.
"""
from pathlib import Path

SOURCE = Path("round21_fixed/round21_scgpt_contextual.py")
s = SOURCE.read_text()

old_control = '    s = values.fillna("").astype(str).str.lower().str.strip()\n'
new_control = '    s = values.astype("string").fillna("").str.lower().str.strip()\n'
if s.count(old_control) != 1:
    raise RuntimeError(f"categorical patch target count={s.count(old_control)}")
s = s.replace(old_control, new_control)

insertion = "\n\ndef load_control_subset(\n"
helper = '''

def read_h5ad_matrix_subset(
    path: Path,
    node_path: str,
    row_indices: np.ndarray,
    col_indices: np.ndarray,
):
    """Read a small row/column subset with explicit float32 sparse dtype."""
    import h5py

    rows = np.asarray(row_indices, dtype=np.int64)
    cols = np.asarray(col_indices, dtype=np.int64)
    if rows.ndim != 1 or cols.ndim != 1:
        raise ValueError("row and column indices must be one-dimensional")
    if np.any(np.diff(rows) < 0):
        raise ValueError("row indices must be sorted")

    def decoded(value):
        if isinstance(value, bytes):
            return value.decode("utf-8")
        return str(value)

    with h5py.File(path, "r") as handle:
        node = handle[node_path]
        if isinstance(node, h5py.Dataset):
            dense = np.asarray(node[rows, :], dtype=np.float32)
            return sparse.csr_matrix(dense[:, cols], dtype=np.float32)

        encoding = decoded(node.attrs.get("encoding-type", ""))
        shape_raw = node.attrs.get("shape")
        if shape_raw is None and "shape" in node:
            shape_raw = node["shape"][...]
        n_rows, n_cols = tuple(int(x) for x in np.asarray(shape_raw).ravel())

        if encoding == "csr_matrix":
            indptr = np.asarray(node["indptr"][...], dtype=np.int64)
            out_indptr = np.zeros(rows.size + 1, dtype=np.int64)
            data_parts = []
            index_parts = []
            cursor = 0
            for j, row in enumerate(rows.tolist()):
                start, stop = int(indptr[row]), int(indptr[row + 1])
                if stop > start:
                    data_parts.append(
                        np.asarray(node["data"][start:stop], dtype=np.float32)
                    )
                    index_parts.append(
                        np.asarray(node["indices"][start:stop], dtype=np.int64)
                    )
                    cursor += stop - start
                out_indptr[j + 1] = cursor
            data = (
                np.concatenate(data_parts)
                if data_parts
                else np.empty(0, dtype=np.float32)
            )
            indices = (
                np.concatenate(index_parts)
                if index_parts
                else np.empty(0, dtype=np.int64)
            )
            full = sparse.csr_matrix(
                (data, indices, out_indptr),
                shape=(rows.size, n_cols),
                dtype=np.float32,
            )
            return full[:, cols].tocsr().astype(np.float32)

        if encoding == "csc_matrix":
            indptr = np.asarray(node["indptr"][...], dtype=np.int64)
            row_map = np.full(n_rows, -1, dtype=np.int64)
            row_map[rows] = np.arange(rows.size, dtype=np.int64)
            data_parts = []
            row_parts = []
            col_parts = []
            for output_col, original_col in enumerate(cols.tolist()):
                start = int(indptr[original_col])
                stop = int(indptr[original_col + 1])
                if stop <= start:
                    continue
                original_rows = np.asarray(
                    node["indices"][start:stop], dtype=np.int64
                )
                local_rows = row_map[original_rows]
                keep = local_rows >= 0
                if not np.any(keep):
                    continue
                data_parts.append(
                    np.asarray(node["data"][start:stop], dtype=np.float32)[keep]
                )
                row_parts.append(local_rows[keep])
                col_parts.append(
                    np.full(int(np.sum(keep)), output_col, dtype=np.int64)
                )
            if data_parts:
                return sparse.coo_matrix(
                    (
                        np.concatenate(data_parts),
                        (np.concatenate(row_parts), np.concatenate(col_parts)),
                    ),
                    shape=(rows.size, cols.size),
                    dtype=np.float32,
                ).tocsr()
            return sparse.csr_matrix((rows.size, cols.size), dtype=np.float32)

        raise ValueError(
            f"unsupported H5AD matrix encoding {encoding!r} at {node_path}"
        )
'''
if s.count(insertion) != 1:
    raise RuntimeError(f"helper insertion point count={s.count(insertion)}")
s = s.replace(insertion, helper + insertion)

old_subset = '''    sub = a[control_idx, np.asarray(cols, dtype=int)].to_memory()
    raw_shape = [int(a.n_obs), int(a.n_vars)]
    obs_columns = [str(x) for x in a.obs.columns]
    var_columns = [str(x) for x in a.var.columns]
    if getattr(a, "file", None) is not None:
        a.file.close()
    x, matrix_meta = select_source_matrix(sub)
'''
new_subset = '''    raw_shape = [int(a.n_obs), int(a.n_vars)]
    obs_columns = [str(x) for x in a.obs.columns]
    var_columns = [str(x) for x in a.var.columns]
    layer_name = next(
        (
            candidate
            for candidate in ("counts", "raw_counts")
            if candidate in a.layers
        ),
        None,
    )
    node_path = f"layers/{layer_name}" if layer_name else "X"
    layer_keys = [str(x) for x in a.layers.keys()]
    if getattr(a, "file", None) is not None:
        a.file.close()
    raw_matrix = read_h5ad_matrix_subset(
        path,
        node_path,
        control_idx,
        np.asarray(cols, dtype=np.int64),
    )
    temporary = ad.AnnData(X=raw_matrix)
    x, matrix_meta = select_source_matrix(temporary)
    matrix_meta["source"] = node_path
'''
if s.count(old_subset) != 1:
    raise RuntimeError(f"subset patch target count={s.count(old_subset)}")
s = s.replace(old_subset, new_subset)

old_layers = '        "layers": list(sub.layers.keys()),\n'
new_layers = '        "layers": layer_keys,\n'
if s.count(old_layers) != 1:
    raise RuntimeError(f"layers patch target count={s.count(old_layers)}")
s = s.replace(old_layers, new_layers)

SOURCE.write_text(s)
print("patched categorical labels and explicit float32 CSR/CSC reader")
