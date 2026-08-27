#!/usr/bin/env bash
set -euo pipefail

mkdir -p data/embeddings data/scgpt results vendor /tmp/round20_transport

# The pinned scGPT source predates the current PyTorch Transformer internals.
# Pinning torch and NumPy restores its documented-era execution semantics;
# this does not change checkpoint weights or the scientific evaluation.
python -m pip install --upgrade pip
python -m pip install 'numpy<2'
python -m pip install --index-url https://download.pytorch.org/whl/cpu 'torch==2.2.2'
python -m pip install 'numpy<2' anndata h5py pandas scipy scikit-learn matplotlib tabulate scanpy awscli

git clone --filter=blob:none https://github.com/bowang-lab/scGPT.git vendor/scGPT
git -C vendor/scGPT checkout cebd6fae655b9c585a4807daa3ac31bb764f06b4
cat > vendor/scGPT/scgpt/__init__.py <<'PY'
import logging, sys
logger = logging.getLogger('scGPT')
if not logger.handlers:
    logger.setLevel(logging.INFO)
    logger.addHandler(logging.StreamHandler(sys.stdout))
PY
cat > vendor/scGPT/scgpt/model/__init__.py <<'PY'
from .model import TransformerModel
PY

# Reconstruct the audited Round 20 evaluation base.
for i in 00 01 02 03; do
  curl -L --fail --retry 8 \
    "https://raw.githubusercontent.com/wzgstein/DAGUAI/9cc0bc29e73f6ce926b788d76cc08b36deb36d6e/round20_transport/part${i}" \
    -o "/tmp/round20_transport/part${i}"
done
cat /tmp/round20_transport/part00 /tmp/round20_transport/part01 /tmp/round20_transport/part02 /tmp/round20_transport/part03 > /tmp/round20_source.b64
test "$(wc -c < /tmp/round20_source.b64)" -eq 15144
echo 'd89606a1e7c3a64053b2f4b0cc19a6ef73d5e48ab16c433041b8e5e4593eb8d1  /tmp/round20_source.b64' | sha256sum -c -
base64 --decode /tmp/round20_source.b64 | gzip -dc > round21_fixed/round20_base.py

# Reconstruct the contextual source already committed in audited chunks.
cat tmp_round21/contextual_transport/part* | base64 --decode | gzip -dc > round21_fixed/round21_scgpt_contextual.py
echo '8ecb07bcb7f48d7b1bab4f9e05642e7b5a616bd114fd41812e13de5f25f2a171  round21_fixed/round20_base.py' | sha256sum -c -
echo '128dd96a1748005717d1399ba94ffde4bfd41c435083c88b4e512c42f1dd0280  round21_fixed/round21_scgpt_contextual.py' | sha256sum -c -

# Backed H5AD/h5py permits only one fancy index in a single read.  Preserve the
# exact rows/columns while applying them sequentially.
python - <<'PY'
from pathlib import Path
p = Path('round21_fixed/round20_base.py')
s = p.read_text()
old = '''def matrix_slice(adata: ad.AnnData, rows: np.ndarray, cols: np.ndarray) -> np.ndarray:\n    x = adata[rows, cols].X\n    if sparse.issparse(x):\n        x = x.toarray()\n    return np.asarray(x, dtype=np.float32)\n'''
new = '''def matrix_slice(adata: ad.AnnData, rows: np.ndarray, cols: np.ndarray) -> np.ndarray:\n    x = adata[rows, :].X\n    if sparse.issparse(x):\n        x = x[:, cols].toarray()\n    else:\n        x = np.asarray(x)[:, cols]\n    return np.asarray(x, dtype=np.float32)\n'''
if s.count(old) != 1:
    raise RuntimeError(f'expected one matrix_slice implementation, found {s.count(old)}')
p.write_text(s.replace(old, new))
print('patched backed-H5AD sequential slicing')
PY
python -m py_compile round21_fixed/round20_base.py round21_fixed/round21_scgpt_contextual.py

# Real control cells and the frozen Round 19/20 evaluation objects.
curl -L --fail --retry 12 --retry-all-errors --retry-delay 8 --connect-timeout 30 --speed-limit 1024 --speed-time 180 \
  https://storage.googleapis.com/scperteval/processed/replogle22k562_processed_complete.h5ad \
  -o data/K562_single_complete.h5ad
curl -L --fail --retry 12 --retry-all-errors --retry-delay 8 --connect-timeout 30 --speed-limit 1024 --speed-time 180 \
  https://storage.googleapis.com/scperteval/processed/replogle22rpe1_processed_complete.h5ad \
  -o data/RPE1_single_complete.h5ad
curl -L --fail --retry 12 --retry-all-errors --retry-delay 8 \
  https://ndownloader.figshare.com/files/35773217 -o data/K562_bulk.h5ad
curl -L --fail --retry 12 --retry-all-errors --retry-delay 8 \
  https://ndownloader.figshare.com/files/35775512 -o data/RPE1_bulk.h5ad
curl -L --fail --retry 8 --retry-all-errors \
  https://raw.githubusercontent.com/broadinstitute/2022_PERISCOPE/17d345046beee46c8fa8601f970b66651232e7ed/common_files/CORUM_humanComplexes.txt \
  -o data/CORUM_humanComplexes.txt
aws s3 sync --no-sign-request s3://czi-scgenept-public/models/pretrained/scgpt data/scgpt
aws s3 cp --no-sign-request \
  's3://czi-scgenept-public/models/gene_embeddings/NCBI+UniProt_embeddings-gpt3.5-ada.pkl' \
  'data/embeddings/NCBI+UniProt_embeddings-gpt3.5-ada.pkl'
aws s3 cp --no-sign-request \
  's3://czi-scgenept-public/models/gene_embeddings/GO_C_gene_embeddings-gpt3.5-ada-concat.pickle' \
  'data/embeddings/GO_C_gene_embeddings-gpt3.5-ada-concat.pickle'
sha256sum data/*.h5ad data/CORUM_humanComplexes.txt data/scgpt/* data/embeddings/* | tee results/INPUT_SHA256.txt

# Frozen reliability rule from Round 19.
python - <<'PY'
import json
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
report = {}
for source, target in [
    (Path('data/K562_bulk.h5ad'), Path('data/K562_bulk_significant.h5ad')),
    (Path('data/RPE1_bulk.h5ad'), Path('data/RPE1_bulk_significant.h5ad')),
]:
    a = ad.read_h5ad(source)
    x = a.X.toarray() if sparse.issparse(a.X) else np.asarray(a.X)
    finite = np.isfinite(x).all(axis=1)
    p = pd.to_numeric(a.obs['energy_test_p_value'], errors='coerce').to_numpy(float)
    significant = np.isfinite(p) & (p < 1e-3)
    keep = finite & significant
    report[source.name] = {
        'rows_before': int(a.n_obs),
        'finite_rows': int(finite.sum()),
        'significant_rows': int(significant.sum()),
        'rows_kept': int(keep.sum()),
        'rule': 'complete finite response row and energy_test_p_value < 0.001',
    }
    if keep.sum() < 200:
        raise RuntimeError(f'too few reliable rows: {source} {keep.sum()}')
    a[keep].copy().write_h5ad(target, compression='gzip')
Path('results/PSEUDOBULK_FILTER.json').write_text(json.dumps(report, indent=2))
print(json.dumps(report, indent=2))
PY
rm data/K562_bulk.h5ad data/RPE1_bulk.h5ad

export PYTHONPATH="$PWD/vendor/scGPT"
/usr/bin/time -v python round21_fixed/round21_scgpt_contextual.py \
  --round20-base round21_fixed/round20_base.py \
  --k562-bulk data/K562_bulk_significant.h5ad \
  --rpe1-bulk data/RPE1_bulk_significant.h5ad \
  --k562-single data/K562_single_complete.h5ad \
  --rpe1-single data/RPE1_single_complete.h5ad \
  --corum data/CORUM_humanComplexes.txt \
  --embedding-dir data/embeddings \
  --scgpt-dir data/scgpt \
  --out results \
  --control-cells 96 \
  --max-length 384 \
  --batch-size 4 \
  --min-occurrences 3 \
  --response-genes 2000 \
  --variance-sample-rows 1200 \
  --folds 5 \
  --random-reps 20 \
  --embedding-pcs 64 2>&1 | tee results/RUN_LOG.txt
