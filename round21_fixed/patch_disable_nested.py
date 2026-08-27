#!/usr/bin/env python3
"""Disable PyTorch nested-tensor/MHA fast paths for layer-hook extraction.

The frozen model weights, tokenization, layers, and outputs are unchanged. This
only forces TransformerEncoderLayer hooks to receive ordinary dense tensors;
PyTorch's NestedTensorImpl does not expose sizes required by the audited hook.
"""
from pathlib import Path

SOURCE = Path("round21_fixed/round21_scgpt_contextual.py")
s = SOURCE.read_text()
old = '''    model.to(device)
    model.eval()
    config["loaded_checkpoint_tensors"] = len(compatible)
'''
new = '''    model.to(device)
    model.eval()
    # Hooks on TransformerEncoderLayer outputs require ordinary dense tensors.
    # Disable inference fast paths that convert padded batches to NestedTensor.
    if hasattr(model.transformer_encoder, "enable_nested_tensor"):
        model.transformer_encoder.enable_nested_tensor = False
    if hasattr(model.transformer_encoder, "use_nested_tensor"):
        model.transformer_encoder.use_nested_tensor = False
    if hasattr(torch.backends, "mha") and hasattr(
        torch.backends.mha, "set_fastpath_enabled"
    ):
        torch.backends.mha.set_fastpath_enabled(False)
    config["nested_tensor_fastpath_disabled"] = True
    config["loaded_checkpoint_tensors"] = len(compatible)
'''
if s.count(old) != 1:
    raise RuntimeError(f"nested fast-path patch target count={s.count(old)}")
SOURCE.write_text(s.replace(old, new))
print("disabled PyTorch nested-tensor/MHA fast paths for contextual hooks")
