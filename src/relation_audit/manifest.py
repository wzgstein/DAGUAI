from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import jsonschema


def sha256_lines(values: list[str]) -> str:
    payload = "\n".join(sorted(dict.fromkeys(values))).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def load_and_validate_manifest(
    path: str | Path,
    schema_path: str | Path,
) -> dict[str, Any]:
    manifest = json.loads(Path(path).read_text(encoding="utf-8"))
    schema = json.loads(Path(schema_path).read_text(encoding="utf-8"))
    jsonschema.Draft202012Validator(schema).validate(manifest)
    return manifest
