#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from relation_audit.manifest import load_and_validate_manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument(
        "--schema",
        type=Path,
        default=Path("schema/experiment_manifest.schema.json"),
    )
    args = parser.parse_args()
    manifest = load_and_validate_manifest(args.manifest, args.schema)
    print(f"valid manifest: {manifest['experiment_id']}")


if __name__ == "__main__":
    main()
