#!/usr/bin/env python3
"""
Inspect the Stage 3 protein test set.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description="Check accessions and group counts in test_proteins.json.")
    parser.add_argument("--input", type=Path, default=Path("test_proteins.json"))
    args = parser.parse_args()

    proteins = json.loads(args.input.read_text(encoding="utf-8"))
    counts = Counter(p.get("_test_group", "<missing>") for p in proteins)

    print(f"Total proteins: {len(proteins)}")
    print("\nGroup counts:")
    for group, count in counts.items():
        print(f"  {group}: {count}")

    print("\nAccessions by group:")
    grouped = defaultdict(list)
    for p in proteins:
        grouped[p.get("_test_group", "<missing>")].append(p.get("accession", "<missing>"))

    for group, accessions in grouped.items():
        print(f"  {group}: {', '.join(accessions)}")

    missing_accession = [i for i, p in enumerate(proteins) if not p.get("accession")]
    if missing_accession:
        print(f"\nWARNING: missing accession at indexes: {missing_accession}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
