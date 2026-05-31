#!/usr/bin/env python3
"""
Summarize Stage 3 protein profile summarization run status.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Any


LIST_FIELDS = [
    "go_annotation_summary",
    "enzyme_and_reaction_summary",
    "domain_family_and_feature_summary",
    "pathway_and_context_summary",
    "tool_prediction_summary",
    "strong_or_curated_information",
    "weak_predicted_or_indirect_information",
    "conflicting_or_inconsistent_information",
    "missing_or_limited_information",
    "review_notes",
]


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    records = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            if line.strip():
                records.append(json.loads(line))
    return records


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize Stage 3 JSONL results.")
    parser.add_argument("--results", type=Path, default=Path("outputs/stage3_results.jsonl"))
    args = parser.parse_args()

    records = read_jsonl(args.results)
    print(f"Records: {len(records)}")

    status_counts = Counter(r.get("status", "<missing>") for r in records)
    print("\nStatus counts:")
    for status, count in status_counts.items():
        print(f"  {status}: {count}")

    group_counts = Counter(r.get("group", "<missing>") for r in records)
    print("\nGroup counts:")
    for group, count in group_counts.items():
        print(f"  {group}: {count}")

    warning_counts = Counter(bool(r.get("validation_warnings")) for r in records)
    print("\nValidation warnings:")
    print(f"  without warnings: {warning_counts[False]}")
    print(f"  with warnings: {warning_counts[True]}")

    parsed = [r.get("llm_response_json") or {} for r in records]
    print("\nAverage number of bullets per summary section:")
    for field in LIST_FIELDS:
        counts = [len(s.get(field) or []) for s in parsed if isinstance(s, dict)]
        avg = sum(counts) / len(counts) if counts else 0
        print(f"  {field}: {avg:.1f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
