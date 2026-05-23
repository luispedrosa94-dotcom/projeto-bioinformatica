#!/usr/bin/env python3
"""
Create a supervisor-friendly manual review CSV from Stage 3 protein profile summaries.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            if line.strip():
                records.append(json.loads(line))
    return records


def as_text(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, list):
        return "\n".join(str(item) for item in value)
    if isinstance(value, dict):
        return json.dumps(value, ensure_ascii=False)
    return str(value)


def get_nested(record: dict[str, Any], path: list[str], default: Any = "") -> Any:
    cur: Any = record
    for key in path:
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return cur


def main() -> int:
    parser = argparse.ArgumentParser(description="Export protein profile summaries to a manual review CSV.")
    parser.add_argument("--results", type=Path, default=Path("outputs/stage3_results.jsonl"))
    parser.add_argument("--output", type=Path, default=Path("outputs/review_sheet.csv"))
    args = parser.parse_args()

    records = read_jsonl(args.results)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    fields = [
        "accession",
        "group",
        "status",
        "input_protein_name",
        "input_reviewed_status",
        "summary_protein_name",
        "summary_organism",
        "source_status",
        "overall_profile_summary",
        "identity_summary",
        "reported_function_summary",
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
        "validation_warnings",
        "completeness_0_3",
        "faithfulness_0_3",
        "clarity_0_3",
        "evidence_separation_0_3",
        "caution_0_3",
        "reviewer_missing_information",
        "reviewer_overclaims_or_errors",
        "reviewer_notes",
    ]

    with args.output.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()

        for record in records:
            summary = record.get("llm_response_json") or {}
            input_protein = record.get("input_protein") or {}

            row = {
                "accession": record.get("accession", ""),
                "group": record.get("group", ""),
                "status": record.get("status", ""),
                "input_protein_name": get_nested(input_protein, ["identity", "protein_name"]),
                "input_reviewed_status": get_nested(input_protein, ["identity", "reviewed_status"]),
                "summary_protein_name": summary.get("protein_name", ""),
                "summary_organism": summary.get("organism", ""),
                "source_status": summary.get("source_status", ""),
                "overall_profile_summary": summary.get("overall_profile_summary", ""),
                "identity_summary": summary.get("identity_summary", ""),
                "reported_function_summary": summary.get("reported_function_summary", ""),
                "go_annotation_summary": as_text(summary.get("go_annotation_summary")),
                "enzyme_and_reaction_summary": as_text(summary.get("enzyme_and_reaction_summary")),
                "domain_family_and_feature_summary": as_text(summary.get("domain_family_and_feature_summary")),
                "pathway_and_context_summary": as_text(summary.get("pathway_and_context_summary")),
                "tool_prediction_summary": as_text(summary.get("tool_prediction_summary")),
                "strong_or_curated_information": as_text(summary.get("strong_or_curated_information")),
                "weak_predicted_or_indirect_information": as_text(summary.get("weak_predicted_or_indirect_information")),
                "conflicting_or_inconsistent_information": as_text(summary.get("conflicting_or_inconsistent_information")),
                "missing_or_limited_information": as_text(summary.get("missing_or_limited_information")),
                "review_notes": as_text(summary.get("review_notes")),
                "validation_warnings": as_text(record.get("validation_warnings")),
                "completeness_0_3": "",
                "faithfulness_0_3": "",
                "clarity_0_3": "",
                "evidence_separation_0_3": "",
                "caution_0_3": "",
                "reviewer_missing_information": "",
                "reviewer_overclaims_or_errors": "",
                "reviewer_notes": "",
            }
            writer.writerow(row)

    print(f"Wrote review sheet: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
