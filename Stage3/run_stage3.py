#!/usr/bin/env python3
"""
Stage 3 protein profile summarization runner.

Purpose:
- Load a 25-protein test set.
- Send one full protein JSON object at a time to a local Ollama model.
- Use a fixed prompt and fresh model request for every protein.
- Save raw and parsed outputs after every protein.
- Perform only minimal structural validation.

This script deliberately does NOT validate biological correctness and does NOT
ask the model to choose a final annotation. Outputs are information summaries
for manual inspection.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
import time
from pathlib import Path
from typing import Any, Iterable

import requests


REQUIRED_FIELDS = [
    "protein_id",
    "protein_name",
    "organism",
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
]

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

STRING_FIELDS = [
    "protein_id",
    "protein_name",
    "organism",
    "source_status",
    "overall_profile_summary",
    "identity_summary",
    "reported_function_summary",
]


def now_iso() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)


def append_jsonl(path: Path, record: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as f:
        f.write(json.dumps(record, ensure_ascii=False) + "\n")


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    records: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as f:
        for line_no, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError:
                print(f"Warning: could not parse JSONL line {line_no} in {path}", file=sys.stderr)
    return records


def remove_empty_values(obj: Any) -> Any:
    """
    Optional context reduction.

    Default behavior is to pass the full protein JSON unchanged. Use --strip-empty
    only if your model has context-window issues. This function removes only
    structurally empty values: None, empty strings, empty lists, and empty dicts.
    """
    if isinstance(obj, dict):
        cleaned = {}
        for key, value in obj.items():
            new_value = remove_empty_values(value)
            if new_value in (None, "", [], {}):
                continue
            cleaned[key] = new_value
        return cleaned

    if isinstance(obj, list):
        cleaned_list = [remove_empty_values(item) for item in obj]
        return [item for item in cleaned_list if item not in (None, "", [], {})]

    return obj


def build_prompt(template: str, protein: dict[str, Any], strip_empty: bool = False) -> str:
    prompt_protein = remove_empty_values(protein) if strip_empty else protein
    protein_json = json.dumps(prompt_protein, ensure_ascii=False, indent=2)
    return template.replace("{protein_json}", protein_json)


def extract_first_json_object(text: str) -> str | None:
    """
    Best-effort fallback if a model wraps JSON in extra text despite instructions.
    This is not a biological retry loop; it only tries to recover a syntactically
    valid JSON object from the raw response.
    """
    start = text.find("{")
    if start == -1:
        return None

    depth = 0
    in_string = False
    escape = False

    for idx in range(start, len(text)):
        char = text[idx]

        if in_string:
            if escape:
                escape = False
            elif char == "\\":
                escape = True
            elif char == '"':
                in_string = False
            continue

        if char == '"':
            in_string = True
        elif char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                return text[start : idx + 1]

    return None


def parse_llm_json(raw: str) -> tuple[dict[str, Any] | None, str | None]:
    try:
        parsed = json.loads(raw)
        if isinstance(parsed, dict):
            return parsed, None
        return None, "Parsed JSON is not an object."
    except json.JSONDecodeError as first_error:
        candidate = extract_first_json_object(raw)
        if not candidate:
            return None, f"JSON parse failed: {first_error}"

        try:
            parsed = json.loads(candidate)
            if isinstance(parsed, dict):
                return parsed, "Recovered JSON object from extra response text."
            return None, "Recovered JSON was not an object."
        except json.JSONDecodeError as second_error:
            return None, f"JSON parse failed: {first_error}; recovery failed: {second_error}"


def validate_summary(parsed: dict[str, Any] | None, accession: str) -> list[str]:
    warnings: list[str] = []

    if parsed is None:
        return ["Response could not be parsed as a JSON object."]

    for field in REQUIRED_FIELDS:
        if field not in parsed:
            warnings.append(f"Missing required field: {field}")

    protein_id = parsed.get("protein_id")
    if protein_id != accession:
        warnings.append(f"protein_id mismatch: expected {accession!r}, got {protein_id!r}")

    for field in LIST_FIELDS:
        if field in parsed and not isinstance(parsed[field], list):
            warnings.append(f"Field {field!r} should be a list.")

    for field in STRING_FIELDS:
        if field in parsed and not isinstance(parsed[field], str):
            warnings.append(f"Field {field!r} should be a string.")
        elif field in parsed and not parsed[field].strip():
            warnings.append(f"Field {field!r} is empty.")

    # Important conceptual guard: this summarization version should not ask
    # the model to produce a final annotation decision. If a model adds this
    # anyway, flag it for reviewer awareness rather than failing the run.
    discouraged_fields = ["recommended_annotation", "confidence", "manual_review_recommended"]
    for field in discouraged_fields:
        if field in parsed:
            warnings.append(
                f"Unexpected field {field!r}: this workflow summarizes protein information, "
                "not final annotation decisions."
            )

    return warnings


def call_ollama(
    *,
    url: str,
    model: str,
    prompt: str,
    schema: dict[str, Any] | None,
    format_mode: str,
    temperature: float,
    timeout: int,
    num_ctx: int | None,
) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "model": model,
        "prompt": prompt,
        "stream": False,
        "think": False,
        "options": {
            "temperature": temperature,
        },
    }

    if num_ctx is not None:
        payload["options"]["num_ctx"] = num_ctx

    if format_mode == "schema":
        if schema is None:
            raise ValueError("format_mode='schema' requires a schema.")
        payload["format"] = schema
    elif format_mode == "json":
        payload["format"] = "json"
    elif format_mode == "none":
        pass
    else:
        raise ValueError(f"Unknown format mode: {format_mode}")

    response = requests.post(url, json=payload, timeout=timeout)
    response.raise_for_status()
    return response.json()


def selected_proteins(
    proteins: list[dict[str, Any]],
    groups: set[str] | None,
    limit: int | None,
    only_accessions: set[str] | None,
) -> list[dict[str, Any]]:
    selected = proteins

    if groups:
        selected = [p for p in selected if p.get("_test_group") in groups]

    if only_accessions:
        selected = [p for p in selected if p.get("accession") in only_accessions]

    if limit is not None:
        selected = selected[:limit]

    return selected


def completed_accessions(results_path: Path, rerun_failed: bool) -> set[str]:
    done: set[str] = set()
    for record in read_jsonl(results_path):
        accession = record.get("accession")
        status = record.get("status")
        if not accession:
            continue

        if rerun_failed and status not in {"success", "success_with_warnings"}:
            continue

        done.add(accession)

    return done


def write_current_json(results_path: Path, current_json_path: Path) -> None:
    records = read_jsonl(results_path)
    write_json(current_json_path, records)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run Stage 3 protein profile summarization over a protein test set with Ollama."
    )
    parser.add_argument("--input", type=Path, default=Path("test_proteins.json"))
    parser.add_argument("--prompt-template", type=Path, default=Path("prompt_template.txt"))
    parser.add_argument("--schema", type=Path, default=Path("schema/protein_profile_summary.schema.json"))
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    parser.add_argument("--model", default="llama3.1:8b")
    parser.add_argument("--ollama-url", default="http://localhost:11435/api/generate")
    parser.add_argument("--format-mode", choices=["schema", "json", "none"], default="schema")
    parser.add_argument("--temperature", type=float, default=0.0)
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--num-ctx", type=int, default=None, help="Optional Ollama num_ctx override.")
    parser.add_argument("--limit", type=int, default=None, help="Run only the first N selected proteins.")
    parser.add_argument(
        "--group",
        action="append",
        default=None,
        help="Only run proteins from this _test_group. Can be used multiple times.",
    )
    parser.add_argument(
        "--accession",
        action="append",
        default=None,
        help="Only run this accession. Can be used multiple times.",
    )
    parser.add_argument("--strip-empty", action="store_true", help="Remove null/empty values from prompt JSON.")
    parser.add_argument("--save-prompts", action="store_true", help="Save the exact prompt sent for each protein.")
    parser.add_argument("--dry-run", action="store_true", help="Build prompts and validate inputs, but do not call Ollama.")
    parser.add_argument("--rerun-failed", action="store_true", help="Skip successful records but rerun failed records.")
    parser.add_argument(
        "--no-keep-input",
        action="store_true",
        help="Do not store the full input protein inside each result record.",
    )
    args = parser.parse_args()

    proteins = load_json(args.input)
    if not isinstance(proteins, list):
        raise ValueError(f"Expected a list of proteins in {args.input}")

    template = args.prompt_template.read_text(encoding="utf-8")
    schema = load_json(args.schema) if args.format_mode == "schema" else None

    output_dir = args.output_dir
    raw_dir = output_dir / "raw_responses"
    prompt_dir = output_dir / "prompts"
    results_jsonl = output_dir / "stage3_results.jsonl"
    current_json = output_dir / "stage3_results.current.json"
    output_dir.mkdir(parents=True, exist_ok=True)
    raw_dir.mkdir(parents=True, exist_ok=True)
    if args.save_prompts or args.dry_run:
        prompt_dir.mkdir(parents=True, exist_ok=True)

    groups = set(args.group) if args.group else None
    accessions = set(args.accession) if args.accession else None
    todo = selected_proteins(proteins, groups=groups, limit=args.limit, only_accessions=accessions)

    done = completed_accessions(results_jsonl, rerun_failed=args.rerun_failed)

    print(f"Input proteins: {len(proteins)}")
    print(f"Selected proteins: {len(todo)}")
    print(f"Already completed/skipped: {len(done)}")
    print(f"Output JSONL: {results_jsonl}")

    for index, protein in enumerate(todo, start=1):
        accession = protein.get("accession")
        if not accession:
            print(f"[{index}/{len(todo)}] Skipping protein without accession.", file=sys.stderr)
            continue

        if accession in done:
            print(f"[{index}/{len(todo)}] Skipping completed accession {accession}")
            continue

        group = protein.get("_test_group")
        print(f"[{index}/{len(todo)}] Processing {accession} ({group})")

        prompt = build_prompt(template, protein, strip_empty=args.strip_empty)

        if args.save_prompts or args.dry_run:
            (prompt_dir / f"{accession}.txt").write_text(prompt, encoding="utf-8")

        result: dict[str, Any] = {
            "accession": accession,
            "group": group,
            "status": "pending",
            "run_timestamp_utc": now_iso(),
            "model": args.model,
            "format_mode": args.format_mode,
            "temperature": args.temperature,
            "input_protein": None if args.no_keep_input else protein,
            "llm_response_raw": None,
            "llm_response_json": None,
            "validation_warnings": [],
            "parse_note": None,
            "error": None,
        }

        if args.dry_run:
            result["status"] = "dry_run"
            result["validation_warnings"] = ["Dry run only; Ollama was not called."]
            append_jsonl(results_jsonl, result)
            write_current_json(results_jsonl, current_json)
            continue

        try:
            response_payload = call_ollama(
                url=args.ollama_url,
                model=args.model,
                prompt=prompt,
                schema=schema,
                format_mode=args.format_mode,
                temperature=args.temperature,
                timeout=args.timeout,
                num_ctx=args.num_ctx,
            )

            raw = response_payload.get("response", "")
            if not isinstance(raw, str):
                raw = json.dumps(raw, ensure_ascii=False)

            (raw_dir / f"{accession}.txt").write_text(raw, encoding="utf-8")

            parsed, parse_note = parse_llm_json(raw)
            warnings = validate_summary(parsed, accession)

            result["llm_response_raw"] = raw
            result["llm_response_json"] = parsed
            result["validation_warnings"] = warnings
            result["parse_note"] = parse_note

            if parsed is None:
                result["status"] = "json_parse_failed"
            elif warnings:
                result["status"] = "success_with_warnings"
            else:
                result["status"] = "success"

            # Keep selected Ollama metadata for reproducibility without bloating the result.
            for key in [
                "created_at",
                "done",
                "total_duration",
                "load_duration",
                "prompt_eval_count",
                "prompt_eval_duration",
                "eval_count",
                "eval_duration",
            ]:
                if key in response_payload:
                    result[f"ollama_{key}"] = response_payload[key]

        except requests.exceptions.ConnectionError as exc:
            result["status"] = "error"
            result["error"] = (
                "Could not connect to Ollama. Check that Ollama is running and that "
                f"--ollama-url is correct. Details: {exc}"
            )
        except requests.exceptions.Timeout as exc:
            result["status"] = "error"
            result["error"] = f"Ollama request timed out after {args.timeout} seconds: {exc}"
        except Exception as exc:
            result["status"] = "error"
            result["error"] = repr(exc)

        append_jsonl(results_jsonl, result)
        write_current_json(results_jsonl, current_json)

        if result["status"] == "error":
            print(f"  ERROR: {result['error']}", file=sys.stderr)
        elif result["validation_warnings"]:
            print(f"  Completed with warnings: {result['validation_warnings']}")
        else:
            print("  Completed successfully.")

        # Small pause prevents accidental hammering if Ollama is running on a small laptop.
        time.sleep(0.2)

    write_current_json(results_jsonl, current_json)
    print(f"Done. Combined JSON written to {current_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
