"""
Inspect the normalized output.

Reads the JSON files produced by normalize.py and prints a friendly summary
of the data — counts by tool, by annotation type, by confidence level, plus
a few example proteins to explore.

The JSON files themselves are already human-readable, so this script no
longer exports CSV — open the JSON files directly in a text editor to see
the raw content.

Usage:
    python scripts/peek.py                          # summary
    python scripts/peek.py --accession A0B9K2       # focus on one protein
    python scripts/peek.py --outputs ../outputs     # custom outputs directory
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--outputs",
        default="outputs",
        help="Directory where annotations.json and proteins.json were written",
    )
    parser.add_argument(
        "--accession",
        default=None,
        help="Show all annotations for one specific UniProt accession (e.g. A0B9K2)",
    )
    args = parser.parse_args()

    outputs = Path(args.outputs)
    annot_path = outputs / "annotations.json"
    prot_path = outputs / "proteins.json"

    if not annot_path.exists():
        print(f"ERROR: {annot_path} not found.", file=sys.stderr)
        print("Did you run scripts/normalize.py first?", file=sys.stderr)
        sys.exit(1)

    annotations = pd.read_json(annot_path, orient="records")
    proteins = pd.read_json(prot_path, orient="records") if prot_path.exists() else None

    # ---- Mode 1: focus on a single protein ----
    if args.accession:
        acc = args.accession.upper()
        sub = annotations[annotations["uniprot_accession"] == acc]
        if sub.empty:
            print(f"No records for {acc}.")
            return

        if proteins is not None:
            protein_row = proteins[proteins["uniprot_accession"] == acc]
            print(f"\n=== Protein {acc} ===")
            if not protein_row.empty:
                row = protein_row.iloc[0]
                print(f"Original ID:       {row['original_id']}")
                print(f"DB source:         {row['db_source']}")
                print(f"Entry name:        {row['entry_name']}")
                print(f"Poorly annotated:  {row['in_poorly_annotated_subset']}")

        print(f"\n{len(sub)} annotation records:\n")
        cols = ["source_tool", "annotation_type", "value", "label",
                "score", "score_type", "confidence_level"]
        with pd.option_context("display.max_colwidth", 70, "display.width", 200):
            print(sub[cols].to_string(index=False))
        return

    # ---- Mode 2: global summary ----
    print("=" * 60)
    print("NORMALIZATION OUTPUT — SUMMARY")
    print("=" * 60)
    print(f"\nTotal proteins:    {len(proteins) if proteins is not None else '(no proteins.json)'}")
    print(f"Total annotations: {len(annotations)}")

    print("\n--- Annotations by tool ---")
    print(annotations.groupby("source_tool").size().to_string())

    print("\n--- Annotations by type ---")
    print(annotations.groupby("annotation_type").size().sort_values(ascending=False).to_string())

    print("\n--- Annotations by confidence level ---")
    print(annotations.groupby("confidence_level").size().to_string())

    print("\n--- Confidence × tool ---")
    print(
        annotations.groupby(["source_tool", "confidence_level"]).size()
        .unstack(fill_value=0).to_string()
    )

    # A few example proteins
    print("\n--- Sample proteins (first 5 with at least 5 annotations) ---")
    counts = annotations.groupby("uniprot_accession").size().sort_values(ascending=False)
    sample_acc = counts[counts >= 5].head(5).index.tolist()
    for acc in sample_acc:
        n = counts[acc]
        types = annotations[annotations["uniprot_accession"] == acc]["annotation_type"].unique()
        print(f"  {acc}: {n} annotations across {len(types)} types — try: --accession {acc}")

    print("\nTip: focus on a single protein with `--accession A0B9K2`")
    print(f"Tip: open {annot_path} directly in a text editor to inspect raw JSON")


if __name__ == "__main__":
    main()
