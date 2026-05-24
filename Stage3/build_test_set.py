"""
Builds test_proteins.json from a curated list of 25 accessions
selected by the supervisor across 5 representative groups.

Usage:
    cd Stage3
    python build_test_set.py

The script auto-discovers protein_profiles.json in (in order):
  1. The current directory (Stage3/)
  2. ../outputs/ (the default location for the consolidated pipeline output)
  3. outputs/ (if launched from the project root)

The output test_proteins.json is written to the current directory.
"""
import json
from pathlib import Path

# ── Supervisor's curated selection ────────────────────────────────────────────

SELECTION = [
    # Group 1 — Reviewed UniProt proteins
    {"accession": "Q46505",  "group": "Reviewed UniProt proteins"},
    {"accession": "Q2LXU2",  "group": "Reviewed UniProt proteins"},
    {"accession": "Q2RI41",  "group": "Reviewed UniProt proteins"},
    {"accession": "P11560",  "group": "Reviewed UniProt proteins"},
    {"accession": "P06131",  "group": "Reviewed UniProt proteins"},

    # Group 2 — Poorly annotated proteins
    {"accession": "B3E8W6",  "group": "Poorly annotated proteins"},
    {"accession": "I9KF72",  "group": "Poorly annotated proteins"},
    {"accession": "G8TMW4",  "group": "Poorly annotated proteins"},
    {"accession": "Q01Q37",  "group": "Poorly annotated proteins"},
    {"accession": "A1ATR7",  "group": "Poorly annotated proteins"},

    # Group 3 — Proteins with EC predictions
    {"accession": "U6EEC5",  "group": "Proteins with EC predictions"},
    {"accession": "G3F6H2",  "group": "Proteins with EC predictions"},
    {"accession": "H1YXH1",  "group": "Proteins with EC predictions"},
    {"accession": "A8UKK6",  "group": "Proteins with EC predictions"},
    {"accession": "H2EIF6",  "group": "Proteins with EC predictions"},

    # Group 4 — Conflicting evidence proteins
    {"accession": "H1XYF8",  "group": "Conflicting evidence proteins"},
    {"accession": "F6D5H3",  "group": "Conflicting evidence proteins"},
    {"accession": "U6EF42",  "group": "Conflicting evidence proteins"},
    {"accession": "K2R6I7",  "group": "Conflicting evidence proteins"},
    {"accession": "F0T7S9",  "group": "Conflicting evidence proteins"},

    # Group 5 — ML-only / weak-evidence proteins
    {"accession": "U6EF83",  "group": "ML-only / weak-evidence proteins"},
    {"accession": "H0UK06",  "group": "ML-only / weak-evidence proteins"},
    {"accession": "W2LNG5",  "group": "ML-only / weak-evidence proteins"},
    {"accession": "L0DDH1",  "group": "ML-only / weak-evidence proteins"},
    {"accession": "J0HCT8",  "group": "ML-only / weak-evidence proteins"},
]

# ── Load profiles ─────────────────────────────────────────────────────────────

candidates = [
    "protein_profiles.json",
    "../outputs/protein_profiles.json",
    "outputs/protein_profiles.json",
]

profiles_path = None
for c in candidates:
    if Path(c).exists():
        profiles_path = c
        break

if profiles_path is None:
    print("ERROR: protein_profiles.json not found.")
    raise SystemExit(1)

print(f"Loading profiles from {profiles_path}...")
with open(profiles_path, encoding="utf-8") as f:
    all_profiles = json.load(f)
print(f"Loaded {len(all_profiles)} proteins")

# ── Build index ───────────────────────────────────────────────────────────────

index = {p["accession"]: p for p in all_profiles}

# ── Extract selected proteins ─────────────────────────────────────────────────

selected = []
missing  = []

for item in SELECTION:
    acc   = item["accession"]
    group = item["group"]
    if acc in index:
        p = dict(index[acc])
        p["_test_group"] = group
        selected.append(p)
    else:
        missing.append(acc)
        print(f"  WARNING: {acc} not found in protein_profiles.json")

# ── Save output ───────────────────────────────────────────────────────────────

out_path = Path("test_proteins.json")
with open(out_path, "w", encoding="utf-8") as f:
    json.dump(selected, f, indent=2, ensure_ascii=False)

print(f"\nSaved {len(selected)} proteins → {out_path}")

if missing:
    print(f"Missing accessions ({len(missing)}): {', '.join(missing)}")

# ── Print summary ─────────────────────────────────────────────────────────────

print("\n=== Selection summary ===")
groups: dict[str, list] = {}
for p in selected:
    g = p["_test_group"]
    groups.setdefault(g, []).append(p)

for group, prots in groups.items():
    print(f"\n{group} ({len(prots)})")
    for p in prots:
        acc   = p["accession"]
        name  = (p.get("identity") or {}).get("protein_name") or "Uncharacterized protein"
        conf  = (p.get("evidence_summary") or {}).get("overall_confidence", "—")
        rev   = (p.get("identity") or {}).get("reviewed_status", "—")
        print(f"  {acc}  [{rev}] [{conf}]  {name[:60]}")
