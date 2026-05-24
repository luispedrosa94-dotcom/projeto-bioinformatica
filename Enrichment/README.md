# Stage 2 — Enrichment

Enrichment stage of the protein annotation pipeline. Queries the UniProt
and InterPro APIs to complement the normalized output from Stage 1 with
curated protein metadata, GO aspect resolution, and protein signatures.

A second script (`consolidate.py`) merges the raw API responses with the
Stage 1 annotations into a single per-protein profile used by the
downstream stages and by the Streamlit explorer.

## Layout

```
.
├── configs/
│   └── default.yaml            # API config and scope
├── scripts/
│   ├── enrich.py               # Fetches UniProt + InterPro
│   └── consolidate.py          # Merges everything into protein_profiles.json
├── src/
│   ├── schema.py               # EnrichmentRecord (Pydantic)
│   ├── clients/
│   │   ├── uniprot.py          # UniProt REST client (parallel, with checkpoint)
│   │   └── interpro.py         # InterPro REST client
│   └── utils/
│       ├── checkpoint.py       # Resumable runs
│       ├── rate_limit.py       # Retry + exponential backoff
│       └── interpro_extract.py # InterPro field extraction
└── tests/
    └── test_clients.py         # Unit tests (mocked, offline)
```

## Setup

This project uses a conda environment. From the repository root:

```bash
# Create and activate the environment (one-time)
conda create -n stage3 python=3.11 -y
conda activate stage3

# Install Stage 2 dependencies
pip install -r Enrichment/requirements.txt
```

Dependencies pulled by `pip`:

- `requests` — HTTP client for UniProt and InterPro APIs
- `pydantic` — schema validation
- `pyyaml` — config file parsing
- `pytest` — test runner

The environment is shared with Stage 1 (Normalization) and Stage 3 (LLM
summarization). If you have already installed it for another stage,
activating it is enough.

## Usage

The stage is split in two steps that must be run in order.

### Step 1 — Fetch API data

```bash
python scripts/enrich.py --config configs/default.yaml
```

Options:
```bash
# Restrict to the 166 poorly annotated proteins (faster, for testing)
python scripts/enrich.py --config configs/default.yaml --scope poorly_annotated

# More parallel workers for the UniProt phase
python scripts/enrich.py --config configs/default.yaml --workers 20

# Skip the InterPro phase
python scripts/enrich.py --config configs/default.yaml --skip-interpro
```

What it does:
- Fetches the full UniProt JSON entry for each protein and caches it in
  `outputs/uniprot_raw/{accession}.json`.
- Fetches the InterPro coverage for each protein and caches it in
  `outputs/interpro_raw/{accession}.json`.
- Resolves the `GO_unknown` records produced in Stage 1 by UPIMAPI and
  eggNOG (which do not distinguish GO aspects) using the GO term → aspect
  map extracted from the UniProt records.
- Writes intermediate JSON files in `outputs/` (see below).
- Saves a checkpoint every 100 proteins, so an interrupted run can be
  resumed by re-running the same command.

### Step 2 — Consolidate

```bash
python scripts/consolidate.py --config configs/default.yaml
```

What it does:
- Reads the raw UniProt files from `outputs/uniprot_raw/` and extracts
  every available field.
- Merges them with the Stage 1 annotations and the InterPro section.
- Writes a single `outputs/protein_profiles.json` with one record per
  protein in a canonical schema (see below).

The consolidation step intentionally extracts every available UniProt
field, even those not currently used downstream. This avoids re-querying
the API when a new field is needed later.

## Tests

```bash
python -m pytest tests/ -v
```

## Outputs

All files are written to `../outputs/` (shared with Stage 1).

| File | Description |
|---|---|
| `uniprot_enrichment.json` | UniProt enrichment records (long format) |
| `interpro_enrichment.json` | InterPro entries (long format) |
| `go_aspect_map.json` | GO term → aspect mapping used for resolution |
| `annotations.json` | Stage 1 annotations with `GO_unknown` resolved |
| `protein_profiles.json` | Consolidated per-protein profile (one entry per protein) |
| `uniprot_raw/{acc}.json` | Cached raw UniProt JSON per protein |
| `interpro_raw/{acc}.json` | Cached raw InterPro JSON per protein |
| `checkpoints/` | Resume state for interrupted runs |

The `uniprot_raw/`, `interpro_raw/`, and `checkpoints/` directories are
not versioned (regenerable by re-running the pipeline). The consolidated
`protein_profiles.json` is also gitignored due to its size (~120 MB).

## Current numbers (full dataset, 1802 proteins)

| Item | Count |
|---|---|
| Proteins consolidated | 1802 |
| UniProt enrichment records | 19403 |
| InterPro enrichment records | 19919 |
| GO terms in aspect map | 900 |

UniProt API throughput: ~1 min for 1802 proteins with 20 parallel workers.
InterPro API throughput: ~6 min for 1802 proteins with 20 parallel workers.

## Consolidated schema

Each entry in `protein_profiles.json` has the following top-level keys:

- `accession`
- `identity` — protein name, gene name, organism, taxonomy, reviewed
  status, annotation score, entry version, public dates
- `function` — description, keywords, subcellular location, catalytic
  activity (with EC numbers, Rhea/CHEBI), pathway, subunit, similarity
- `go_annotations` — molecular function, biological process, cellular
  component (each with GO id, label, evidence, source, confidence)
- `enzymatic` — EC numbers with sources and confidence
- `domains` — COG, KOG, Pfam, TIGRFAM, SMART, COG categories
- `features` — domains, active sites, binding sites (with positions and
  ligands)
- `pathways` — KEGG KO, pathways, modules
- `sequence` — sequence value, length, molecular weight, CRC64
- `references` — literature references with PMID, title, authors, etc.
- `cross_references` — links to PDB, InterPro, Pfam, etc.
- `interpro` — InterPro signatures and integrated entries
- `evidence_summary` — tools used, annotation count, confidence
  distribution, overall confidence, poorly-annotated subset flag
- `provenance` — source, UniProt URL, raw file location

## Notes

- The InterPro client uses the `/entry/all/protein/UniProt/{accession}`
  endpoint to fetch every signature (Pfam, PANTHER, SUPERFAMILY, CATH,
  PROSITE, etc.) plus the integrated InterPro entries in a single
  request.
- `GO_unknown` resolution from UniProt avoids an additional API call:
  UniProt already returns GO terms with their aspect prefix (`F:`, `P:`,
  `C:`).
- Both clients implement exponential backoff with retry. Rate limits are
  handled in `rate_limit.py`.
