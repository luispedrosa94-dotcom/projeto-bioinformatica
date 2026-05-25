# App — Streamlit Explorer

Interactive dashboard for exploring the consolidated protein profile
dataset produced by Stages 1 and 2 of the pipeline. Shows per-protein
identity, GO annotations, enzymatic data, domains, InterPro entries,
sequence info, references, and the Stage 3 LLM summary when available.

## Setup

Activate the shared conda environment and install the app dependencies:

```bash
conda activate stage3
pip install streamlit pandas plotly
```

## Run

Launch from the **repository root** (not from inside `App/`):

```bash
cd /path/to/projeto-bioinformatica
streamlit run App/app.py --server.port 8501
```

The app auto-discovers `protein_profiles.json` in `outputs/` and the
Stage 3 results in `Stage3/outputs/stage3_results.jsonl` (when present).

Open the URL printed in the terminal (usually `http://localhost:8501`)
to access the dashboard.

## What it shows

Four top-level tabs:

- **Overview** — summary metrics, distributions, and aggregated charts
  across the whole dataset (1802 proteins).
- **Table** — filterable table of proteins with all enrichment fields.
- **Protein detail** — per-protein view selected by accession, with
  nine sub-tabs (Summary, GO & enzyme, Domains & pathways, InterPro,
  Sequence, References & xrefs, Origin/meta, Raw JSON) and an inline
  Stage 3 LLM summary block above them for the 25 test-set proteins.
- **Export** — download filtered subsets as CSV or JSON.

## Inputs

| Path | Required | Notes |
|---|---|---|
| `../outputs/03_consolidated/protein_profiles.json` | yes | Consolidated profile from Stage 2 |
| `../Stage3/outputs/stage3_results.jsonl` | optional | LLM summaries; if absent, the Stage 3 block is hidden |

If `protein_profiles.json` is missing, the sidebar shows a file uploader
as a fallback.

## Layout

​```
App/
└── app.py    # single-file Streamlit application (~1400 lines)
​```

The app is a single Python module; no separate package structure or
build step is needed.