# Integrative Enrichment of Protein Functional Annotations Using Large Language Models

Bioinformatics pipeline for the integration, normalization, and enrichment of
protein functional annotation results from multiple computational tools.
 
---

## Pipeline Overview

```
Raw tool outputs
      │
      ▼
┌─────────────────┐
│  Stage 1        │  Normalization
│  Normalization  │  Reads 6 tool outputs → canonical long-format schema
└────────┬────────┘  234 031 records, 1 802 proteins
         │
         ▼
┌─────────────────┐
│  Stage 2        │  Enrichment
│  Enrichment     │  UniProt API + STRING API → curated metadata
└────────┬────────┘  14 353 UniProt records, 171 STRING terms
         │
         ▼
┌─────────────────┐
│  Stage 3        │  [Planned] LLM Harmonization
│  Harmonization  │  Claude API → protein-level functional summaries
└─────────────────┘
```

## Tools in scope

| Tool | Type | Output |
|---|---|---|
| UPIMAPI | Homology (DIAMOND vs UniProt) | GO, EC, Pfam, descriptions |
| eggNOG-mapper | Orthology | GO, KEGG, COG, Pfam, EC |
| reCOGnizer | Domain search (CDD) | COG, KOG, Pfam, TIGR, SMART, EC, KEGG |
| DeepFRI | Structure + ML | GO_MF |
| DeepGO2 | ML | GO_BP, GO_MF, GO_CC |
| CLEAN | ML | EC numbers |

**Out of scope:** Foldseek, ColabFold (per supervisor decision).

## Repository structure

```
projeto-bioinformatica/
├── README.md                  ← this file
├── data/
│   └── raw_outputs/           ← tool output files (input to Stage 1)
├── outputs/                   ← pipeline outputs (written by Stages 1 & 2)
│   ├── annotations.json       ← 234 031 normalized + enriched records
│   ├── proteins.json          ← 1 802 unique proteins
│   ├── uniprot_enrichment.json
│   ├── string_enrichment.json
│   └── go_aspect_map.json
├── Normalization/             ← Stage 1 (see Normalization/README.md)
└── Enrichment/                ← Stage 2 (see Enrichment/README.md)
```

## Quick start

```bash
# Stage 1 — Normalization
cd Normalization
pip install -r requirements.txt
python -m pytest tests/ -v          # 60 tests
python scripts/normalize.py --config configs/default.yaml

# Stage 2 — Enrichment
cd ../Enrichment
pip install -r requirements.txt
python -m pytest tests/ -v          # 21 tests
python scripts/enrich.py --config configs/default.yaml
```

## Key design decisions

- **Long format** — one record per annotation, not one row per protein. Makes
  merging and downstream processing trivial.
- **Categorical confidence level** — every annotation carries `high/medium/low`
  alongside the raw score, giving a uniform signal across incompatible scales
  (e-value vs ML confidence).
- **JSON output** — human-readable, no extra dependencies (vs parquet).
- **GO aspect deferred** — UPIMAPI and eggNOG mix BP/MF/CC; aspect is resolved
  in Stage 2 using UniProt data.
- **Foldseek/ColabFold excluded** — per supervisor decision; schema enums
  retain their entries for future re-inclusion without migration.
>>>>>>> master
