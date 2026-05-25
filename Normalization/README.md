# Stage 1 — Normalization

Normalization layer for the protein annotation enrichment pipeline. Reads
heterogeneous tool outputs and produces two JSON files in a canonical
long-format schema, ready for Stage 2 (Enrichment) and Stage 3 (LLM summarization).

## Where this fits in the pipeline

Stage 1 is the entry point. It takes the raw output of six annotation
tools (see `data/raw_outputs/`) and produces two JSON files used by
every downstream stage:

- `outputs/01_normalization/annotations.json` — long-format record of every annotation
- `outputs/01_normalization/proteins.json` — one record per unique protein

Stage 2 (`Enrichment/`) reads these files and adds UniProt + InterPro
data. Stage 3 (`Stage3/`) consumes the consolidated profiles produced
by Stage 2. See the root README for the full pipeline overview.

## Layout

```
.
├── configs/
│   └── default.yaml         # Pipeline config (paths, scope, thresholds)
├── scripts/
│   ├── normalize.py         # Orchestrator entry point
│   └── peek.py              # Exploratory data analysis helper
├── src/
│   ├── schema.py            # Pydantic schema (AnnotationRecord, ProteinRecord)
│   ├── config.py            # YAML loader + scope filter + thresholds
│   ├── parsers/
│   │   ├── base.py          # BaseParser interface
│   │   ├── upimapi.py       # UPIMAPI parser
│   │   ├── eggnog.py        # eggNOG-mapper parser
│   │   ├── recognizer.py    # reCOGnizer parser
│   │   ├── deepfri.py       # DeepFRI parser
│   │   ├── deepgo2.py       # DeepGO2 parser
│   │   └── clean.py         # CLEAN parser
│   └── utils/
│       ├── ids.py           # UniProt ID parsing (handles all variants)
│       └── confidence.py    # Categorical confidence mapping rules
└── tests/
    └── test_parsers.py      # 60 unit tests
```

## Setup

This project uses a conda environment shared across all stages. From the
repository root:

```bash
conda activate stage3
pip install -r Normalization/requirements.txt
```

Dependencies: `pandas`, `pydantic`, `pyyaml`, `pytest`.

## Usage

```bash
python -m pytest tests/ -v
python scripts/normalize.py --config configs/default.yaml
```

Outputs written to `../outputs/01_normalization/`:
- `annotations.json` — one record per annotation (long format)
- `proteins.json` — one record per unique protein

## Scope decisions (from supervision meeting)

- **Tools in scope:** UPIMAPI, eggNOG-mapper, reCOGnizer, DeepFRI, DeepGO2,
  CLEAN — six homology- and ML-based tools.
- **Tools out of scope (this iteration):** Foldseek, ColabFold (structural
  evidence). Their support remains in the schema enums for future
  re-inclusion without migration.
- **Confidence model:** every annotation carries a categorical
  `confidence_level` (high / medium / low / unknown), in addition to the raw
  numeric `score`. This implements the supervisor's principle of "not
  treating all annotation evidence as equivalent".

## Canonical schema

**`annotations.json`** — long format, one row per piece of evidence:

| field | description |
|---|---|
| `uniprot_accession` | canonical key (e.g. `A0B9K2`) |
| `original_id` | ID as it appeared in the source file |
| `source_tool` | which tool produced this row |
| `annotation_type` | enum: `GO_BP`, `GO_MF`, `GO_CC`, `GO_unknown`, `EC`, `pfam`, `protein_description`, ... |
| `value` | the annotation itself |
| `label` | human-readable label (e.g. GO term name) |
| `score` | numeric confidence/significance (raw, kept for auditing) |
| `score_type` | enum: `evalue`, `bitscore`, `pident`, `confidence`, `plddt`, `none` |
| **`confidence_level`** | enum: `high`, `medium`, `low`, `unknown` — the abstraction layer |
| `evidence_rank` | preserves order when a tool emits multiple ranked hits |
| `raw_extras` | tool-specific metadata not fitting the canonical fields |

**`proteins.json`** — one row per unique accession.

## Confidence-level mapping (starting thresholds)

These are starting values from the supervision meeting. The YAML config
exposes them under `confidence_thresholds.<tool>` for empirical refinement.

| Tool | Score type | HIGH | MEDIUM | LOW |
|---|---|---|---|---|
| UPIMAPI / eggNOG / reCOGnizer | e-value | ≤ 1e-50 | 1e-50 to 1e-10 | > 1e-10 |
| DeepFRI / DeepGO2 | ML score (0..1) | ≥ 0.7 | 0.3 to 0.7 | < 0.3 |
| CLEAN | ML score (0..1) | ≥ 0.5 | 0.1 to 0.5 | < 0.1 |

**Empirical observation on starting thresholds:** with the full dataset of
1 802 proteins, UPIMAPI yields 98.8% of annotations in HIGH — too coarse
to discriminate. Recommended action: review the empirical e-value distribution
per tool with the supervisor and tighten thresholds before the final pipeline
run. The YAML config supports this without code changes.

## Current status

| Tool | Parser | Records | Status |
|---|---|---|---|
| UPIMAPI | ✓ | 10 561 | Extracts protein description, GO, EC, Pfam, function CC, family |
| eggNOG-mapper | ✓ | 54 833 | Extracts GO, KEGG, COG, Pfam, EC, eggNOG OGs |
| reCOGnizer | ✓ | 26 662 | Extracts COG, KOG, Pfam, TIGR, SMART, PRK, EC, KEGG KO |
| DeepFRI | ✓ | 8 660 | Extracts GO_MF predictions with ML confidence |
| DeepGO2 | ✓ | 131 072 | Extracts GO_BP, GO_MF, GO_CC predictions with ML confidence |
| CLEAN | ✓ | 2 243 | Extracts EC number predictions with ML confidence |
| **Total** | | **234 031** | |

## Full pipeline output (scope=all, default thresholds)

- **1 802** unique proteins
- **234 031** annotation records across 25 annotation types

Confidence distribution by tool:

| Tool | HIGH | MEDIUM | LOW |
|---|---|---|---|
| UPIMAPI | 10 439 | 120 | 2 |
| eggNOG | 50 957 | 3 306 | 570 |
| reCOGnizer | 11 038 | 10 026 | 5 598 |
| DeepFRI | 2 426 | 954 | 5 280 |
| DeepGO2 | 12 036 | 34 413 | 84 623 |
| CLEAN | 501 | 137 | 1 605 |

## Design decisions for the final article

1. **Long format** instead of wide — each annotation is independent; merges
   and enrichment become append operations.
2. **JSON output** instead of parquet — no extra dependencies (`pyarrow`),
   human-readable, and compatible with all Python versions including 3.13.
3. **Categorical `confidence_level`** alongside raw `score` — gives downstream
   stages (notably the LLM harmonizer) a uniform ordinal signal across
   incompatible score scales (e-value vs ML confidence) while preserving
   the raw value for auditing.
4. **GO aspect deferred:** UPIMAPI and eggNOG mix BP/MF/CC, so records are
   tagged `GO_unknown` and the aspect is resolved in Stage 2 via the UniProt API.
5. **Tool-specific confidence thresholds in YAML:** the heterogeneity of
   score scales is acknowledged at the config layer, not hidden in code.
6. **Foldseek / ColabFold excluded:** per supervisor decision; schema enums
   retain their entries for future re-inclusion without migration.
