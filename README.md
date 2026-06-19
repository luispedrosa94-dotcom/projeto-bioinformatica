# Integrative Enrichment of Protein Functional Annotations Using Large Language Models

Bioinformatics pipeline for the integration, normalization, and enrichment
of protein functional annotation results from multiple computational tools.
The pipeline processes 1802 proteins through three sequential stages and
provides an interactive web dashboard for inspecting the results.

**Author:** Luís Pedrosa
**Institution:** Centre of Biological Engineering, University of Minho

## Pipeline

The pipeline has three stages plus an interactive explorer:

- **Stage 1 — Normalization** (`Normalization/`): reads the raw output of
  six annotation tools (UPIMAPI, eggNOG-mapper, reCOGnizer, DeepFRI,
  DeepGO2, CLEAN) and converts them into a canonical long-format schema
  with categorical confidence levels.
- **Stage 2 — Enrichment** (`Enrichment/`): queries the UniProt and
  InterPro APIs to complement the normalized records with curated
  metadata, GO aspect resolution, and protein signatures. A consolidation
  step then merges everything into a single per-protein profile
  (`protein_profiles.json`).
- **Stage 3 — LLM summarization** (`Stage3/`): runs a local Ollama model
  over a curated 25-protein test set, producing a structured,
  schema-validated summary per protein. The LLM acts as a reviewer's
  assistant, not as a decision-maker.
- **App — Streamlit explorer** (`App/`): interactive dashboard over the
  consolidated dataset, including the Stage 3 LLM summaries when
  available.

Each stage has its own `README.md` with details, configuration, and
commands.

## Repository structure

```
projeto-bioinformatica/
├── README.md                  this file
├── App/                       Streamlit explorer
├── Normalization/             Stage 1 code, tests, configs
├── Enrichment/                Stage 2 code, tests, configs
├── Stage3/                    Stage 3 toolkit, prompt, test set, outputs
├── data/
│   └── raw_outputs/           input from the six annotation tools
├── outputs/                   pipeline outputs (Stages 1 + 2)
│   ├── 01_normalization/      Stage 1 outputs (annotations, proteins)
│   ├── 02_enrichment/         Stage 2 API outputs (UniProt, InterPro, GO aspect map)
│   ├── 03_consolidated/       Stage 2 consolidated profile (protein_profiles.json)
│   └── caches/                raw API responses and checkpoints (gitignored)
├── article/                   intermediate report PDF
└── .gitignore
```

The `outputs/` folder holds the shared intermediate files (Stages 1 and
2 write here, the App and Stage 3 read from here). Large generated
artefacts such as `protein_profiles.json` and the raw API caches are not
versioned (they are regenerable by re-running the pipeline).

## Setup

The project uses a single conda environment shared by all stages.

```bash
conda create -n stage3 python=3.11 -y
conda activate stage3

pip install -r Normalization/requirements.txt
pip install -r Enrichment/requirements.txt
pip install -r Stage3/requirements.txt
pip install streamlit pandas plotly
```

## Quickstart

To reproduce the full pipeline from the raw tool outputs to the Streamlit
dashboard, run each stage in order from the repository root.

```bash
# Stage 1 — Normalization
python Normalization/scripts/normalize.py --config Normalization/configs/default.yaml

# Stage 2 — Enrichment (UniProt + InterPro)
python Enrichment/scripts/enrich.py --config Enrichment/configs/default.yaml
python Enrichment/scripts/consolidate.py --config Enrichment/configs/default.yaml

# Stage 3 — LLM summarization (requires Ollama running locally)
cd Stage3
python run_stage3.py \
  --input test_proteins.json \
  --model llama3.1 \
  --num-ctx 32768 \
  --save-prompts

# Run with Qwen3.6-35B (non-thinking mode)
python run_stage3.py \
  --input test_proteins.json \
  --model qwen3.6:35b \
  --num-ctx 32768 \
  --output-dir outputs/qwen3 \
  --save-prompts
cd ..

# App — Streamlit dashboard
streamlit run App/app.py --server.port 8501
```

See each stage's `README.md` for available options and configuration.

## Tools used in Stage 1

| Tool | Type | Outputs extracted |
|---|---|---|
| UPIMAPI | Homology (DIAMOND vs UniProt) | GO, EC, Pfam, descriptions |
| eggNOG-mapper | Orthology | GO, KEGG, COG, Pfam, EC |
| reCOGnizer | Domain search (CDD) | COG, KOG, Pfam, TIGR, SMART, EC, KEGG |
| DeepFRI | Structure + ML | GO_MF |
| DeepGO2 | ML | GO_BP, GO_MF, GO_CC |
| CLEAN | ML | EC numbers |

Foldseek and ColabFold are out of scope in this iteration; the schema
enums retain entries for future re-inclusion without migration.

## Current dataset

| Item | Count |
|---|---:|
| Proteins | 1,802 |
| Stage 1 annotation records | 234,031 |
| UniProt enrichment records (Stage 2) | 19,403 |
| InterPro enrichment records (Stage 2) | 19,919 |
| GO terms in aspect map | 900 |
| Stage 3 test set proteins | 25 |
| Stage 3 success rate (latest run) | 23/25 |

## Design decisions

- **Long format throughout** — one record per annotation, not one row per
  protein. Merging and downstream processing become append operations.
- **Categorical confidence level** — every annotation carries
  `high / medium / low / unknown` alongside the raw score, giving a
  uniform ordinal signal across incompatible scales (e-value vs ML
  confidence).
- **JSON output** — human-readable, no extra dependencies.
- **GO aspect resolution deferred to Stage 2** — UPIMAPI and eggNOG mix
  BP/MF/CC; the aspect is resolved using the GO term → aspect map
  extracted from UniProt records (no extra API call needed).
- **LLM as reviewer's assistant** — Stage 3 does not pick a single
  annotation. It produces a structured summary that a curator can scan
  quickly.
