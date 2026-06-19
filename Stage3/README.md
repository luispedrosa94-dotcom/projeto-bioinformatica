# Stage 3 Protein Profile Summarization Toolkit

This folder implements a lightweight Stage 3 workflow:

```text
test_proteins.json
        ↓
one fresh Ollama request per protein
        ↓
structured protein information summary
        ↓
raw + parsed JSON outputs
        ↓
manual review CSV
```

## Goal

The goal is **not** to make the LLM choose a final biological annotation.

The goal is to test whether a local LLM can read one complete enriched protein profile and summarize all relevant information about that specific protein in a clear, cautious, reviewable way.

The summary should cover what is present in the record:

- identity and source status;
- reported protein name, gene name, organism, and metadata;
- curated or literature-backed function descriptions if present;
- GO terms and their evidence/source types;
- EC numbers, catalytic activities, reactions, Rhea/CHEBI references if present;
- domains, families, motifs, binding sites, active sites, and sequence features;
- pathways, KEGG, or contextual information;
- computational predictions from tools such as eggNOG, reCOGnizer, DeepGO2, DeepFRI, CLEAN, InterPro, or similar tools;
- strong/curated information;
- weak, predicted, indirect, sparse, missing, or conflicting information.

The output should help a supervisor or biologist inspect the protein faster. It should not replace biological curation.

## Contents

```text
Stage3/
  test_proteins.json
  build_test_set.py
  prompt_template.txt
  run_stage3.py
  make_review_sheet.py
  summarize_results.py
  check_test_set.py
  manual_review_rubric.csv
  requirements.txt
  schema/
    protein_profile_summary.schema.json
  docs/
    methodology.md
  outputs/
    raw_responses/
    prompts/
```


Expected design:

```text
25 proteins total
5 reviewed UniProt proteins
5 poorly annotated proteins
5 proteins with EC predictions
5 conflicting-evidence proteins
5 ML-only / weak-evidence proteins
```

These groups are still useful, but now they test whether the LLM can summarize different types of protein profiles, not whether it can recommend an annotation.

> **Note:** All examples below use `llama3.1` as the default model. To run with Qwen3.6-35B, replace `--model llama3.1` with `--model qwen3.6:35b` and use `--output-dir outputs/qwen3`. See the "Running with Qwen3.6-35B" section below for details.

## Dry run


This validates input handling and saves prompts without calling Ollama:

```bash
python run_stage3.py \
  --input test_proteins.json \
  --model llama3.1 \
  --dry-run \
  --save-prompts
```

## Run one protein first

Recommended before the full run:

```bash
python run_stage3.py \
  --input test_proteins.json \
  --model llama3.1 \
  --accession Q46505 \
  --output-dir outputs_test \
  --format-mode schema \
  --save-prompts
```

## Run the full Stage 3 summarization experiment

```bash
python run_stage3.py \
  --input test_proteins.json \
  --model llama3.1 \
  --output-dir outputs \
  --format-mode schema \
  --save-prompts
```

The runner saves after every protein:

```text
outputs/stage3_results.jsonl
outputs/stage3_results.current.json
outputs/raw_responses/<accession>.txt
outputs/prompts/<accession>.txt
```

## Run one group only

```bash
python run_stage3.py \
  --input test_proteins.json \
  --model llama3.1 \
  --group "ML-only / weak-evidence proteins" \
  --output-dir outputs_ml_only
```

## If JSON schema mode fails

Older Ollama versions or some models may behave better with plain JSON mode:

```bash
python run_stage3.py \
  --input test_proteins.json \
  --model llama3.1:8b \
  --format-mode json
```

If needed, disable format control entirely:

```bash
python run_stage3.py \
  --input test_proteins.json \
  --model llama3.1:8b \
  --format-mode none
```

## If the model context window is too small

Default behavior sends the full protein JSON. If a model fails because the prompt is too long, try:

```bash
python run_stage3.py \
  --input test_proteins.json \
  --model llama3.1:8b \
  --strip-empty
```

This removes only null values, empty strings, empty lists, and empty dictionaries. It does not intentionally remove biological evidence.


## Create the manual review sheet

```bash
python make_review_sheet.py \
  --results outputs/stage3_results.jsonl \
  --output outputs/review_sheet.csv
```

Open `outputs/review_sheet.csv` in Excel, Google Sheets, LibreOffice, or similar.

## Summarize results

```bash
python summarize_results.py --results outputs/stage3_results.jsonl
```

## Manual review rubric

Use `manual_review_rubric.csv`.

Suggested criteria:

| Criterion | Question |
|---|---|
| Completeness | Did the summary cover the main information categories present in the profile? |
| Faithfulness | Did the summary stay within the input JSON? |
| Clarity | Is the profile easier to understand than the raw JSON? |
| Evidence separation | Did it distinguish curated/reported evidence from predictions, weak evidence, and context? |
| Caution | Did it avoid becoming a final annotation or overclaiming? |

Suggested score:

| Score | Meaning |
|---:|---|
| 0 | Poor / wrong / unsupported |
| 1 | Weak / incomplete / vague |
| 2 | Useful but needs correction |
| 3 | Good, faithful, clear, and useful |

## Recommended interpretation

The outputs should be treated as **reviewable protein profile summaries**, not final biological annotations.

## Results (latest run)

The latest run (llama3.1 with `num_ctx=32768`, current prompt template, full
25-protein test set with InterPro v2 input) produced:

- **23/25** proteins completed successfully without warnings.
- **2/25** completed with one minor warning each (empty
  `reported_function_summary` on `H0UK06`, `W2LNG5` — both ML-only
  proteins where no curated function exists in the input record).
- **0/25** errors or timeouts.

### Prompt iteration — removing concrete examples

An earlier version of the prompt (rules 12–14) included six concrete example
IDs (`PF01257`, `EC 1.12.1.3`, `EC 1.6.5.3`, `GO:0051539`, `GO:0016491`,
`IPR036188`, `K18330`, `PubMed PMID`) intended to illustrate the format of
review notes and conflict descriptions. Systematic validation against the
input JSONs revealed that the model was copying these IDs into outputs of
proteins where the IDs did not appear in the input: **35 leaks across the
25 proteins**.

The most affected example was `PubMed PMID`, which leaked into nearly every
review note (20 of 25 proteins) regardless of whether the input contained
literature references at all. `PF01257` leaked into all 5 poorly annotated
proteins.

The prompt was revised to replace concrete examples with abstract
descriptions (e.g. *"a specific EC number reported in the input"* instead
of *"EC 1.12.1.3"*). Re-running with the revised prompt eliminated all
leaks — **35 → 0** — confirmed by the same cross-check script.

### Average bullets per summary section

| Section | Bullets |
|---|---:|
| `go_annotation_summary` | 4.6 |
| `enzyme_and_reaction_summary` | 1.8 |
| `domain_family_and_feature_summary` | 2.4 |
| `pathway_and_context_summary` | 1.8 |
| `tool_prediction_summary` | 2.2 |
| `strong_or_curated_information` | 1.7 |
| `weak_predicted_or_indirect_information` | 1.9 |
| `conflicting_or_inconsistent_information` | 0.2 |
| `missing_or_limited_information` | 1.8 |
| `review_notes` | 3.2 |

The drop in `conflicting_or_inconsistent_information` (0.9 → 0.2) is
explained by the leak removal: most "conflicts" reported in the previous
run referenced leaked example IDs and were therefore not real conflicts.
The remaining 0.2 are real conflicts present in the input.

Prompt token range: ~5,200 (smallest, `I9KF72` — poorly annotated) to
~28,700 (largest, `P06131` — reviewed UniProt with full InterPro coverage).
Maximum context usage: 88% of `num_ctx=32768` — no truncation.

---

## Running on tesla server (project-specific notes)

These notes apply to the bioinformatics pipeline server (`tesla.di.uminho.pt`).

### Ollama port

The default `--ollama-url` in `run_stage3.py` points to port `11435`. The Ollama
container in this project listens on port **`11434`**. Use the explicit flag:

```bash
--ollama-url http://localhost:11434/api/generate
```

### Conda environment

The required dependencies (`requests`) are already in the `stage3` conda
environment used by Stages 1 and 2. Activate it before running:

```bash
conda activate stage3
```

If the environment is missing dependencies, install them locally:

```bash
pip install -r requirements.txt
```

### Recommended first run (single protein)

Quick sanity check before doing the full 25:

```bash
cd /home/lpedrosa/projeto-bioinformatica/Stage3

python run_stage3.py \
  --input test_proteins.json \
  --model llama3.1 \
  --ollama-url http://localhost:11434/api/generate \
  --accession Q46505 \
  --output-dir outputs_test \
  --save-prompts
```

### Full 25-protein run

After the single-protein run looks correct:

```bash
python run_stage3.py \
  --input test_proteins.json \
  --model llama3.1 \
  --ollama-url http://localhost:11434/api/generate \
  --output-dir outputs \
  --save-prompts
```

The runner writes after each protein, so partial progress is preserved on
SSH disconnects, server reboots, or timeouts.


### Running with Qwen3.6-35B

Qwen3.6-35B is a sparse Mixture-of-Experts model. It was run in non-thinking
mode with Q4_K_M quantization:

```bash
python run_stage3.py \
  --input test_proteins.json \
  --model qwen3.6:35b \
  --ollama-url http://localhost:11434/api/generate \
  --output-dir outputs/qwen3 \
  --num-ctx 32768 \
  --timeout 600 \
  --save-prompts
```

The longer `--timeout 600` accounts for the slower model.

### Resuming after a partial / failed run

The runner skips successful accessions automatically. To also retry the
records that failed, add `--rerun-failed`:

```bash
python run_stage3.py --input test_proteins.json --model llama3.1 --rerun-failed
```

