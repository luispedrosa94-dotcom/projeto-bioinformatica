# Methodology: Stage 3 Protein Profile Summarization

## Objective

Stage 3 evaluates whether a local LLM can summarize complete enriched protein profiles into readable, evidence-aware summaries for manual inspection.

The objective is not to make the model assign or recommend a final biological annotation. Instead, the model receives all available information for one protein and summarizes what the profile contains.

## Input

The input is `test_proteins.json`, a representative set of 25 proteins. Each protein is processed independently. The full protein JSON object is inserted into a fixed prompt.

The protein profiles may contain identity metadata, UniProt information, function descriptions, GO annotations, EC numbers, catalytic activities, domain/family annotations, sequence features, pathway information, computational predictions, source details, confidence categories, and provenance.

## LLM task

For each protein, the LLM produces a structured JSON summary with sections for:

- identity and source status;
- reported function information;
- GO annotations;
- enzyme and reaction information;
- domain, family, and feature information;
- pathway, STRING, KEGG, or contextual information;
- tool predictions;
- strong or curated information;
- weak, predicted, or indirect information;
- conflicts or inconsistencies;
- missing or limited information;
- human review notes.

## Independence between proteins

Each protein is sent as a fresh Ollama request. Previous protein outputs are not passed into the next request. This avoids cross-protein contamination.

## Output handling

For every protein, the runner saves:

- raw model response;
- parsed JSON response when possible;
- minimal structural validation warnings;
- the input protein unless `--no-keep-input` is used;
- selected Ollama timing/count metadata.

The runner saves after every protein using JSONL so partial progress is preserved if a later request fails.

## Validation

Validation is intentionally minimal. The script checks whether:

- the response is parseable JSON;
- required fields are present;
- `protein_id` matches the input accession;
- list fields are lists;
- string fields are non-empty strings.

The script does not validate biological correctness. Manual review is required.

## Manual review

The review sheet supports simple scoring for completeness, faithfulness, clarity, evidence separation, and caution. The goal is to assess whether the summary helps a human understand the protein profile faster and more safely than reading raw JSON.

## Limitations

The LLM may omit details, misread evidence strength, over-summarize large records, or phrase predicted information too strongly. Therefore, outputs should be treated as review aids, not as final biological claims.
