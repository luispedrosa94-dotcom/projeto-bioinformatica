# Raw tool outputs

This folder contains the raw output of six annotation tools run on the
1802-protein dataset, plus a curated subset file. It is the input to
Stage 1 (Normalization).

These files were produced outside this repository by a previous step of
the project. They are versioned here so the pipeline is fully reproducible
from this single repo.

**Do not modify the files in this folder.** Stage 1 parses them as
read-only inputs.

## Subfolders

### upimapi/
Output of UPIMAPI (DIAMOND-based homology search against UniProt).

- `UPIMAPI_results.tsv` — main results table (annotations per protein)
- `uniprotinfo.tsv` — UniProt entry information for the hits
- `aligned.blast`, `unaligned.blast` — raw BLAST alignment files
- `valid_ids.txt`, `not_valid_ids.txt` — accession validation lists

### eggnogmapper_results/
Output of eggNOG-mapper (orthology-based annotation).

- `eggnog_mapper_results.emapper.annotations` — final annotations
- `eggnog_mapper_results.emapper.hits` — HMMER hits
- `eggnog_mapper_results.emapper.seed_orthologs` — seed orthologs

### recognizer_results/
Output of reCOGnizer (CDD-based domain search).

- `reCOGnizer_results.tsv` and `.xlsx` — combined results
- `COG_report.tsv`, `KOG_report.tsv`, `Pfam_report.tsv`, `TIGR_report.tsv`, `SMART_report.tsv`, `PRK_report.tsv`, `NCBI_Curated_report.tsv` — per-database reports
- `COG_quantification.tsv`, `KOG_quantification.tsv` — quantification summaries

### deepfri/
Output of DeepFRI (structure-based ML predictions, molecular function).

- `DeepFRI_MF_predictions.csv` — final predictions
- `DeepFRI_MF_pred_scores.json` — full per-GO score matrix (~43 MB)

### deepgo2/
Output of DeepGO2 (ML predictions across all three GO aspects).

- `subset_f_preds_bp.tsv` — biological process predictions
- `subset_f_preds_mf.tsv` — molecular function predictions
- `subset_f_preds_cc.tsv` — cellular component predictions

### clean/
Output of CLEAN (ML predictions for EC numbers).

- `subset_f_maxsep.csv` — EC predictions per protein

## The proteinas_poorly_annotated_sorted file

A long-format TSV containing 484 protein-tool rows covering 166 unique
proteins identified as poorly annotated in UniProt. This is a curated
subset of the 1802 proteins, intended to highlight where the pipeline
adds most value.

Five of the 25 proteins in the Stage 3 test set (the "Poorly annotated
proteins" group) are drawn from this file: B3E8W6, I9KF72, G8TMW4,
Q01Q37, A1ATR7.

## Total size

About 65 MB across 29 versioned files. The DeepFRI scores file
(`DeepFRI_MF_pred_scores.json`, 43 MB) is the largest.