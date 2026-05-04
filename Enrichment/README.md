# Stage 2 — Enrichment

Enrichment stage of the protein annotation pipeline. Queries the **UniProt**
and **STRING** APIs to complement the normalized outputs from Stage 1 with
curated protein metadata, GO aspect resolution, and protein interaction data.

## Layout

```
.
├── configs/
│   └── default.yaml          # API config (batch sizes, delays, scope)
├── scripts/
│   └── enrich.py             # Orchestrator entry point
├── src/
│   ├── schema.py             # EnrichmentRecord (Pydantic)
│   ├── clients/
│   │   ├── uniprot.py        # UniProt REST API client (parallel, with checkpoint)
│   │   └── string_db.py      # STRING API client
│   └── utils/
│       ├── checkpoint.py     # Save/load progress for resumable runs
│       └── rate_limit.py     # Retry + exponential backoff helpers
└── tests/
    └── test_clients.py       # 21 unit tests (mocked, offline)
```

## Usage

```bash
pip install -r requirements.txt
python -m pytest tests/ -v
python scripts/enrich.py --config configs/default.yaml
```

Options:
```bash
# Run only on poorly annotated proteins (faster, for testing)
python scripts/enrich.py --config configs/default.yaml --scope poorly_annotated

# Skip STRING queries
python scripts/enrich.py --config configs/default.yaml --skip-string

# More parallel workers for UniProt (faster, ~30s vs ~1min)
python scripts/enrich.py --config configs/default.yaml --workers 20
```

## What each API provides

### UniProt REST API (`/uniprotkb/{accession}`)

Queried individually per protein using 10 parallel threads (checkpoint saved
every 100 proteins for resumable runs). Provides:

- **Reviewed status** — SwissProt (manually curated) vs TrEMBL (automatic)
- **Protein name** — canonical recommended name (or submitted name for TrEMBL)
- **Gene name** — primary gene symbol
- **GO terms** — with aspect already assigned (GO_BP / GO_MF / GO_CC) via the
  `F:` / `P:` / `C:` prefixes in the UniProt response
- **EC numbers** — curated enzymatic classification
- **Subcellular location** — cell compartment information
- **Function description** — free-text functional summary
- **Keywords** — controlled vocabulary tags (e.g. ATP synthesis, Transport)

### GO aspect resolution (from UniProt data)

The `GO_unknown` records produced in Stage 1 by UPIMAPI and eggNOG (which do
not distinguish GO aspects) are resolved using the GO term → aspect map
extracted from the UniProt records fetched above. No additional API call is
needed — UniProt already returns GO terms with their aspect.

### STRING API

Provides protein interaction context for the full protein set:

- **Identifier mapping** — maps UniProt accessions to STRING IDs
- **Functional enrichment** — GO terms, KEGG pathways overrepresented in the
  protein set vs genome background (set-level analysis)

> Note: protein-protein interaction edges are not expected for metagenomic
> datasets where proteins come from different organisms.

## Outputs (written to `../outputs/`)

| File | Description |
|---|---|
| `annotations.json` | Stage 1 annotations updated — GO_unknown terms resolved to GO_BP/MF/CC |
| `uniprot_enrichment.json` | EnrichmentRecords from UniProt (14 353 records) |
| `string_enrichment.json` | EnrichmentRecords from STRING (171 enrichment terms) |
| `go_aspect_map.json` | GO term → aspect mapping used for resolution (audit trail) |
| `checkpoints/uniprot.json` | UniProt fetch progress (can be deleted after completion) |

## Pipeline results (full dataset, 1 802 proteins)

| Step | Result |
|---|---|
| UniProt records | 14 353 across 9 types |
| GO_unknown resolved | 9 457 / 29 866 total |
| GO_unknown still unresolved | 20 409 (TrEMBL proteins without GO annotations in UniProt) |
| STRING proteins mapped | 589 / 1 802 |
| STRING enrichment terms | 171 |

UniProt records by type:

| Type | Count |
|---|---|
| keyword | 6 556 |
| GO_MF | 2 214 |
| reviewed_status | 1 802 |
| GO_BP | 1 028 |
| GO_CC | 784 |
| protein_name | 709 |
| gene_name | 465 |
| subcellular_location | 462 |
| function_description | 333 |

## Design decisions

1. **Individual UniProt queries** (`/uniprotkb/{acc}`) instead of batch
   endpoints — batch endpoints consistently returned 400 errors due to URL
   length limits with the required fields. Individual queries with 10 parallel
   workers achieve comparable throughput (~1 min for 1 802 proteins).
2. **GO aspect resolution from UniProt data** — UniProt already returns GO
   terms with aspect prefixes (`F:`, `P:`, `C:`), so no additional API is
   needed to resolve `GO_unknown` terms. This keeps the pipeline limited to
   the two APIs specified in the project proposal (UniProt + STRING).
3. **STRING set-level enrichment** — the enrichment endpoint is more
   statistically meaningful for metagenomic datasets than per-protein queries,
   since proteins from different organisms cannot form interaction networks.
4. **Checkpoint system** — UniProt progress is saved every 100 proteins.
   If the pipeline is interrupted, re-running the same command resumes from
   where it stopped automatically.
