"""
Stage 2b — Consolidation

Reads the raw UniProt JSON files (outputs/uniprot_raw/{acc}.json) saved by
enrich.py and extracts ALL available fields — nothing is filtered or lost.

This approach guarantees that if a new field is needed in the future, it can
be extracted without re-querying the UniProt API.

Schema per protein:
{
  "accession": "A0LIY7",
  "identity": {
    protein_name, alternative_names, gene_name, gene_name_synonyms,
    reviewed_status, annotation_score, protein_existence, entry_name,
    entry_version, first_public_date, last_annotation_update,
    organism { scientific_name, common_name, taxon_id, lineage }
  },
  "function": {
    description, keywords, subcellular_location,
    catalytic_activity [ { reaction, ec_number, rhea_ids, chebi_ids, evidences } ],
    pathway, subunit, similarity
  },
  "go_annotations": {
    molecular_function, biological_process, cellular_component
    — each: go_id, label, evidence_code, evidence_source, sources, confidence
  },
  "enzymatic": {
    ec_numbers [ { ec_id, sources, confidence } ]
  },
  "domains": { cog, kog, pfam, tigrfam, smart, cog_categories },
  "features": {
    domains, active_sites, binding_sites
    — each: description, position_start, position_end, evidences, ligand (binding)
  },
  "pathways": { kegg_ko, kegg_pathways, kegg_modules },
  "sequence": { value, length, mol_weight, crc64 },
  "references": [ { id, pmid, citation_type, title, authors, date, journal,
                    volume, pages, reference_positions, submission_database } ],
  "cross_references": { pdb, interpro, pfam, ... },
  "evidence_summary": { tools, annotation_count, confidence_distribution,
                        overall_confidence, in_poorly_annotated_subset },
  "provenance": { source, uniprot_url, db_source, raw_file }
}
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from collections import defaultdict
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def resolve_path(base: Path, p: str) -> Path:
    pp = Path(p)
    return pp if pp.is_absolute() else (base / pp).resolve()


def _dedupe_ordered(items: list) -> list:
    seen = set()
    out = []
    for x in items:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


def _best_conf(records: list[dict], annotation_type: str) -> str | None:
    order = {"high": 0, "medium": 1, "low": 2, "unknown": 3}
    matches = [r for r in records if r.get("annotation_type") == annotation_type]
    if not matches:
        return None
    matches.sort(key=lambda r: order.get(r.get("confidence_level", "unknown"), 3))
    return matches[0]["value"]


def _all_values_with_meta(records: list[dict], annotation_type: str) -> list[dict]:
    order = {"high": 0, "medium": 1, "low": 2, "unknown": 3}
    seen: dict[str, dict] = {}
    for r in records:
        if r.get("annotation_type") != annotation_type:
            continue
        val  = r["value"]
        tool = r.get("source_tool", "").replace("SourceTool.", "")
        conf = r.get("confidence_level", "unknown").replace("ConfidenceLevel.", "")
        if val not in seen:
            seen[val] = {"value": val, "label": r.get("label"), "sources": [], "confidence": conf}
        if tool and tool not in seen[val]["sources"]:
            seen[val]["sources"].append(tool)
        if order.get(conf, 3) < order.get(seen[val]["confidence"], 3):
            seen[val]["confidence"] = conf
        if not seen[val]["label"] and r.get("label"):
            seen[val]["label"] = r.get("label")
    return sorted(seen.values(), key=lambda x: order.get(x["confidence"], 3))


def _overall_confidence(records: list[dict]) -> str:
    counts = defaultdict(int)
    for r in records:
        c = r.get("confidence_level", "unknown").replace("ConfidenceLevel.", "")
        counts[c] += 1
    total = sum(counts.values())
    if total == 0:
        return "unknown"
    if counts["high"] / total >= 0.5:
        return "high"
    if (counts["high"] + counts["medium"]) / total >= 0.5:
        return "medium"
    return "low"


def _evidences(evs: list[dict]) -> list[dict]:
    return [{"code": e.get("evidenceCode", ""), "source": e.get("source", ""), "id": e.get("id", "")} for e in evs or []]


def _extract_uniprot(entry: dict, acc: str, ann: list[dict]) -> dict:
    """
    Extract ALL fields from a raw UniProt JSON entry.
    Reads directly from the raw response — nothing is filtered.
    """

    # ── Entry audit ────────────────────────────────────────────────────
    audit = entry.get("entryAudit", {})

    # ── Identity ───────────────────────────────────────────────────────
    entry_type = entry.get("entryType", "").lower()
    is_reviewed = "reviewed" in entry_type and "unreviewed" not in entry_type

    organism = entry.get("organism", {})

    # Protein names — recommended + ALL alternative + ALL submitted WITH evidences
    prot_desc = entry.get("proteinDescription", {})
    protein_name = None
    protein_name_evidences = []
    alternative_names = []

    rec_name = prot_desc.get("recommendedName", {})
    if rec_name:
        full = rec_name.get("fullName", {})
        protein_name = full.get("value", "") or None
        protein_name_evidences = _evidences(full.get("evidences", []))
        for alt in prot_desc.get("alternativeNames", []):
            alt_full = alt.get("fullName", {})
            alt_val = alt_full.get("value", "")
            if alt_val:
                alternative_names.append({
                    "name": alt_val,
                    "evidences": _evidences(alt_full.get("evidences", [])),
                    "short_names": [
                        {"value": s.get("value", ""), "evidences": _evidences(s.get("evidences", []))}
                        for s in alt.get("shortNames", [])
                    ],
                })

    # Submitted names (TrEMBL) — use last as primary (highest quality), others as alternatives
    submitted = prot_desc.get("submissionNames", [])
    if not protein_name and submitted:
        # Use last submitted name as primary (most recent submission)
        primary_sub = submitted[-1].get("fullName", {})
        protein_name = primary_sub.get("value", "") or None
        protein_name_evidences = _evidences(primary_sub.get("evidences", []))
        # All others go to alternative_names
        for sub in submitted[:-1]:
            sub_full = sub.get("fullName", {})
            sub_val = sub_full.get("value", "")
            if sub_val:
                alternative_names.append({
                    "name": sub_val,
                    "evidences": _evidences(sub_full.get("evidences", [])),
                    "short_names": [],
                })

    if not protein_name:
        protein_name = _best_conf(ann, "protein_description")

    # Gene names — with evidences for each ORF name
    gene_name_uniprot = None
    gene_synonyms = []
    gene_name_synonym_details = []
    for gene in entry.get("genes", []):
        if not gene_name_uniprot:
            gn = gene.get("geneName", {})
            gene_name_uniprot = gn.get("value", "") or None
        for s in gene.get("synonyms", []):
            gene_synonyms.append(s.get("value", ""))
        for o in gene.get("orfNames", []):
            orf_val = o.get("value", "")
            if orf_val:
                gene_synonyms.append(orf_val)
                gene_name_synonym_details.append({
                    "value": orf_val,
                    "type": "orf_name",
                    "evidences": _evidences(o.get("evidences", [])),
                })

    gene_name_inferred = _best_conf(ann, "preferred_gene_name") if not gene_name_uniprot else None

    identity = {
        "protein_name": protein_name,
        "protein_name_evidences": protein_name_evidences,
        "alternative_names": alternative_names,
        "gene_name": gene_name_uniprot,
        "gene_name_inferred": {
            "value": gene_name_inferred,
            "note": "inferred by homology from tools — not directly from UniProt for this protein",
        } if gene_name_inferred else None,
        "gene_name_synonyms": _dedupe_ordered(gene_synonyms),
        "gene_name_synonym_details": gene_name_synonym_details,
        "reviewed_status": "reviewed" if is_reviewed else "unreviewed",
        "entry_type": entry.get("entryType", "") or None,
        "annotation_score": entry.get("annotationScore"),
        "protein_existence": entry.get("proteinExistence", ""),
        "entry_name": entry.get("uniProtkbId", ""),
        "entry_version": audit.get("entryVersion"),
        "sequence_version": audit.get("sequenceVersion"),
        "first_public_date": audit.get("firstPublicDate"),
        "last_annotation_update": audit.get("lastAnnotationUpdateDate"),
        "last_sequence_update": audit.get("lastSequenceUpdateDate"),
        "organism": {
            "scientific_name": organism.get("scientificName", "") or None,
            "common_name": organism.get("commonName", "") or None,
            "taxon_id": str(organism.get("taxonId", "")) or None,
            "lineage": organism.get("lineage", []),
            "evidences": _evidences(organism.get("evidences", [])),
        },
        "uniprot_id": entry.get("extraAttributes", {}).get("uniParcId") or None,
    }

    # ── Function (from comments) ────────────────────────────────────────
    function_descriptions = []
    catalytic_activities = []
    pathways = []
    subunits = []
    similarities = []
    subcellular_locations = []

    # Keywords — with evidences and category
    keywords = []
    for kw in entry.get("keywords", []):
        kw_name = kw.get("name", "")
        if kw_name:
            keywords.append({
                "id": kw.get("id", "") or None,
                "name": kw_name,
                "category": kw.get("category", "") or None,
                "evidences": _evidences(kw.get("evidences", [])),
            })

    for comment in entry.get("comments", []):
        ctype = comment.get("commentType", "")

        if ctype == "FUNCTION":
            for text in comment.get("texts", []):
                val = text.get("value", "")
                if val:
                    function_descriptions.append({
                        "value": val,
                        "evidences": _evidences(text.get("evidences", [])),
                    })

        elif ctype == "CATALYTIC ACTIVITY":
            reaction = comment.get("reaction", {})
            reaction_name = reaction.get("name", "")
            ec_number = reaction.get("ecNumber", "")
            rhea_ids, chebi_ids = [], []
            for xref in reaction.get("reactionCrossReferences", []):
                db = xref.get("database", "")
                xid = xref.get("id", "")
                if db == "Rhea":
                    rhea_ids.append(xid)
                elif db == "ChEBI":
                    chebi_ids.append(xid)
            if reaction_name or ec_number:
                catalytic_activities.append({
                    "reaction": reaction_name or None,
                    "ec_number": ec_number or None,
                    "rhea_ids": rhea_ids,
                    "chebi_ids": chebi_ids,
                    "evidences": _evidences(reaction.get("evidences", [])),
                })

        elif ctype == "PATHWAY":
            for text in comment.get("texts", []):
                val = text.get("value", "")
                if val:
                    pathways.append({"value": val, "evidences": _evidences(text.get("evidences", []))})

        elif ctype == "SUBUNIT":
            for text in comment.get("texts", []):
                val = text.get("value", "")
                if val:
                    subunits.append({"value": val, "evidences": _evidences(text.get("evidences", []))})

        elif ctype == "SIMILARITY":
            for text in comment.get("texts", []):
                val = text.get("value", "")
                if val:
                    similarities.append({"value": val, "evidences": _evidences(text.get("evidences", []))})

        elif ctype == "SUBCELLULAR LOCATION":
            for loc_entry in comment.get("subcellularLocations", []):
                loc_obj = loc_entry.get("location", {})
                loc_val = loc_obj.get("value", "")
                if loc_val:
                    subcellular_locations.append({
                        "location": loc_val,
                        "location_id": loc_obj.get("id", "") or None,
                        "topology": loc_entry.get("topology", {}).get("value", "") or None,
                        "evidences": _evidences(loc_obj.get("evidences", [])),
                    })

    function = {
        "description": function_descriptions,
        "keywords": keywords,
        "subcellular_location": subcellular_locations,
        "catalytic_activity": catalytic_activities,
        "pathway": pathways,
        "subunit": subunits,
        "similarity": similarities,
    }

    # ── GO terms ────────────────────────────────────────────────────────
    aspect_map = {"F:": "molecular_function", "P:": "biological_process", "C:": "cellular_component"}
    go_from_uniprot: dict[str, list[dict]] = defaultdict(list)

    for xref in entry.get("uniProtKBCrossReferences", []):
        if xref.get("database") != "GO":
            continue
        go_id   = xref.get("id", "")
        props   = {p["key"]: p["value"] for p in xref.get("properties", [])}
        go_name = props.get("GoTerm", "")
        evidence = props.get("GoEvidenceType", "")
        prefix  = go_name[:2] if len(go_name) >= 2 else ""
        aspect  = aspect_map.get(prefix)
        if aspect and go_id:
            ev_parts = evidence.split(":", 1)
            go_from_uniprot[aspect].append({
                "go_id": go_id,
                "label": go_name[2:].strip(),
                "evidence_code": ev_parts[0] if ev_parts else None,
                "evidence_source": ev_parts[1] if len(ev_parts) > 1 else None,
            })

    order = {"high": 0, "medium": 1, "low": 2, "unknown": 3}

    # Source type per tool — what approach was used to identify the annotation
    _SOURCE_TYPE = {
        "upimapi":    "homology",
        "eggnog":     "orthology",
        "deepfri":    "ml_structure",
        "deepgo2":    "ml_sequence",
        "clean":      "ml_enzyme",
        "recognizer": "domain_search",
        "uniprot":    "curated",
    }

    def _build_go(aspect_ann_type: str, aspect_key: str) -> list[dict]:
        seen: dict[str, dict] = {}
        for r in ann:
            if r.get("annotation_type") != aspect_ann_type:
                continue
            go_id = r["value"]
            tool  = r.get("source_tool", "").replace("SourceTool.", "")
            conf  = r.get("confidence_level", "unknown").replace("ConfidenceLevel.", "")
            score = r.get("score")
            if go_id not in seen:
                seen[go_id] = {
                    "go_id": go_id,
                    "label": r.get("label"),
                    "evidence_code": None,
                    "evidence_source": None,
                    "sources": [],
                    "source_details": {},
                    "confidence": conf,
                }
            if tool and tool not in seen[go_id]["sources"]:
                seen[go_id]["sources"].append(tool)
                seen[go_id]["source_details"][tool] = {
                    "source_type": _SOURCE_TYPE.get(tool, "unknown"),
                    "confidence": conf,
                    "score": score,
                }
            if order.get(conf, 3) < order.get(seen[go_id]["confidence"], 3):
                seen[go_id]["confidence"] = conf

        for go in go_from_uniprot.get(aspect_key, []):
            go_id = go["go_id"]
            if go_id not in seen:
                seen[go_id] = {
                    "go_id": go_id,
                    "label": go["label"],
                    "evidence_code": go["evidence_code"],
                    "evidence_source": go["evidence_source"],
                    "sources": ["uniprot"],
                    "source_details": {
                        "uniprot": {
                            "source_type": "curated",
                            "evidence_code": go["evidence_code"],
                            "evidence_source": go["evidence_source"],
                        }
                    },
                    "confidence": "high",
                }
            else:
                if "uniprot" not in seen[go_id]["sources"]:
                    seen[go_id]["sources"].append("uniprot")
                    seen[go_id]["source_details"]["uniprot"] = {
                        "source_type": "curated",
                        "evidence_code": go["evidence_code"],
                        "evidence_source": go["evidence_source"],
                    }
                seen[go_id]["confidence"] = "high"
                if not seen[go_id]["evidence_code"]:
                    seen[go_id]["evidence_code"] = go["evidence_code"]
                if not seen[go_id]["evidence_source"]:
                    seen[go_id]["evidence_source"] = go["evidence_source"]
                if not seen[go_id]["label"] and go["label"]:
                    seen[go_id]["label"] = go["label"]

        return sorted(seen.values(), key=lambda x: order.get(x["confidence"], 3))

    go_annotations = {
        "molecular_function": _build_go("GO_MF", "molecular_function"),
        "biological_process": _build_go("GO_BP", "biological_process"),
        "cellular_component": _build_go("GO_CC", "cellular_component"),
    }

    # ── EC numbers ──────────────────────────────────────────────────────
    ec_items = _all_values_with_meta(ann, "EC")
    ec_seen = {item["value"] for item in ec_items}

    for act in catalytic_activities:
        ec = act.get("ec_number")
        if ec and ec not in ec_seen:
            ec_items.append({"value": ec, "label": None, "sources": ["uniprot"], "confidence": "high"})
            ec_seen.add(ec)
        elif ec:
            existing = next(x for x in ec_items if x["value"] == ec)
            if "uniprot" not in existing["sources"]:
                existing["sources"].append("uniprot")
                existing["confidence"] = "high"

    _SOURCE_TYPE_EC = {
        "upimapi":    "homology",
        "eggnog":     "orthology",
        "deepfri":    "ml_structure",
        "deepgo2":    "ml_sequence",
        "clean":      "ml_enzyme",
        "recognizer": "domain_search",
        "uniprot":    "curated",
    }

    enzymatic = {
        "ec_numbers": [
            {
                "ec_id": item["value"],
                "sources": item["sources"],
                "source_details": {
                    src: {"source_type": _SOURCE_TYPE_EC.get(src, "unknown")}
                    for src in item["sources"]
                },
                "confidence": item["confidence"],
            }
            for item in ec_items
        ]
    }

    # ── Features — ALL types ─────────────────────────────────────────────
    feature_domains, active_sites, binding_sites, other_features = [], [], [], []

    for feature in entry.get("features", []):
        ftype = feature.get("type", "")
        location = feature.get("location", {})
        feat = {
            "type": ftype,
            "description": feature.get("description", "") or None,
            "position_start": location.get("start", {}).get("value"),
            "position_end": location.get("end", {}).get("value"),
            "position_start_modifier": location.get("start", {}).get("modifier"),
            "position_end_modifier": location.get("end", {}).get("modifier"),
            "feature_id": feature.get("featureId", "") or None,
            "evidences": _evidences(feature.get("evidences", [])),
        }
        if ftype == "Domain":
            feature_domains.append(feat)
        elif ftype == "Active site":
            active_sites.append(feat)
        elif ftype == "Binding site":
            ligand = feature.get("ligand", {})
            feat["ligand_name"] = ligand.get("name", "") or None
            feat["ligand_id"]   = ligand.get("id", "") or None
            feat["ligand_note"] = ligand.get("note", "") or None
            feat["ligand_chebi"] = [
                x.get("id") for x in feature.get("featureCrossReferences", [])
                if x.get("database") == "ChEBI"
            ]
            binding_sites.append(feat)
        else:
            # Signal, Chain, Propeptide, Transit peptide, etc.
            other_features.append(feat)

    features = {
        "domains":        feature_domains,
        "active_sites":   active_sites,
        "binding_sites":  binding_sites,
        "other_features": other_features,
    }

    # ── Sequence ────────────────────────────────────────────────────────
    seq = entry.get("sequence", {})
    sequence = {
        "value":      seq.get("value", "") or None,
        "length":     seq.get("length"),
        "mol_weight": seq.get("molWeight"),
        "crc64":      seq.get("crc64", "") or None,
        "md5":        seq.get("md5", "") or None,
    } if seq else None

    # ── References — with referenceNumber, evidences and referenceComments ─
    references = []
    for ref in entry.get("references", []):
        citation = ref.get("citation", {})
        cite_type = citation.get("citationType", "")
        pmid = citation.get("id", "") if cite_type == "journal article" else None
        ref_id = pmid or citation.get("id", "")

        # Strain / tissue / plasmid info from referenceComments
        ref_comments = [
            {
                "type": rc.get("type", ""),
                "value": rc.get("value", ""),
                "evidences": _evidences(rc.get("evidences", [])),
            }
            for rc in ref.get("referenceComments", [])
        ]

        references.append({
            "reference_number": ref.get("referenceNumber"),
            "id": ref_id or None,
            "pmid": pmid or None,
            "citation_type": cite_type,
            "title": citation.get("title", "") or None,
            "authors": citation.get("authors", []),
            "authoring_group": citation.get("authoringGroup", []),
            "date": citation.get("publicationDate", "") or None,
            "journal": citation.get("journal", "") or None,
            "volume": citation.get("volume", "") or None,
            "pages": citation.get("firstPage", "") or None,
            "reference_positions": ref.get("referencePositions", []),
            "reference_comments": ref_comments,
            "submission_database": citation.get("submissionDatabase", "") or None,
            "evidences": _evidences(ref.get("evidences", [])),
        })

    # ── Cross-references — with all properties ───────────────────────────
    cross_refs: dict[str, list[dict]] = defaultdict(list)
    for xref in entry.get("uniProtKBCrossReferences", []):
        db = xref.get("database", "")
        xid = xref.get("id", "")
        if db and xid and db != "GO":
            props = {p["key"]: p["value"] for p in xref.get("properties", [])}
            cross_refs[db].append({
                "id": xid,
                "properties": props if props else None,
            })

    # ── Extra attributes ─────────────────────────────────────────────────
    extra = entry.get("extraAttributes", {})

    return {
        "identity": identity,
        "function": function,
        "go_annotations": go_annotations,
        "enzymatic": enzymatic,
        "features": features,
        "sequence": sequence,
        "references": references,
        "cross_references": dict(cross_refs),
        "extra_attributes": {
            "count_by_comment_type": extra.get("countByCommentType", {}),
            "count_by_feature_type": extra.get("countByFeatureType", {}),
        } if extra else None,
    }


def consolidate(
    proteins_path: Path,
    annotations_path: Path,
    raw_dir: Path,
    output_path: Path,
) -> None:
    log = logging.getLogger("consolidate")

    log.info("Loading annotation data...")
    with open(proteins_path) as f:
        proteins = json.load(f)
    with open(annotations_path) as f:
        annotations = json.load(f)

    log.info("Loaded %d proteins, %d annotations", len(proteins), len(annotations))

    ann_by_acc: dict[str, list[dict]] = defaultdict(list)
    for r in annotations:
        ann_by_acc[r["uniprot_accession"]].append(r)

    consolidated: list[dict] = []
    missing_raw = 0

    for prot in proteins:
        acc = prot["uniprot_accession"]
        ann = ann_by_acc.get(acc, [])

        # Read raw UniProt JSON
        raw_file = raw_dir / f"{acc}.json"
        uniprot_data = {}
        if raw_file.exists():
            try:
                with open(raw_file, encoding="utf-8") as f:
                    raw_entry = json.load(f)
                uniprot_data = _extract_uniprot(raw_entry, acc, ann)
            except Exception as e:
                log.warning("Failed to parse raw for %s: %s", acc, e)
                missing_raw += 1
        else:
            missing_raw += 1

        # Fallback identity from proteins.json
        identity = uniprot_data.get("identity", {
            "protein_name": _best_conf(ann, "protein_description"),
            "alternative_names": [],
            "gene_name": None,
            "gene_name_inferred": {
                "value": _best_conf(ann, "preferred_gene_name"),
                "note": "inferred by homology from tools — not directly from UniProt for this protein",
            },
            "gene_name_synonyms": [],
            "reviewed_status": prot.get("db_source", "unknown"),
            "annotation_score": None,
            "protein_existence": None,
            "entry_name": prot.get("entry_name"),
            "entry_version": None,
            "first_public_date": None,
            "last_annotation_update": None,
            "last_sequence_update": None,
            "organism": {"scientific_name": None, "common_name": None, "taxon_id": None, "lineage": []},
        })

        # Domains from Stage 1
        def _domain_list(atype):
            items = _all_values_with_meta(ann, atype)
            return [{"id": i["value"], "label": i.get("label"), "sources": i["sources"], "confidence": i["confidence"]} for i in items]

        domains = {
            "cog":     _domain_list("cog"),
            "kog":     _domain_list("kog"),
            "pfam":    _domain_list("pfam"),
            "tigrfam": _domain_list("tigrfam"),
            "smart":   _domain_list("smart"),
            "cog_categories": _dedupe_ordered([r["value"] for r in ann if r.get("annotation_type") == "cog_category"]),
        }

        # Pathways from Stage 1
        pathways = {
            "kegg_ko":       _dedupe_ordered([r["value"] for r in ann if r.get("annotation_type") == "KEGG_ko"]),
            "kegg_pathways": _dedupe_ordered([r["value"] for r in ann if r.get("annotation_type") == "KEGG_pathway"]),
            "kegg_modules":  _dedupe_ordered([r["value"] for r in ann if r.get("annotation_type") == "KEGG_module"]),
        }

        # Evidence summary
        tools = _dedupe_ordered([r["source_tool"].replace("SourceTool.", "") for r in ann])
        conf_counts: dict[str, int] = defaultdict(int)
        for r in ann:
            c = r.get("confidence_level", "unknown").replace("ConfidenceLevel.", "")
            conf_counts[c] += 1

        evidence_summary = {
            "tools": tools,
            "annotation_count": len(ann),
            "confidence_distribution": dict(conf_counts),
            "overall_confidence": _overall_confidence(ann),
            "in_poorly_annotated_subset": prot.get("in_poorly_annotated_subset", False),
        }

        provenance = {
            "source": "UniProtKB",
            "uniprot_url": f"https://www.uniprot.org/uniprotkb/{acc}/entry",
            "db_source": prot.get("db_source", "unknown"),
            "raw_file": f"{acc}.json",
        }

        consolidated.append({
            "accession": acc,
            "identity": identity,
            "function": uniprot_data.get("function", {
                "description": [], "keywords": [], "subcellular_location": [],
                "catalytic_activity": [], "pathway": [], "subunit": [], "similarity": [],
            }),
            "go_annotations": uniprot_data.get("go_annotations", {
                "molecular_function": [], "biological_process": [], "cellular_component": []
            }),
            "enzymatic": uniprot_data.get("enzymatic", {"ec_numbers": _all_values_with_meta(ann, "EC")}),
            "domains": domains,
            "features": uniprot_data.get("features", {"domains": [], "active_sites": [], "binding_sites": []}),
            "pathways": pathways,
            "sequence": uniprot_data.get("sequence"),
            "references": uniprot_data.get("references", []),
            "cross_references": uniprot_data.get("cross_references", {}),
            "extra_attributes": uniprot_data.get("extra_attributes"),
            "evidence_summary": evidence_summary,
            "provenance": provenance,
        })

    log.info("Consolidated %d proteins (%d without raw UniProt data)", len(consolidated), missing_raw)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(consolidated, f, indent=2, ensure_ascii=False)
    log.info("Written → %s", output_path)

    # Summary
    def _count(field_path):
        count = 0
        for r in consolidated:
            val = r
            for key in field_path:
                val = val.get(key, {}) if isinstance(val, dict) else {}
            if val:
                count += 1
        return count

    print("\n=== Consolidation summary ===")
    print(f"Proteins:                  {len(consolidated)}")
    print(f"With protein name:         {sum(1 for r in consolidated if r['identity']['protein_name'])}")
    print(f"With organism:             {sum(1 for r in consolidated if r['identity']['organism']['scientific_name'])}")
    print(f"With annotation score:     {sum(1 for r in consolidated if r['identity']['annotation_score'])}")
    print(f"With GO terms:             {sum(1 for r in consolidated if any(r['go_annotations'].values()))}")
    print(f"With EC number:            {sum(1 for r in consolidated if r['enzymatic']['ec_numbers'])}")
    print(f"With catalytic activity:   {sum(1 for r in consolidated if r['function']['catalytic_activity'])}")
    print(f"With function description: {sum(1 for r in consolidated if r['function']['description'])}")
    print(f"With features:             {sum(1 for r in consolidated if any(r['features'].values()))}")
    print(f"With references:           {sum(1 for r in consolidated if r['references'])}")
    print(f"With sequence:             {sum(1 for r in consolidated if r['sequence'])}")
    print(f"Without raw UniProt data:  {missing_raw}")

    conf_dist: dict[str, int] = defaultdict(int)
    for r in consolidated:
        conf_dist[r["evidence_summary"]["overall_confidence"]] += 1
    print(f"\nOverall confidence:")
    for c in ["high", "medium", "low", "unknown"]:
        print(f"  {c}: {conf_dist[c]}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/default.yaml")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    cfg_path = Path(args.config)
    with open(cfg_path) as f:
        cfg = yaml.safe_load(f)
    base = cfg_path.parent
    output_root = resolve_path(base, cfg["output_root"])

    consolidate(
        proteins_path    = resolve_path(base, cfg["proteins_path"]),
        annotations_path = resolve_path(base, cfg["annotations_path"]),
        raw_dir          = output_root / "uniprot_raw",
        output_path      = output_root / "protein_profiles.json",
    )


if __name__ == "__main__":
    main()
