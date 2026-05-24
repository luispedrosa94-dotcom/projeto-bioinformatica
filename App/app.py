"""
Protein Profiles Explorer — Streamlit dashboard

Reads protein_profiles.json (from ../outputs/) and provides interactive
exploration of the protein dataset.

Run from the project root:
    cd /path/to/projeto-bioinformatica
    streamlit run App/app.py --server.port 8501

The app expects to be launched from the project root so it can find
protein_profiles.json in ../outputs/ relative to its own location.

"""

import json
import io
import zipfile
import math
from pathlib import Path

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

# ── Default data file lookup paths ────────────────────────────────────────────
# Resolved relative to this script, not to the current working directory,
# so the app can be launched from anywhere (e.g. streamlit run App/app.py).
SCRIPT_DIR = Path(__file__).parent
DEFAULT_DATA_PATHS = [
    SCRIPT_DIR / "protein_profiles.json",
    SCRIPT_DIR / "protein_profiles.zip",
    SCRIPT_DIR.parent / "outputs" / "protein_profiles.json",
    SCRIPT_DIR.parent / "outputs" / "protein_profiles.zip",
]

# ── Stage 3 LLM results lookup paths ──────────────────────────────────────────
DEFAULT_STAGE3_PATHS = [
    SCRIPT_DIR / "stage3_results.jsonl",
    SCRIPT_DIR.parent / "Stage3" / "outputs" / "stage3_results.jsonl",
]


@st.cache_data(show_spinner=False)
def load_stage3_results() -> dict[str, dict]:
    """Try to load Stage 3 LLM results from known paths.

    Returns a dict accession -> record. Empty dict if no file found.
    """
    for candidate in DEFAULT_STAGE3_PATHS:
        if candidate.exists():
            index: dict[str, dict] = {}
            try:
                with open(candidate, encoding="utf-8") as f:
                    for line in f:
                        line = line.strip()
                        if not line:
                            continue
                        rec = json.loads(line)
                        acc = rec.get("accession")
                        if acc:
                            index[acc] = rec
                return index
            except Exception:
                return {}
    return {}

# ── Page config ───────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="Protein Profiles Explorer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── Custom CSS ────────────────────────────────────────────────────────────────
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;600&family=IBM+Plex+Sans:wght@300;400;500;600&display=swap');
    
    html, body, [class*="css"] { font-family: 'IBM Plex Sans', sans-serif; }
    
    .metric-card {
        background: #f8f9fa;
        border: 1px solid #e9ecef;
        border-radius: 10px;
        padding: 22px 24px;
        text-align: center;
        box-shadow: 0 1px 3px rgba(0,0,0,0.04);
        transition: box-shadow 0.2s, transform 0.2s;
    }
    .metric-card:hover {
        box-shadow: 0 4px 12px rgba(0,0,0,0.08);
        transform: translateY(-1px);
    }
    .metric-card .label { font-size: 13px; color: #6c757d; text-transform: uppercase; letter-spacing: 0.06em; margin-bottom: 6px; font-weight: 500; }
    .metric-card .value { font-size: 34px; font-weight: 700; color: #111827; font-family: 'IBM Plex Mono', monospace; line-height: 1.2; }
    .metric-card .sub   { font-size: 12px; color: #28a745; margin-top: 4px; font-weight: 500; }
    
    .hero {
        background: linear-gradient(135deg, #f0fdf4 0%, #e0f2fe 100%);
        border: 1px solid #d1fae5;
        border-radius: 12px;
        padding: 28px 32px;
        margin-bottom: 24px;
    }
    .hero h1 { font-size: 32px; font-weight: 700; margin: 0; color: #064e3b; }
    .hero p  { color: #065f46; margin: 6px 0 0; font-size: 14px; }
    
    .section-title { font-size: 18px; font-weight: 600; color: #1f2937; margin: 24px 0 12px; border-bottom: 2px solid #e5e7eb; padding-bottom: 6px; }
    
    .badge {
        display: inline-block;
        padding: 2px 10px;
        border-radius: 12px;
        font-size: 12px;
        font-weight: 500;
        margin: 2px;
        background: #e0f2fe;
        color: #0369a1;
    }
    .badge-green  { background: #d1fae5; color: #065f46; }
    .badge-yellow { background: #fef9c3; color: #854d0e; }
    .badge-red    { background: #fee2e2; color: #991b1b; }
    .badge-gray   { background: #f3f4f6; color: #374151; }
    
    .evidence-row {
        background: #f9fafb;
        border-left: 3px solid #6366f1;
        padding: 8px 12px;
        margin: 4px 0;
        border-radius: 0 4px 4px 0;
        font-size: 13px;
    }
    .evidence-row .eco { font-family: 'IBM Plex Mono', monospace; font-size: 11px; color: #6b7280; }
    
    .tip-box {
        background: #eff6ff;
        border: 1px solid #bfdbfe;
        border-radius: 6px;
        padding: 10px 14px;
        font-size: 13px;
        color: #1e40af;
        margin-bottom: 16px;
    }
    
    div[data-testid="stSidebarNav"] { display: none; }
</style>
""", unsafe_allow_html=True)

# ── Data loading ──────────────────────────────────────────────────────────────

@st.cache_data(show_spinner="Loading proteins…")
def load_data(raw_bytes: bytes, filename: str) -> list[dict]:
    if filename.endswith(".zip"):
        with zipfile.ZipFile(io.BytesIO(raw_bytes)) as zf:
            json_names = [n for n in zf.namelist() if n.endswith(".json")]
            if not json_names:
                st.error("No JSON file found inside ZIP.")
                return []
            with zf.open(json_names[0]) as f:
                return json.load(f)
    else:
        return json.loads(raw_bytes.decode("utf-8"))


@st.cache_data(show_spinner="Building table…")
def flatten_proteins(proteins: list[dict]) -> pd.DataFrame:
    rows = []
    for p in proteins:
        ident = p.get("identity", {})
        org   = ident.get("organism", {})
        ev    = p.get("evidence_summary", {})
        seq   = p.get("sequence") or {}
        go    = p.get("go_annotations", {})
        n_go  = sum(len(v) for v in go.values())
        conf_dist = ev.get("confidence_distribution", {})
        rows.append({
            "accession":          p.get("accession", ""),
            "protein_name":       ident.get("protein_name") or "—",
            "gene_name":          ident.get("gene_name") or (
                                      (ident.get("gene_name_inferred") or {}).get("value") or None),
            "organism":           org.get("scientific_name") or "Unknown",
            "taxon_id":           org.get("taxon_id"),
            "reviewed_status":    ident.get("reviewed_status", "unreviewed"),
            "annotation_score":   ident.get("annotation_score"),
            "protein_existence":  ident.get("protein_existence", ""),
            "overall_confidence": ev.get("overall_confidence", "unknown"),
            "annotation_count":   ev.get("annotation_count", 0),
            "n_go":               n_go,
            "n_ec":               len((p.get("enzymatic") or {}).get("ec_numbers", [])),
            "n_tools":            len(ev.get("tools", [])),
            "poorly_annotated":   ev.get("in_poorly_annotated_subset", False),
            "seq_length":         seq.get("length"),
            "mol_weight":         seq.get("mol_weight"),
            "conf_high":          conf_dist.get("high", 0),
            "conf_medium":        conf_dist.get("medium", 0),
            "conf_low":           conf_dist.get("low", 0),
        })
    return pd.DataFrame(rows)


# ── Sidebar ───────────────────────────────────────────────────────────────────

with st.sidebar:
    st.markdown("### 🧬 Protein report")
    st.markdown("**Upload protein JSON or ZIP**")
    uploaded = st.file_uploader("", type=["json", "zip"], label_visibility="collapsed")

    # Try bundled file
    proteins_raw = None
    bundle_name  = None
    if uploaded:
        proteins_raw = uploaded.read()
        bundle_name  = uploaded.name
    else:
        for candidate in DEFAULT_DATA_PATHS:
            if candidate.exists():
                bundle_name = candidate.name
                with open(candidate, "rb") as f:
                    proteins_raw = f.read()
                st.success(f"Loaded bundled file:\n{candidate.name}")
                break
    if proteins_raw is None:
        st.info("Upload a protein_profiles.json or .zip file to begin.")
        st.stop()

    all_proteins = load_data(proteins_raw, bundle_name)
    df_all = flatten_proteins(all_proteins)

    st.markdown("---")
    st.markdown("### Filters")

    search_q = st.text_input("Search", placeholder="accession, gene, organism, GO/domain…")

    rev_opts    = ["reviewed", "unreviewed"]
    rev_default = rev_opts
    rev_sel     = st.multiselect("Reviewed status", rev_opts, default=rev_default)

    conf_opts    = ["high", "medium", "low", "unknown"]
    conf_default = conf_opts
    conf_sel     = st.multiselect("Overall confidence", conf_opts, default=conf_default)

    tool_opts = sorted({t for tools in df_all["n_tools"] for t in []})
    # Collect all unique tools
    all_tools = sorted({t for p in all_proteins for t in p.get("evidence_summary", {}).get("tools", [])})
    tool_sel  = st.multiselect("Evidence tools/sources", all_tools, default=[])

    all_orgs = sorted(df_all["organism"].dropna().unique())
    org_sel  = st.multiselect("Organism", all_orgs, default=[])

    seq_min = int(df_all["seq_length"].dropna().min()) if df_all["seq_length"].notna().any() else 0
    seq_max = int(df_all["seq_length"].dropna().max()) if df_all["seq_length"].notna().any() else 10000
    seq_range = st.slider("Sequence length", seq_min, seq_max, (seq_min, seq_max))

    score_min = float(df_all["annotation_score"].dropna().min()) if df_all["annotation_score"].notna().any() else 0.0
    score_max = float(df_all["annotation_score"].dropna().max()) if df_all["annotation_score"].notna().any() else 5.0
    score_range = st.slider("Annotation score", score_min, score_max, (score_min, score_max))

    go_min = st.number_input("Minimum GO annotations", min_value=0, value=0, step=1)


# ── Filtering ─────────────────────────────────────────────────────────────────

def apply_filters(df: pd.DataFrame, proteins: list[dict]) -> tuple[pd.DataFrame, list[dict]]:
    mask = pd.Series([True] * len(df), index=df.index)

    if search_q:
        q = search_q.lower()
        mask &= (
            df["accession"].str.lower().str.contains(q, na=False) |
            df["protein_name"].str.lower().str.contains(q, na=False) |
            df["organism"].str.lower().str.contains(q, na=False) |
            df["gene_name"].astype(str).str.lower().str.contains(q, na=False)
        )

    if rev_sel:
        mask &= df["reviewed_status"].isin(rev_sel)

    if conf_sel:
        mask &= df["overall_confidence"].isin(conf_sel)

    if org_sel:
        mask &= df["organism"].isin(org_sel)

    if tool_sel:
        def has_tools(p):
            tools = p.get("evidence_summary", {}).get("tools", [])
            return all(t in tools for t in tool_sel)
        acc_with_tools = {p["accession"] for p in proteins if has_tools(p)}
        mask &= df["accession"].isin(acc_with_tools)

    mask &= (df["seq_length"].isna() | ((df["seq_length"] >= seq_range[0]) & (df["seq_length"] <= seq_range[1])))
    mask &= (df["annotation_score"].isna() | ((df["annotation_score"] >= score_range[0]) & (df["annotation_score"] <= score_range[1])))
    mask &= df["n_go"] >= go_min

    filtered_df  = df[mask].copy()
    filtered_acc = set(filtered_df["accession"])
    filtered_prot = [p for p in proteins if p["accession"] in filtered_acc]
    return filtered_df, filtered_prot


df_filtered, proteins_filtered = apply_filters(df_all, all_proteins)

# ── Helper functions ──────────────────────────────────────────────────────────

def conf_badge(c: str) -> str:
    color = {"high": "green", "medium": "yellow", "low": "red"}.get(c, "gray")
    return f'<span class="badge badge-{color}">{c}</span>'


def tool_badge(t: str) -> str:
    return f'<span class="badge">{t}</span>'

def _interpro_render_section(interpro: dict) -> None:
    """Render the complete InterPro section for one protein."""
    entries = interpro.get("entries", []) or []
    summary = interpro.get("summary", {}) or {}

    if not entries:
        st.info("This protein has no InterPro matches.")
        return

    # ── Summary metrics ──────────────────────────────────────────────────
    total       = summary.get("total_entries", 0)
    n_ipr       = summary.get("total_interpro_integrated", 0)
    n_unint     = summary.get("total_unintegrated", 0)
    by_src      = summary.get("by_source_database", {}) or {}
    by_type     = summary.get("by_type", {}) or {}
    member_dbs  = summary.get("member_databases_used", []) or []
    go_counts   = summary.get("go_terms_count", {}) or {}
    go_terms    = summary.get("go_terms_list", []) or []

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Total entries", total)
    m2.metric("IPR canónicos", n_ipr)
    m3.metric("Member-db signatures", total - n_ipr)
    m4.metric("Unintegrated", n_unint)

    # Badges with source DBs found
    if by_src:
        st.markdown("**Source databases found:**")
        st.markdown(
            " ".join(f'<span class="badge">{src} ({n})</span>'
                     for src, n in sorted(by_src.items(), key=lambda x: -x[1])),
            unsafe_allow_html=True,
        )

    # Badges with types found
    if by_type:
        st.markdown("**Entry types:**")
        st.markdown(
            " ".join(f'<span class="badge badge-gray">{t} ({n})</span>'
                     for t, n in sorted(by_type.items(), key=lambda x: -x[1])),
            unsafe_allow_html=True,
        )

    # Member databases used inside IPR entries
    if member_dbs:
        st.markdown("**Member databases used (inside IPR entries):**")
        st.markdown(
            " ".join(f'<span class="badge">{db}</span>' for db in member_dbs),
            unsafe_allow_html=True,
        )

    # Aggregated GO terms across all entries
    if any(go_counts.values()):
        st.markdown(
            f"**GO terms aggregated:** "
            f'<span class="badge badge-green">BP: {go_counts.get("biological_process", 0)}</span> '
            f'<span class="badge badge-green">MF: {go_counts.get("molecular_function", 0)}</span> '
            f'<span class="badge badge-green">CC: {go_counts.get("cellular_component", 0)}</span>',
            unsafe_allow_html=True,
        )

    st.markdown("---")

    # ── Filters ──────────────────────────────────────────────────────────
    f1, f2 = st.columns(2)
    with f1:
        view_filter = st.radio(
            "Show:",
            ["All", "Only IPR canonical", "Only signatures (non-IPR)"],
            horizontal=True,
            key=f"ipro_view_{id(interpro)}",
        )
    with f2:
        sources_available = sorted({e.get("source_database") or "unknown" for e in entries})
        source_filter = st.selectbox(
            "Source database:",
            options=["All"] + sources_available,
            key=f"ipro_src_{id(interpro)}",
        )

    # Apply filters
    filtered = entries
    if view_filter == "Only IPR canonical":
        filtered = [e for e in filtered if e.get("source_database") == "interpro"]
    elif view_filter == "Only signatures (non-IPR)":
        filtered = [e for e in filtered if e.get("source_database") != "interpro"]
    if source_filter != "All":
        filtered = [e for e in filtered if e.get("source_database") == source_filter]

    st.caption(f"Showing {len(filtered)} of {total} entries")

    if not filtered:
        st.info("No entries match the current filters.")
        return

    # ── Compact table view ───────────────────────────────────────────────
    table_rows = []
    for e in filtered:
        # Get the first protein location's positions (usually only one)
        positions = "—"
        score = "—"
        prot_blocks = e.get("proteins", []) or []
        if prot_blocks:
            locs = prot_blocks[0].get("entry_protein_locations", []) or []
            if locs:
                first_loc = locs[0]
                fragments = first_loc.get("fragments", []) or []
                if fragments:
                    parts = []
                    for fr in fragments:
                        parts.append(f"{fr.get('start','?')}-{fr.get('end','?')}")
                    positions = ", ".join(parts)
                sc = first_loc.get("score")
                if sc is not None and sc != 0:
                    score = f"{sc:.2e}" if isinstance(sc, (int, float)) and sc != 0 else str(sc)
                elif sc == 0:
                    score = "0"

        # GO terms count for this entry
        go_dict = e.get("go_terms", {}) or {}
        n_go = sum(len(v or []) for v in go_dict.values())

        table_rows.append({
            "Accession":      e.get("accession") or "—",
            "Name":           e.get("name") or "—",
            "Source DB":      e.get("source_database") or "—",
            "Type":           e.get("type") or "—",
            "Integrated to":  e.get("integrated") or "—",
            "Positions":      positions,
            "Score":          score,
            "GO terms":       n_go,
        })

    df_table = pd.DataFrame(table_rows)
    st.dataframe(df_table, use_container_width=True, hide_index=True)

    # ── Detailed cards ───────────────────────────────────────────────────
    st.markdown("### Detailed view")
    st.caption("Click each entry to expand.")

    for e in filtered:
        acc = e.get("accession") or "?"
        name = e.get("name") or "(no name)"
        src = e.get("source_database") or "?"
        etype = e.get("type") or "?"
        is_ipr = (src == "interpro")
        badge_class = "badge-green" if is_ipr else "badge"

        header = (
            f"**`{acc}`** — {name}"
            f" · {src}/{etype}"
        )
        with st.expander(header):
            # Identity row
            cols = st.columns(3)
            cols[0].markdown(f"**Source DB**\n\n{src}")
            cols[1].markdown(f"**Type**\n\n{etype}")
            integrated = e.get("integrated")
            cols[2].markdown(
                f"**Integrated to**\n\n"
                + (f"`{integrated}`" if integrated else "—")
            )

            # Description
            desc = e.get("description")
            if desc:
                if isinstance(desc, list):
                    desc_text = "\n\n".join(desc)
                else:
                    desc_text = str(desc)
                st.markdown(f"**Description**\n\n{desc_text}")

            # Member databases (only for IPR entries)
            member_dbs_dict = e.get("member_databases") or {}
            if member_dbs_dict:
                st.markdown("**Member databases inside this IPR:**")
                for db_name, sigs in member_dbs_dict.items():
                    if isinstance(sigs, list) and sigs:
                        sigs_str = ", ".join(
                            f"`{s.get('accession','?')}` ({s.get('name','?')})"
                            for s in sigs
                        )
                        st.markdown(f"- **{db_name}**: {sigs_str}")

            # GO terms
            go_dict = e.get("go_terms") or {}
            go_lines = []
            for cat, terms in go_dict.items():
                if terms:
                    for term in terms:
                        gid = term.get("go_id")
                        gname = term.get("name") or "(no name)"
                        if gid:
                            go_lines.append(f"- `{gid}` ({cat}) — {gname}")
            if go_lines:
                st.markdown("**GO terms:**")
                st.markdown("\n".join(go_lines))

            # Per-protein locations
            prot_blocks = e.get("proteins") or []
            if prot_blocks:
                pb = prot_blocks[0]
                length = pb.get("protein_length")
                organism = pb.get("organism")
                in_af = pb.get("in_alphafold")
                in_bf = pb.get("in_bfvd")
                meta_line_parts = []
                if length is not None:
                    meta_line_parts.append(f"protein length: {length} aa")
                if organism:
                    meta_line_parts.append(f"organism: {organism}")
                if in_af is not None:
                    meta_line_parts.append(
                        f"AlphaFold: {'yes' if in_af else 'no'}"
                    )
                if in_bf is not None:
                    meta_line_parts.append(
                        f"BFVD: {'yes' if in_bf else 'no'}"
                    )
                if meta_line_parts:
                    st.caption(" · ".join(meta_line_parts))

                locs = pb.get("entry_protein_locations") or []
                if locs:
                    st.markdown("**Protein locations:**")
                    loc_rows = []
                    for li, loc in enumerate(locs, 1):
                        fragments = loc.get("fragments", []) or []
                        frag_strs = []
                        for fr in fragments:
                            s = fr.get("start", "?")
                            ed = fr.get("end", "?")
                            dc = fr.get("dc-status", "")
                            frag_strs.append(f"{s}-{ed}" + (f" [{dc}]" if dc and dc != "CONTINUOUS" else ""))
                        loc_rows.append({
                            "#":            li,
                            "Fragments":    "; ".join(frag_strs),
                            "Representative": "yes" if loc.get("representative") else "no",
                            "Model":        loc.get("model") or "—",
                            "Score":        f"{loc.get('score'):.2e}" if isinstance(loc.get("score"), (int, float)) and loc.get("score") not in (0, None) else str(loc.get("score") if loc.get("score") is not None else "—"),
                        })
                        # Subfamily (PANTHER)
                        sub = loc.get("subfamily")
                        if sub:
                            sub_acc = sub.get("accession", "?")
                            sub_name = sub.get("name", "?")
                            st.caption(f"Subfamily: `{sub_acc}` — {sub_name}")
                    st.dataframe(pd.DataFrame(loc_rows), use_container_width=True, hide_index=True)

            # Hierarchy
            hier = e.get("hierarchy") or {}
            if hier:
                with st.expander("🌳 Hierarchy", expanded=False):
                    st.json(hier, expanded=False)

            # Cross-references
            xrefs = e.get("cross_references") or {}
            if xrefs:
                with st.expander("🔗 Cross-references", expanded=False):
                    st.json(xrefs, expanded=False)

            # Literature
            lit = e.get("literature") or {}
            if lit:
                with st.expander("📚 Literature", expanded=False):
                    st.json(lit, expanded=False)




# ── Stage 3 LLM Summary render ────────────────────────────────────────────────
def _llm_summary_render_section(record: dict) -> None:
    """Render the Stage 3 LLM summary for a given protein record."""
    if not record:
        return

    parsed = record.get("llm_response_json") or {}
    status = record.get("status", "")
    model = record.get("model", "—")
    group = record.get("_test_group") or record.get("test_group", "")
    warnings = record.get("validation_warnings") or []
    pt = record.get("ollama_prompt_eval_count")
    rt = record.get("ollama_eval_count")
    dur = record.get("ollama_total_duration_seconds")

    st.markdown("### 🤖 Stage 3 LLM Summary")

    meta_bits = []
    if status:
        meta_bits.append(f"**Status:** {status}")
    if model:
        meta_bits.append(f"**Model:** {model}")
    if group:
        meta_bits.append(f"**Group:** {group}")
    if pt:
        meta_bits.append(f"**Prompt:** {pt:,} tokens")
    if rt:
        meta_bits.append(f"**Response:** {rt:,} tokens")
    if dur:
        meta_bits.append(f"**Time:** {dur:.1f}s")
    if meta_bits:
        st.markdown(" · ".join(meta_bits))

    if warnings:
        st.warning("⚠ " + " ".join(warnings))

    # Overall profile
    overall = parsed.get("overall_profile_summary", "").strip()
    if overall:
        st.markdown("**Overall profile**")
        st.markdown(overall)

    # Identity + reported function side by side
    identity = parsed.get("identity_summary", "").strip()
    function = parsed.get("reported_function_summary", "").strip()
    if identity or function:
        c1, c2 = st.columns(2)
        with c1:
            if identity:
                st.markdown("**Identity**")
                st.markdown(identity)
        with c2:
            if function:
                st.markdown("**Reported function**")
                st.markdown(function)

    # Bullet sections
    list_sections = [
        ("go_annotation_summary",                   "GO annotations"),
        ("enzyme_and_reaction_summary",             "Enzyme & reactions"),
        ("domain_family_and_feature_summary",       "Domains, families & features"),
        ("pathway_and_context_summary",             "Pathways & context"),
        ("tool_prediction_summary",                 "Tool predictions"),
        ("strong_or_curated_information",           "Strong / curated information"),
        ("weak_predicted_or_indirect_information",  "Weak / predicted / indirect information"),
        ("conflicting_or_inconsistent_information", "Conflicting / inconsistent information"),
        ("missing_or_limited_information",          "Missing / limited information"),
    ]
    for key, label in list_sections:
        bullets = parsed.get(key) or []
        if not bullets:
            continue
        with st.expander(f"{label} ({len(bullets)})", expanded=False):
            for b in bullets:
                st.markdown(f"- {b}")

    # Review notes — always shown as a final highlighted list
    review = parsed.get("review_notes") or []
    if review:
        st.markdown("**📝 Review notes**")
        for r in review:
            st.markdown(f"- {r}")



# ── Header ────────────────────────────────────────────────────────────────────

st.markdown("""
<div class="hero">
    <h1>🧬 Protein Profiles Explorer</h1>
    <p>Interactive Streamlit dashboard for nested protein-profile JSON reports.</p>
</div>
""", unsafe_allow_html=True)

# ── Metric cards ──────────────────────────────────────────────────────────────

total   = len(df_all)
n_filt  = len(df_filtered)
n_orgs  = df_filtered["organism"].nunique()
n_rev   = (df_filtered["reviewed_status"] == "reviewed").sum()
n_seq   = df_filtered["seq_length"].notna().sum()
mean_sc = df_filtered["annotation_score"].mean()
mean_ln = df_filtered["seq_length"].mean()

# InterPro global statistics
n_with_interpro = 0
n_total_ipr     = 0
member_dbs_set  = set()
for p in proteins_filtered:
    ipro = p.get("interpro", {}) or {}
    summ = ipro.get("summary", {}) or {}
    if summ.get("total_entries", 0) > 0:
        n_with_interpro += 1
    n_total_ipr += summ.get("total_interpro_integrated", 0)
    for db in summ.get("member_databases_used", []) or []:
        member_dbs_set.add(db)

pct_interpro = (n_with_interpro / n_filt * 100) if n_filt > 0 else 0

# First row — UniProt-derived metrics
st.markdown("##### Dataset overview")
c1, c2, c3, c4, c5, c6 = st.columns(6)
for col, label, val, sub in [
    (c1, "Proteins",      f"{n_filt:,}", f"of {total:,} total"),
    (c2, "Organisms",     f"{n_orgs:,}", ""),
    (c3, "Reviewed",      f"{n_rev:,}",  "(SwissProt)"),
    (c4, "With sequence", f"{n_seq:,}",  ""),
    (c5, "Mean score",    f"{mean_sc:.2f}" if not math.isnan(mean_sc) else "—", "annotation"),
    (c6, "Mean length",   f"{int(mean_ln):,}" if not math.isnan(mean_ln or float('nan')) else "—", "aa"),
]:
    col.markdown(f"""
    <div class="metric-card">
        <div class="label">{label}</div>
        <div class="value">{val}</div>
        <div class="sub">{sub}</div>
    </div>""", unsafe_allow_html=True)

st.markdown("<br>", unsafe_allow_html=True)

# Second row — InterPro / Pipeline metrics
st.markdown("##### InterPro coverage")
i1, i2, i3 = st.columns(3)
for col, label, val, sub in [
    (i1, "With InterPro",    f"{n_with_interpro:,}", f"{pct_interpro:.1f}% of filtered"),
    (i2, "IPR entries",      f"{n_total_ipr:,}",     "across all proteins"),
    (i3, "Member DBs used",  f"{len(member_dbs_set):,}", ", ".join(sorted(member_dbs_set)[:5]) + ("…" if len(member_dbs_set) > 5 else "")),
]:
    col.markdown(f"""
    <div class="metric-card">
        <div class="label">{label}</div>
        <div class="value">{val}</div>
        <div class="sub">{sub}</div>
    </div>""", unsafe_allow_html=True)
    
st.markdown("<br>", unsafe_allow_html=True)

# ── Tabs ──────────────────────────────────────────────────────────────────────

tab_overview, tab_table, tab_detail, tab_export = st.tabs([
    "Overview", "Protein table", "Protein detail", "Data export"
])


# ══════════════════════════════════════════════════════════════════════════════
# TAB 1 — OVERVIEW
# ══════════════════════════════════════════════════════════════════════════════
with tab_overview:

    col_left, col_mid, col_right = st.columns(3)

    # Top organisms
    with col_left:
        st.markdown('<div class="section-title">Top organisms</div>', unsafe_allow_html=True)
        top_orgs = df_filtered["organism"].value_counts().head(10).reset_index()
        top_orgs.columns = ["organism", "count"]
        fig_org = px.bar(
            top_orgs, x="count", y="organism", orientation="h",
            color_discrete_sequence=["#3b82f6"],
            height=320,
        )
        fig_org.update_layout(margin=dict(l=0, r=0, t=0, b=0), yaxis_title="", xaxis_title="count",
                              plot_bgcolor="white", paper_bgcolor="white", font_size=11)
        fig_org.update_yaxes(autorange="reversed")
        st.plotly_chart(fig_org, use_container_width=True)

    # Confidence pie
    with col_mid:
        st.markdown('<div class="section-title">Confidence</div>', unsafe_allow_html=True)
        conf_counts = df_filtered["overall_confidence"].value_counts()
        fig_conf = go.Figure(go.Pie(
            labels=conf_counts.index,
            values=conf_counts.values,
            hole=0.45,
            marker_colors=["#3b82f6", "#93c5fd", "#ef4444"],
            textinfo="percent",
        ))
        fig_conf.update_layout(margin=dict(l=0, r=0, t=0, b=0), height=320,
                               showlegend=True, legend=dict(orientation="v"),
                               paper_bgcolor="white")
        st.plotly_chart(fig_conf, use_container_width=True)

    # Evidence tools
    with col_right:
        st.markdown('<div class="section-title">Evidence tools</div>', unsafe_allow_html=True)
        tool_counts: dict[str, int] = {}
        for p in proteins_filtered:
            for t in p.get("evidence_summary", {}).get("tools", []):
                tool_counts[t] = tool_counts.get(t, 0) + 1
        if tool_counts:
            df_tools = pd.DataFrame(sorted(tool_counts.items(), key=lambda x: x[1]), columns=["tool", "count"])
            fig_tools = px.bar(df_tools, x="count", y="tool", orientation="h",
                               color_discrete_sequence=["#6366f1"], height=320)
            fig_tools.update_layout(margin=dict(l=0, r=0, t=0, b=0), yaxis_title="", xaxis_title="count",
                                    plot_bgcolor="white", paper_bgcolor="white", font_size=11)
            st.plotly_chart(fig_tools, use_container_width=True)

    # Scatter: annotation score vs sequence length
    st.markdown('<div class="section-title">Annotation score vs. sequence length</div>', unsafe_allow_html=True)
    scatter_df = df_filtered[df_filtered["seq_length"].notna() & df_filtered["annotation_score"].notna()].copy()
    if not scatter_df.empty:
        fig_sc = px.scatter(
            scatter_df, x="seq_length", y="annotation_score",
            color="overall_confidence",
            color_discrete_map={"high": "#3b82f6", "medium": "#93c5fd", "low": "#ef4444", "unknown": "#9ca3af"},
            hover_data=["accession", "protein_name", "organism"],
            labels={"seq_length": "sequence_length", "annotation_score": "annotation_score"},
            height=320,
            opacity=0.7,
        )
        fig_sc.update_layout(margin=dict(l=0, r=0, t=0, b=0),
                             plot_bgcolor="white", paper_bgcolor="white", font_size=11)
        st.plotly_chart(fig_sc, use_container_width=True)




    # ── New section: InterPro analysis ─────────────────────────────────────
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown('<div class="section-title">InterPro analysis</div>', unsafe_allow_html=True)

    # Collect InterPro entries across all filtered proteins
    entry_counter: dict[str, dict] = {}  # accession -> info + count
    source_db_counter: dict[str, int] = {}
    type_counter: dict[str, int] = {}
    for p in proteins_filtered:
        ipro = p.get("interpro", {}) or {}
        for e in ipro.get("entries", []) or []:
            acc = e.get("accession")
            if not acc:
                continue
            if acc not in entry_counter:
                entry_counter[acc] = {
                    "accession":       acc,
                    "name":            e.get("name") or "—",
                    "source_database": e.get("source_database") or "unknown",
                    "type":            e.get("type") or "unknown",
                    "count":           0,
                }
            entry_counter[acc]["count"] += 1
            src = e.get("source_database") or "unknown"
            source_db_counter[src] = source_db_counter.get(src, 0) + 1
            tp = e.get("type") or "unknown"
            type_counter[tp] = type_counter.get(tp, 0) + 1

    col_ipro_l, col_ipro_r = st.columns(2)

    with col_ipro_l:
        st.markdown("**Top 20 entries (most frequent across dataset)**")
        if entry_counter:
            top_entries = sorted(entry_counter.values(), key=lambda x: -x["count"])[:20]
            df_top = pd.DataFrame(top_entries)
            df_top["label"] = df_top["accession"] + " — " + df_top["name"].str.slice(0, 40)
            fig_top = px.bar(
                df_top.iloc[::-1],
                x="count", y="label", orientation="h",
                color="source_database",
                hover_data=["accession", "name", "type", "source_database"],
                height=520,
            )
            fig_top.update_layout(margin=dict(l=0, r=0, t=0, b=0),
                                  yaxis_title="", xaxis_title="proteins",
                                  plot_bgcolor="white", paper_bgcolor="white",
                                  font_size=11, legend_title_text="Source DB")
            st.plotly_chart(fig_top, use_container_width=True)
        else:
            st.info("No InterPro entries in filtered proteins.")

    with col_ipro_r:
        st.markdown("**Source databases — annotations distribution**")
        if source_db_counter:
            df_src = pd.DataFrame(
                sorted(source_db_counter.items(), key=lambda x: -x[1]),
                columns=["source_database", "count"],
            )
            fig_src = px.bar(
                df_src.iloc[::-1],
                x="count", y="source_database", orientation="h",
                color_discrete_sequence=["#8b5cf6"],
                height=260,
            )
            fig_src.update_layout(margin=dict(l=0, r=0, t=0, b=0),
                                  yaxis_title="", xaxis_title="annotations",
                                  plot_bgcolor="white", paper_bgcolor="white",
                                  font_size=11)
            st.plotly_chart(fig_src, use_container_width=True)
        else:
            st.info("No InterPro source databases in filtered proteins.")

        st.markdown("**Entry types**")
        if type_counter:
            df_tp = pd.DataFrame(
                sorted(type_counter.items(), key=lambda x: -x[1]),
                columns=["type", "count"],
            )
            fig_tp = px.pie(
                df_tp, names="type", values="count",
                hole=0.4,
                height=220,
                color_discrete_sequence=px.colors.qualitative.Set2,
            )
            fig_tp.update_layout(margin=dict(l=0, r=0, t=0, b=0),
                                 paper_bgcolor="white", font_size=11,
                                 legend=dict(orientation="v", y=0.5))
            st.plotly_chart(fig_tp, use_container_width=True)

    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown('<div class="section-title">Stage 1 tool consensus</div>', unsafe_allow_html=True)
    st.caption("How many proteins were annotated by each combination of tools.")

    from collections import Counter as _Counter
    tools_per_protein = []
    for p in proteins_filtered:
        n_tools = len(p.get("evidence_summary", {}).get("tools", []))
        tools_per_protein.append(n_tools)
    if tools_per_protein:
        consensus_counts = _Counter(tools_per_protein)
        df_cons = pd.DataFrame(
            sorted(consensus_counts.items()),
            columns=["n_tools", "n_proteins"],
        )
        df_cons["label"] = df_cons["n_tools"].astype(str) + " tool(s)"
        fig_cons = px.bar(
            df_cons, x="label", y="n_proteins",
            color_discrete_sequence=["#10b981"],
            text="n_proteins",
            height=300,
        )
        fig_cons.update_traces(textposition="outside")
        fig_cons.update_layout(margin=dict(l=0, r=0, t=20, b=0),
                               xaxis_title="Number of Stage 1 tools that annotated",
                               yaxis_title="Number of proteins",
                               plot_bgcolor="white", paper_bgcolor="white",
                               font_size=11)
        st.plotly_chart(fig_cons, use_container_width=True)

# ══════════════════════════════════════════════════════════════════════════════
# TAB 2 — PROTEIN TABLE
# ══════════════════════════════════════════════════════════════════════════════
with tab_table:
    st.markdown("### Filtered proteins")
    st.markdown('<div class="tip-box">💡 Tip: use the sidebar filters, then click column headers to sort the table.</div>',
                unsafe_allow_html=True)

    display_cols = ["accession", "protein_name", "gene_name", "organism", "reviewed_status", "overall_confidence"]
    st.dataframe(
        df_filtered[display_cols].reset_index(drop=True),
        use_container_width=True,
        height=520,
        column_config={
            "accession":          st.column_config.TextColumn("Accession"),
            "protein_name":       st.column_config.TextColumn("Protein name"),
            "gene_name":          st.column_config.TextColumn("Gene name"),
            "organism":           st.column_config.TextColumn("Organism"),
            "reviewed_status":    st.column_config.TextColumn("Reviewed status"),
            "overall_confidence": st.column_config.TextColumn("Overall confidence"),
        }
    )

    st.markdown("""
    <div class="tip-box">
        Click a row in the table to reveal its origin, metadata, evidence summary, and database cross-references.
    </div>""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════════════════
# TAB 3 — PROTEIN DETAIL
# ══════════════════════════════════════════════════════════════════════════════
with tab_detail:
    st.markdown("### Inspect one protein")

    if "sel_acc" not in st.session_state:
        st.session_state["sel_acc"] = ""

    acc_list = [""] + sorted(df_filtered["accession"].tolist())
    prev = st.session_state["sel_acc"]
    idx  = acc_list.index(prev) if prev in acc_list else 0

    sel_acc = st.selectbox("Protein accession", acc_list, index=idx, key="sel_acc_box")
    st.session_state["sel_acc"] = sel_acc

    if not sel_acc:
        st.info("Select a protein accession above to inspect its full profile.")
    else:
        prot = next((p for p in all_proteins if p["accession"] == sel_acc), None)
        if prot is None:
            st.error("Protein not found.")
        else:
            ident = prot.get("identity", {})
            ev    = prot.get("evidence_summary", {})
            seq   = prot.get("sequence") or {}

            # Header
            pname = ident.get("protein_name") or "Uncharacterized protein"
            ename = ident.get("entry_name", "")
            st.markdown(f"## {pname}")
            st.markdown(f"<small>Accession: **{sel_acc}** · Entry: **{ename}**</small>", unsafe_allow_html=True)

            # Quick metric cards
            qc1, qc2, qc3, qc4 = st.columns(4)
            for col, label, val in [
                (qc1, "Reviewed status",    ident.get("reviewed_status", "—")),
                (qc2, "Annotation score",   str(ident.get("annotation_score", "—"))),
                (qc3, "Sequence length",    f"{seq.get('length', '—')} aa"),
                (qc4, "Overall confidence", ev.get("overall_confidence", "—")),
            ]:
                col.markdown(f"""
                <div class="metric-card">
                    <div class="label">{label}</div>
                    <div class="value" style="font-size:20px;">{val}</div>
                </div>""", unsafe_allow_html=True)

            uniprot_url = prot.get("provenance", {}).get("uniprot_url", "")
            if uniprot_url:
                st.link_button("Open UniProt entry", uniprot_url)

            st.markdown("<br>", unsafe_allow_html=True)

            # Stage 3 LLM summary (only for proteins in the test set)
            _stage3_index = load_stage3_results()
            _stage3_record = _stage3_index.get(sel_acc) if _stage3_index else None
            if _stage3_record:
                _llm_summary_render_section(_stage3_record)
                st.markdown("<br>", unsafe_allow_html=True)
            # Sub-tabs
            st_sum, st_go, st_dom, st_ipro, st_seq, st_refs, st_origin, st_raw = st.tabs([
                "Summary", "GO & enzyme", "Domains & pathways", "InterPro",
                "Sequence", "References & xrefs", "Origin/meta", "Raw JSON"
            ])

            # ── Summary ───────────────────────────────────────────────────
            with st_sum:
                left, right = st.columns(2)

                with left:
                    st.markdown("**Identity**")
                    id_items = {
                        "Gene name":        ident.get("gene_name") or (
                            f"{(ident.get('gene_name_inferred') or {}).get('value', '—')} *(inferred)*"),
                        "Protein existence": ident.get("protein_existence", "—"),
                        "Entry version":    ident.get("entry_version", "—"),
                        "First public":     ident.get("first_public_date", "—"),
                        "Last updated":     ident.get("last_annotation_update", "—"),
                    }
                    for k, v in id_items.items():
                        st.markdown(f"- **{k}:** {v}")

                    org = ident.get("organism", {})
                    if org.get("scientific_name"):
                        st.markdown("**Organism**")
                        st.markdown(f"- **Name:** *{org['scientific_name']}*")
                        if org.get("common_name"):
                            st.markdown(f"- **Common:** {org['common_name']}")
                        st.markdown(f"- **TaxID:** {org.get('taxon_id', '—')}")
                        if org.get("lineage"):
                            st.markdown(f"- **Lineage:** {' › '.join(org['lineage'])}")
                        if org.get("evidences"):
                            st.markdown("- **Organism evidences:**")
                            for ev in org["evidences"]:
                                st.markdown(f'  - `{ev.get("code","")}` · {ev.get("source","")} · `{ev.get("id","")}`')

                    fn = prot.get("function", {})
                    if fn.get("description"):
                        st.markdown("**Function description**")
                        for d in fn["description"]:
                            st.markdown(f"> {d.get('value', '')}")

                    if fn.get("similarity"):
                        st.markdown("**Similarity**")
                        for s in fn["similarity"]:
                            st.markdown(f"- {s.get('value', '')}")

                with right:
                    fn = prot.get("function", {})
                    if fn.get("keywords"):
                        st.markdown("**Keywords**")
                        for kw in fn["keywords"]:
                            if isinstance(kw, dict):
                                kw_name = kw.get("name", "")
                                kw_id   = kw.get("id", "")
                                kw_cat  = kw.get("category", "")
                                kw_evs  = kw.get("evidences", [])
                                ev_str  = " ".join(f"`{e.get('code','')}:{e.get('source','')}:{e.get('id','')}`" for e in kw_evs)
                                st.markdown(f'<span class="badge badge-gray">{kw_name}</span> <small>`{kw_id}` · {kw_cat} {ev_str}</small>', unsafe_allow_html=True)
                            else:
                                st.markdown(f'<span class="badge badge-gray">{kw}</span>', unsafe_allow_html=True)

                    st.markdown("**Tools / sources**")
                    tools_html = " ".join(f'<span class="badge">{t}</span>' for t in ev.get("tools", []))
                    st.markdown(tools_html, unsafe_allow_html=True)

                    st.markdown("**Evidence confidence counts**")
                    cd = ev.get("confidence_distribution", {})
                    for level, color in [("high", "green"), ("medium", "yellow"), ("low", "red")]:
                        count = cd.get(level, 0)
                        st.markdown(f'<span class="badge badge-{color}">{level}: {count}</span>', unsafe_allow_html=True)

                    if fn.get("subcellular_location"):
                        st.markdown("**Subcellular location**")
                        for loc in fn["subcellular_location"]:
                            lv = loc.get("location", "")
                            lid = loc.get("location_id", "")
                            st.markdown(f"- {lv}" + (f" `{lid}`" if lid else ""))

                    if fn.get("pathway"):
                        st.markdown("**Pathways (UniProt)**")
                        for pw in fn["pathway"]:
                            st.markdown(f"- {pw.get('value', '')}")

                    if fn.get("subunit"):
                        st.markdown("**Subunit structure**")
                        for su in fn["subunit"]:
                            st.markdown(f"- {su.get('value', '')}")

                    if fn.get("catalytic_activity"):
                        st.markdown("**Catalytic activity**")
                        for ca in fn["catalytic_activity"]:
                            ec = ca.get("ec_number", "")
                            rxn = ca.get("reaction", "")
                            rhea = ", ".join(ca.get("rhea_ids", []))
                            st.markdown(f"- **EC:** `{ec}` — {rxn}")
                            if rhea:
                                st.markdown(f"  - Rhea: {rhea}")

                    if prot.get("evidence_summary", {}).get("in_poorly_annotated_subset"):
                        st.warning("⚠️ This protein is in the **poorly annotated** subset.")

            # ── GO & Enzyme ───────────────────────────────────────────────
            with st_go:
                go_ann = prot.get("go_annotations", {})

                for aspect, label, color in [
                    ("molecular_function", "Molecular Function", "#0ea5e9"),
                    ("biological_process", "Biological Process", "#22c55e"),
                    ("cellular_component", "Cellular Component", "#a855f7"),
                ]:
                    terms = go_ann.get(aspect, [])
                    if not terms:
                        continue
                    st.markdown(f"**{label}** ({len(terms)} terms)")
                    rows_data = []
                    for t in terms:
                        sd = t.get("source_details", {})
                        sources_str = ", ".join(
                            f"{s} ({sd[s].get('source_type','')}, score={sd[s].get('score','—')})"
                            if s in sd else s
                            for s in t.get("sources", [])
                        )
                        rows_data.append({
                            "GO ID":      t.get("go_id", ""),
                            "Label":      t.get("label") or "—",
                            "Confidence": t.get("confidence", "—"),
                            "Sources":    sources_str,
                            "Evidence":   f"{t.get('evidence_code','—')} / {t.get('evidence_source','—')}",
                        })
                    st.dataframe(pd.DataFrame(rows_data), use_container_width=True, hide_index=True, height=200)

                ec_numbers = (prot.get("enzymatic") or {}).get("ec_numbers", [])
                if ec_numbers:
                    st.markdown("**EC Numbers**")
                    ec_rows = []
                    for ec in ec_numbers:
                        sd = ec.get("source_details", {})
                        sources_str = ", ".join(
                            f"{s} ({sd[s].get('source_type','')})" if s in sd else s
                            for s in ec.get("sources", [])
                        )
                        ec_rows.append({
                            "EC ID":      ec.get("ec_id", ""),
                            "Confidence": ec.get("confidence", "—"),
                            "Sources":    sources_str,
                        })
                    st.dataframe(pd.DataFrame(ec_rows), use_container_width=True, hide_index=True)

            # ── Domains & Pathways ─────────────────────────────────────────
            with st_dom:
                left_d, right_d = st.columns(2)

                with left_d:
                    doms = prot.get("domains", {})
                    for dtype in ["pfam", "cog", "kog", "tigrfam", "smart"]:
                        items = doms.get(dtype, [])
                        if items:
                            st.markdown(f"**{dtype.upper()}**")
                            for item in items:
                                label = item.get("label") or ""
                                conf  = item.get("confidence", "")
                                srcs  = ", ".join(item.get("sources", []))
                                st.markdown(f"- `{item['id']}` {label} — *{conf}* ({srcs})")

                    cog_cats = doms.get("cog_categories", [])
                    if cog_cats:
                        st.markdown("**COG categories**")
                        st.markdown(", ".join(f"`{c}`" for c in cog_cats))

                    feats = prot.get("features", {})
                    for ftype, label in [("domains", "UniProt domains"), ("active_sites", "Active sites"),
                                         ("binding_sites", "Binding sites"), ("other_features", "Other features")]:
                        items = feats.get(ftype, [])
                        if items:
                            st.markdown(f"**{label}**")
                            for item in items:
                                desc   = item.get("description") or item.get("type", "")
                                start  = item.get("position_start", "?")
                                end    = item.get("position_end", "?")
                                s_mod  = item.get("position_start_modifier", "")
                                e_mod  = item.get("position_end_modifier", "")
                                fid    = item.get("feature_id", "")
                                pos_str = f"[{start}–{end}]"
                                if s_mod and s_mod != "EXACT":
                                    pos_str = f"[{s_mod} {start}–{e_mod} {end}]"
                                evs = item.get("evidences", [])
                                ev_str = " ".join(f"`{e.get('code','')}:{e.get('source','')}:{e.get('id','')}`" for e in evs)
                                st.markdown(f"- {desc} {pos_str}" + (f" `{fid}`" if fid else "") + (f" {ev_str}" if ev_str else ""))

                with right_d:
                    pw = prot.get("pathways", {})
                    if pw.get("kegg_ko"):
                        st.markdown("**KEGG KO**")
                        for k in pw["kegg_ko"]:
                            st.markdown(f"- `{k}`")
                    if pw.get("kegg_pathways"):
                        st.markdown("**KEGG Pathways**")
                        for k in pw["kegg_pathways"]:
                            st.markdown(f"- `{k}`")
                    if pw.get("kegg_modules"):
                        st.markdown("**KEGG Modules**")
                        for k in pw["kegg_modules"]:
                            st.markdown(f"- `{k}`")

            # ── InterPro ──────────────────────────────────────────────────
            with st_ipro:
                interpro_data = prot.get("interpro", {}) or {}
                _interpro_render_section(interpro_data)

            # ── Sequence ──────────────────────────────────────────────────
            with st_seq:
                if seq:
                    s1, s2, s3 = st.columns(3)
                    s1.metric("Length", f"{seq.get('length', '—')} aa")
                    s2.metric("Mol. weight", f"{seq.get('mol_weight', '—')} Da")
                    s3.metric("CRC64", seq.get("crc64", "—"))
                    if seq.get("md5"):
                        st.caption(f"MD5: `{seq['md5']}`")
                    if seq.get("value"):
                        st.markdown("**Amino acid sequence**")
                        st.code(seq["value"], language=None)
                else:
                    st.info("No sequence data available for this protein.")

            # ── References & Xrefs ────────────────────────────────────────
            with st_refs:
                refs = prot.get("references", [])
                if refs:
                    st.markdown("**References**")
                    for ref in refs:
                        num   = ref.get("reference_number", "")
                        title = ref.get("title") or "*(no title)*"
                        auth  = ", ".join(ref.get("authors", [])[:3])
                        if len(ref.get("authors", [])) > 3:
                            auth += " et al."
                        date  = ref.get("date", "")
                        pmid  = ref.get("pmid")
                        sub_db = ref.get("submission_database", "")
                        ref_comments = ref.get("reference_comments", [])

                        with st.expander(f"[{num}] {title[:80]}…" if len(title) > 80 else f"[{num}] {title}"):
                            if auth:   st.markdown(f"**Authors:** {auth}")
                            auth_group = ref.get("authoring_group", [])
                            if auth_group: st.markdown(f"**Authoring group:** {', '.join(auth_group)}")
                            if date:   st.markdown(f"**Date:** {date}")
                            if pmid:   st.markdown(f"**PMID:** [{pmid}](https://pubmed.ncbi.nlm.nih.gov/{pmid})")
                            if sub_db: st.markdown(f"**Submission DB:** {sub_db}")
                            for rc in ref_comments:
                                st.markdown(f"**{rc.get('type','')}:** {rc.get('value','')}")
                                for ev in rc.get("evidences", []):
                                    st.caption(f"{ev.get('code','')} · {ev.get('source','')} · {ev.get('id','')}")
                            for pos in ref.get("reference_positions", []):
                                st.caption(pos)
                            ref_evs = ref.get("evidences", [])
                            if ref_evs:
                                st.markdown("**Reference evidences:**")
                                for ev in ref_evs:
                                    st.caption(f"{ev.get('code','')} · {ev.get('source','')} · {ev.get('id','')}") 

                xrefs = prot.get("cross_references", {})
                if xrefs:
                    st.markdown("**Cross-references**")
                    xref_cols = st.columns(3)
                    for i, (db, items) in enumerate(sorted(xrefs.items())):
                        col = xref_cols[i % 3]
                        col.markdown(f"**{db}**")
                        for item in items:
                            xid = item.get("id", item) if isinstance(item, dict) else item
                            props = item.get("properties", {}) if isinstance(item, dict) else {}
                            prop_str = " · ".join(f"{k}: {v}" for k, v in props.items() if v and v != "-")
                            col.markdown(f"- `{xid}`" + (f" {prop_str}" if prop_str else ""))

            # ── Origin / Meta ─────────────────────────────────────────────
            with st_origin:
                st.markdown("### Origin, metadata & cross-references")
                st.caption("This panel is tied to the currently selected protein record.")

                o1, o2, o3 = st.columns(3)
                prov = prot.get("provenance", {})
                o1.markdown("**Origin**")
                o1.markdown(f"`{prov.get('source', '—')}`")
                o2.markdown("**DB source**")
                o2.markdown(f"`{prov.get('db_source', '—')}`")
                o3.markdown("**Cross-ref DBs**")
                xrefs = prot.get("cross_references", {})
                o3.markdown(", ".join(f"`{db}`" for db in sorted(xrefs.keys())) or "—")

                st.markdown("**Entry type**")
                st.markdown(f"`{ident.get('entry_type', '—')}`")

                st.markdown("**Entry audit**")
                a1, a2, a3, a4 = st.columns(4)
                a1.metric("Entry version",    ident.get("entry_version", "—"))
                a2.metric("Sequence version", ident.get("sequence_version", "—"))
                a3.metric("First public",     ident.get("first_public_date", "—"))
                a4.metric("Last updated",     ident.get("last_annotation_update", "—"))

                if ident.get("uniprot_id"):
                    st.markdown(f"**UniParc ID:** `{ident['uniprot_id']}`")

                if ident.get("alternative_names"):
                    st.markdown("**Alternative protein names**")
                    for alt in ident["alternative_names"]:
                        name = alt.get("name", "")
                        evs  = alt.get("evidences", [])
                        ev_str = " ".join(f"`{e.get('code','')}:{e.get('source','')}:{e.get('id','')}`" for e in evs)
                        st.markdown(f"- {name} {ev_str}")

                if ident.get("protein_name_evidences"):
                    st.markdown("**Protein name evidences**")
                    for ev in ident["protein_name_evidences"]:
                        st.markdown(f'<div class="evidence-row"><span class="eco">{ev.get("code","")}</span> · {ev.get("source","")} · {ev.get("id","")}</div>',
                                    unsafe_allow_html=True)

                if ident.get("gene_name_synonym_details"):
                    st.markdown("**ORF names / gene synonyms**")
                    for g in ident["gene_name_synonym_details"]:
                        evs = g.get("evidences", [])
                        ev_str = " ".join(f"`{e.get('code','')}:{e.get('source','')}:{e.get('id','')}`" for e in evs)
                        st.markdown(f"- `{g.get('value','')}` ({g.get('type','')}) {ev_str}")

                extra = prot.get("extra_attributes")
                if extra:
                    st.markdown("**Extra attributes**")
                    st.json(extra)

            # ── Raw JSON ──────────────────────────────────────────────────
            with st_raw:
                st.markdown("**Raw protein profile JSON**")
                st.json(prot, expanded=False)


# ══════════════════════════════════════════════════════════════════════════════
# TAB 4 — DATA EXPORT
# ══════════════════════════════════════════════════════════════════════════════
with tab_export:
    st.markdown("### Export filtered data")

    export_cols = ["accession", "protein_name", "gene_name", "organism", "taxon_id",
                   "reviewed_status", "annotation_score", "overall_confidence",
                   "annotation_count", "n_go", "n_ec", "n_tools", "poorly_annotated",
                   "seq_length", "mol_weight", "conf_high", "conf_medium", "conf_low"]

    st.dataframe(df_filtered[export_cols].reset_index(drop=True), use_container_width=True, height=400)

    csv = df_filtered[export_cols].to_csv(index=False).encode("utf-8")
    st.download_button(
        label="⬇️ Download as CSV",
        data=csv,
        file_name="protein_profiles_filtered.csv",
        mime="text/csv",
    )

    json_out = json.dumps(proteins_filtered, indent=2, ensure_ascii=False).encode("utf-8")
    st.download_button(
        label="⬇️ Download filtered JSON",
        data=json_out,
        file_name="protein_profiles_filtered.json",
        mime="application/json",
    )
