"""
Protein Profiles Explorer — Streamlit dashboard
Reads protein_profiles.json and provides interactive exploration.

Run with:
    pip install streamlit pandas plotly
    streamlit run app.py
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
        border-radius: 8px;
        padding: 16px 20px;
        text-align: center;
    }
    .metric-card .label { font-size: 12px; color: #6c757d; text-transform: uppercase; letter-spacing: 0.05em; margin-bottom: 4px; }
    .metric-card .value { font-size: 28px; font-weight: 600; color: #212529; font-family: 'IBM Plex Mono', monospace; }
    .metric-card .sub   { font-size: 11px; color: #28a745; margin-top: 2px; }
    
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
        for candidate in ["protein_profiles.json", "protein_profiles.zip"]:
            if Path(candidate).exists():
                bundle_name = candidate
                with open(candidate, "rb") as f:
                    proteins_raw = f.read()
                st.success(f"Loaded bundled file:\n{candidate}")
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

c1, c2, c3, c4, c5, c6 = st.columns(6)
for col, label, val, sub in [
    (c1, "Proteins",      f"{n_filt:,}", f"↑ of {total:,}"),
    (c2, "Organisms",     f"{n_orgs:,}", ""),
    (c3, "Reviewed",      f"{n_rev:,}",  ""),
    (c4, "With sequence", f"{n_seq:,}",  ""),
    (c5, "Mean score",    f"{mean_sc:.2f}" if not math.isnan(mean_sc) else "—", ""),
    (c6, "Mean length",   f"{int(mean_ln):,} aa" if not math.isnan(mean_ln or float('nan')) else "—", ""),
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

    search_acc = st.text_input("Search accession or protein name", placeholder="e.g. X5DZ82 or SusD", key="search_acc")
    if search_acc:
        q = search_acc.lower()
        filtered_accs = [
            a for a in df_filtered["accession"].tolist()
            if q in a.lower() or q in (df_filtered.loc[df_filtered["accession"] == a, "protein_name"].values[0] or "").lower()
        ]
        acc_list = [""] + sorted(filtered_accs)

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

            # Sub-tabs
            st_sum, st_go, st_dom, st_seq, st_refs, st_origin, st_raw = st.tabs([
                "Summary", "GO & enzyme", "Domains & pathways",
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
