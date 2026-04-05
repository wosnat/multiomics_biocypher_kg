"""UI components for the cluster extraction review app."""

import re
import subprocess
import webbrowser
from datetime import datetime
from pathlib import Path
from typing import Optional

import pandas as pd
import streamlit as st
import yaml

from multiomics_kg.extraction.cluster.run_manager import RunManager

# Synthesis fields shown below the 3-column view
SYNTHESIS_FIELDS = [
    "id", "name", "functional_description",
    "behavioral_description", "peak_time_hours", "period_hours",
]

# Confidence level color badges
CONFIDENCE_COLORS = {
    "very_high": "#2e7d32",  # green
    "high": "#1565c0",  # blue
    "medium": "#f57f17",  # amber
    "low": "#c62828",  # red
}

# Source label colors
SOURCE_COLORS = {
    "table": "#1565c0",
    "visual": "#6a1b9a",
    "semantic": "#e65100",
}

# Review status options
STATUS_OPTIONS = ["unreviewed", "approve", "edit", "reject", "flag-issue"]

# Issue classification options
ISSUE_OPTIONS = [
    "wrong_cluster",
    "hallucinated",
    "low_info",
    "cross_contamination",
    "partial",
    "source_missing",
]

# Failing stage options
FAILING_STAGE_OPTIONS = [
    "stage1_visual",
    "stage1_semantic",
    "stage1_table",
    "stage2_synthesis",
    "stage3_validation",
    "merge",
]

# Verdict badge colors
VERDICT_COLORS = {
    "pass": ("#2e7d32", "white"),
    "warn": ("#f57f17", "black"),
    "fail": ("#c62828", "white"),
    "none": ("#9e9e9e", "white"),
}


def _confidence_badge(confidence: str) -> str:
    """Return HTML badge for confidence level."""
    color = CONFIDENCE_COLORS.get(confidence, "#9e9e9e")
    return f'<span style="background:{color};color:white;padding:1px 6px;border-radius:3px;font-size:0.75em;">{confidence}</span>'


def _source_badge(source: str) -> str:
    """Return HTML badge for source type."""
    color = SOURCE_COLORS.get(source, "#616161")
    return f'<span style="background:{color};color:white;padding:1px 6px;border-radius:3px;font-size:0.75em;">{source}</span>'


def _verdict_badge(verdict: str) -> str:
    """Return HTML badge for verdict."""
    bg, fg = VERDICT_COLORS.get(verdict, ("#9e9e9e", "white"))
    return f'<span style="background:{bg};color:{fg};padding:2px 8px;border-radius:4px;font-size:0.85em;font-weight:bold;">{verdict}</span>'


def _render_path_fields(fields: dict) -> None:
    """Render extracted fields for one path."""
    if not fields:
        st.caption("No data from this path.")
        return
    for field, value in sorted(fields.items()):
        if value and value != "not described in paper":
            label = field.replace("_", " ").title()
            st.markdown(f"**{label}:** {value}")


def _render_path_quotes(stage1_cluster: dict, source: str) -> None:
    """Render supporting quotes from a specific source."""
    quotes = stage1_cluster.get("supporting_quotes", [])
    source_quotes = [q for q in quotes if isinstance(q, dict) and q.get("source") == source]
    paper_pdf_path = st.session_state.get("paper_pdf_path")
    if source_quotes:
        with st.expander(f"Quotes ({len(source_quotes)})"):
            for q in source_quotes:
                text = q.get("quote", "")
                loc = q.get("location", "")
                st.markdown(f"> {text[:400]}{'...' if len(text) > 400 else ''}")
                if loc:
                    st.caption(f"Location: {loc}")
                    pdf_link = _pdf_page_link(loc, paper_pdf_path)
                    if pdf_link:
                        st.markdown(pdf_link)


def render_merge_view(cluster_key: str, stage1_cluster: dict, stage2_cluster: dict) -> None:
    """Render 3-column path view for one cluster, then synthesis result below."""
    paper_pdf_path = st.session_state.get("paper_pdf_path")

    # Separate stage1 data by source
    path_data = {"table": {}, "visual": {}, "semantic": {}}

    for field, entries in stage1_cluster.items():
        if field in ("supporting_quotes", "gene_count", "retrieved_passages"):
            continue
        if isinstance(entries, list):
            for entry in entries:
                if isinstance(entry, dict) and "source" in entry:
                    source = entry["source"]
                    if source in path_data:
                        path_data[source][field] = entry.get("value", "")

    # 3 columns
    col_table, col_visual, col_semantic = st.columns(3)

    with col_table:
        st.markdown("**Table** (very_high)")
        _render_path_fields(path_data["table"])
        _render_path_quotes(stage1_cluster, "table")

    with col_visual:
        st.markdown("**Visual** (high)")
        _render_path_fields(path_data["visual"])
        _render_path_quotes(stage1_cluster, "visual")

    with col_semantic:
        st.markdown("**Semantic** (medium)")
        _render_path_fields(path_data["semantic"])
        # Show retrieved RAG passages
        retrieved = stage1_cluster.get("retrieved_passages", [])
        if retrieved:
            with st.expander(f"Retrieved passages ({len(retrieved)})"):
                for p in retrieved:
                    score = p.get("relevance_score", "?")
                    text = p.get("text", "")
                    loc = p.get("location", "")
                    st.markdown(f"**[{score}]**")
                    st.markdown(f"> {text[:500]}{'...' if len(text) > 500 else ''}")
                    if loc:
                        pdf_link = _pdf_page_link(loc, paper_pdf_path)
                        if pdf_link:
                            st.markdown(pdf_link)
        else:
            # Fallback: show semantic supporting_quotes
            _render_path_quotes(stage1_cluster, "semantic")

    # Synthesis result below
    st.markdown("---")
    st.markdown("#### Synthesis Result")
    for field in SYNTHESIS_FIELDS:
        value = stage2_cluster.get(field, "")
        label = field.replace("_", " ").title()
        if value is not None and value != "":
            st.markdown(f"**{label}:** {value}")
        else:
            st.markdown(f"**{label}:** *(empty)*")


def _pdf_page_link(location: str, paper_pdf_path: Optional[str]) -> Optional[str]:
    """If location looks like a page number and PDF path is available, return a link."""
    if not location or not paper_pdf_path:
        return None
    page_match = re.search(r'(\d+)', str(location))
    if page_match:
        page_num = page_match.group(1)
        return f"[Open PDF page {page_num}](file://{paper_pdf_path}#page={page_num})"
    return None


def render_quotes(cluster_key: str, stage1_cluster: dict) -> None:
    """Render supporting quotes with source labels."""
    quotes = stage1_cluster.get("supporting_quotes", [])
    if not quotes:
        st.info("No supporting quotes available.")
        return

    paper_pdf_path = st.session_state.get("paper_pdf_path")

    for i, q in enumerate(quotes):
        if isinstance(q, dict):
            quote_text = q.get("quote", "")
            source = q.get("source", "")
            location = q.get("location", "")
            relevance = q.get("relevance_score", "")

            header_parts = []
            if source:
                header_parts.append(_source_badge(source))
            if location:
                header_parts.append(f"Location: {location}")
            if relevance:
                header_parts.append(f"Relevance: {relevance:.2f}" if isinstance(relevance, float) else f"Relevance: {relevance}")

            if header_parts:
                st.markdown(" ".join(header_parts), unsafe_allow_html=True)
            st.markdown(f"> {quote_text[:500]}{'...' if len(quote_text) > 500 else ''}")

            # PDF page link
            pdf_link = _pdf_page_link(location, paper_pdf_path)
            if pdf_link:
                st.markdown(pdf_link)
        else:
            st.markdown(f"> {q}")

        if i < len(quotes) - 1:
            st.markdown("")


def render_verdict(cluster_key: str, stage3_cluster: dict) -> None:
    """Render stage 3 validation verdict with badge."""
    if not stage3_cluster:
        st.warning("No validation data available.")
        return

    verdict = stage3_cluster.get("verdict", "none")
    explanation = stage3_cluster.get("explanation", "")

    st.markdown(f"Verdict: {_verdict_badge(verdict)}", unsafe_allow_html=True)

    # Show individual checks
    checks = [
        ("enrichment_correct", "Enrichment"),
        ("genes_correct", "Genes"),
        ("direction_correct", "Direction"),
        ("behavioral_correct", "Behavioral"),
        ("hallucination", "Hallucination"),
        ("missing_info", "Missing info"),
    ]
    check_parts = []
    for key, label in checks:
        val = stage3_cluster.get(key, "")
        if val:
            icon = {"yes": "Y", "no": "N", "uncertain": "?", "none": "-"}.get(val, val)
            check_parts.append(f"{label}: {icon}")
    if check_parts:
        st.markdown("&nbsp;&nbsp;" + " | ".join(check_parts))

    if explanation:
        st.markdown(f"*{explanation}*")


def render_diff_view(rm: RunManager, cluster_key: str) -> None:
    """Render a diff view comparing current run with a previous run for one cluster."""
    runs = rm.list_runs()
    if len(runs) < 2:
        st.info("No previous run to compare.")
        return

    current_run = runs[-1]
    previous_runs = runs[:-1]

    # Selectbox to pick which previous run to compare against
    prev_labels = [r.name for r in previous_runs]
    default_idx = len(prev_labels) - 1  # most recent previous
    selected_label = st.selectbox(
        "Compare with run",
        prev_labels,
        index=default_idx,
        key=f"diff_run_{cluster_key}",
    )
    selected_prev = previous_runs[prev_labels.index(selected_label)]

    # Load stage 2 data from both runs
    current_stage2 = rm.read_stage(current_run, 2)
    prev_stage2 = rm.read_stage(selected_prev, 2)

    current_cluster = current_stage2.get(cluster_key, {})
    prev_cluster = prev_stage2.get(cluster_key, {})

    diff_fields = [
        "id", "name", "functional_description",
        "behavioral_description", "peak_time_hours", "period_hours",
    ]

    for field in diff_fields:
        old_val = str(prev_cluster.get(field, ""))
        new_val = str(current_cluster.get(field, ""))
        changed = old_val != new_val

        label = field.replace("_", " ").title()
        if changed:
            st.markdown(
                f'<span style="color:red;font-weight:bold;">CHANGED - {label}</span>',
                unsafe_allow_html=True,
            )
            col_old, col_new = st.columns(2)
            with col_old:
                st.markdown(f"**Old:** {old_val or '(empty)'}")
            with col_new:
                st.markdown(f"**New:** {new_val or '(empty)'}")
        else:
            st.markdown(
                f'<span style="color:gray;">{label}: {new_val or "(empty)"}</span>',
                unsafe_allow_html=True,
            )


def render_review_controls(
    cluster_key: str,
    stage2_cluster: dict,
    current_review: dict,
    run_dir: Path,
    rm: RunManager,
    stage1: dict,
) -> None:
    """Render review controls: status, issues, notes, editable fields, save."""
    prefix = f"review_{cluster_key}"

    # Current values
    cur_status = current_review.get("status", "unreviewed")
    cur_issues = current_review.get("issues", [])
    cur_failing = current_review.get("failing_stages", [])
    cur_notes = current_review.get("notes", "")
    cur_edits = current_review.get("edited_fields", {})

    col1, col2 = st.columns(2)
    with col1:
        status = st.selectbox(
            "Review status",
            STATUS_OPTIONS,
            index=STATUS_OPTIONS.index(cur_status) if cur_status in STATUS_OPTIONS else 0,
            key=f"{prefix}_status",
        )
    with col2:
        issues = st.multiselect(
            "Issues",
            ISSUE_OPTIONS,
            default=[i for i in cur_issues if i in ISSUE_OPTIONS],
            key=f"{prefix}_issues",
        )

    failing_stages = st.multiselect(
        "Failing stages",
        FAILING_STAGE_OPTIONS,
        default=[f for f in cur_failing if f in FAILING_STAGE_OPTIONS],
        key=f"{prefix}_failing",
    )

    notes = st.text_area(
        "Notes",
        value=cur_notes,
        key=f"{prefix}_notes",
    )

    # Editable fields only when status is "edit"
    edited_fields = {}
    if status == "edit":
        st.markdown("**Edit description fields:**")
        editable_keys = ["functional_description", "behavioral_description", "id"]
        for fk in editable_keys:
            default_val = cur_edits.get(fk, stage2_cluster.get(fk, ""))
            edited_fields[fk] = st.text_area(
                fk.replace("_", " ").title(),
                value=default_val,
                key=f"{prefix}_edit_{fk}",
            )

    # Re-run buttons
    rerun_col1, rerun_col2, rerun_col3 = st.columns(3)
    with rerun_col1:
        if st.button("Re-synthesize (Stage 2+3)", key=f"{prefix}_resynth"):
            paperconfig_path = st.session_state.get("paperconfig_path")
            entry_key = st.session_state.get("entry_key")
            if paperconfig_path and entry_key:
                from multiomics_kg.extraction.cluster.pipeline import run_pipeline
                with st.spinner("Re-running stages 2+3..."):
                    run_pipeline(Path(paperconfig_path), table_key=entry_key, from_stage=2)
                st.success("Re-synthesis complete.")
                st.rerun()
            else:
                st.error("Missing paperconfig_path or entry_key in session state.")
    with rerun_col2:
        if st.button("Full re-extract", key=f"{prefix}_full_reextract"):
            paperconfig_path = st.session_state.get("paperconfig_path")
            entry_key = st.session_state.get("entry_key")
            if paperconfig_path and entry_key:
                from multiomics_kg.extraction.cluster.pipeline import run_pipeline
                with st.spinner("Running full re-extraction..."):
                    run_pipeline(Path(paperconfig_path), table_key=entry_key)
                st.success("Full re-extraction complete.")
                st.rerun()
            else:
                st.error("Missing paperconfig_path or entry_key in session state.")
    with rerun_col3:
        pass  # spacer

    # Save button
    if st.button("Save review", key=f"{prefix}_save", type="primary"):
        # Read existing stage4 data
        stage4 = rm.read_stage(run_dir, 4)
        input_hash = rm.compute_input_hash(cluster_key, stage1)

        review_entry = {
            "status": status,
            "issues": issues,
            "failing_stages": failing_stages,
            "notes": notes,
            "input_hash": input_hash,
            "reviewed_in_run": run_dir.name,
            "reviewed_at": datetime.now().isoformat(timespec="seconds"),
        }
        if status == "edit" and edited_fields:
            review_entry["edited_fields"] = edited_fields

        stage4[cluster_key] = review_entry
        rm.write_stage(run_dir, 4, stage4)
        st.success(f"Saved review for cluster {cluster_key}")
        st.rerun()


def render_source_access(paper_dir: Path, paperconfig_path: Path) -> None:
    """Render source file access: inline tables for CSV/XLSX, open links for PDFs."""
    if not paper_dir.exists():
        st.warning(f"Paper directory not found: {paper_dir}")
        return

    # Collect files, excluding hidden dirs and __pycache__
    files = sorted(
        f for f in paper_dir.iterdir()
        if f.is_file() and not f.name.startswith(".")
    )

    if not files:
        st.info("No files found in paper directory.")
        return

    # Group by type
    csvs = [f for f in files if f.suffix.lower() in (".csv", ".tsv")]
    excels = [f for f in files if f.suffix.lower() in (".xlsx", ".xls")]
    pdfs = [f for f in files if f.suffix.lower() == ".pdf"]
    yamls = [f for f in files if f.suffix.lower() in (".yaml", ".yml")]
    others = [f for f in files if f not in csvs + excels + pdfs + yamls]

    # Paperconfig
    if paperconfig_path.exists():
        with st.expander("paperconfig.yaml", expanded=False):
            with open(paperconfig_path) as fh:
                st.code(fh.read(), language="yaml")

    # CSVs
    for csv_file in csvs:
        with st.expander(f"CSV: {csv_file.name}", expanded=False):
            try:
                sep = "\t" if csv_file.suffix.lower() == ".tsv" else ","
                df = pd.read_csv(csv_file, sep=sep, nrows=200)
                st.dataframe(df, use_container_width=True)
                if len(df) == 200:
                    st.caption("Showing first 200 rows")
            except Exception as e:
                st.error(f"Could not read {csv_file.name}: {e}")

    # Excel files
    for xls_file in excels:
        with st.expander(f"Excel: {xls_file.name}", expanded=False):
            try:
                xf = pd.ExcelFile(xls_file)
                sheet_names = xf.sheet_names
                if len(sheet_names) > 1:
                    selected_sheet = st.selectbox(
                        "Sheet",
                        sheet_names,
                        key=f"sheet_{xls_file.name}",
                    )
                else:
                    selected_sheet = sheet_names[0]
                df = pd.read_excel(xf, sheet_name=selected_sheet, nrows=200)
                st.dataframe(df, use_container_width=True)
                if len(df) == 200:
                    st.caption("Showing first 200 rows")
            except Exception as e:
                st.error(f"Could not read {xls_file.name}: {e}")

    # PDFs
    for pdf_file in pdfs:
        col1, col2 = st.columns([3, 1])
        with col1:
            st.markdown(f"PDF: **{pdf_file.name}**")
        with col2:
            if st.button(f"Open", key=f"open_{pdf_file.name}"):
                try:
                    webbrowser.open(f"file://{pdf_file.resolve()}")
                except Exception:
                    subprocess.Popen(["xdg-open", str(pdf_file)])

    # Other YAML files (not paperconfig)
    for yf in yamls:
        if yf.resolve() != paperconfig_path.resolve():
            with st.expander(f"YAML: {yf.name}", expanded=False):
                with open(yf) as fh:
                    st.code(fh.read(), language="yaml")

    # Other files
    if others:
        st.markdown("**Other files:**")
        for f in others:
            col1, col2 = st.columns([3, 1])
            with col1:
                st.markdown(f"`{f.name}` ({f.suffix})")
            with col2:
                if st.button("Open", key=f"open_other_{f.name}"):
                    try:
                        subprocess.Popen(["xdg-open", str(f)])
                    except Exception as e:
                        st.error(f"Could not open: {e}")
