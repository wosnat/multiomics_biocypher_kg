"""Streamlit app for reviewing cluster extraction results.

Launch with:
    uv run streamlit run multiomics_kg/review/cluster_review_app.py
"""

import csv
import io
import sys
from pathlib import Path

# Ensure project root is on sys.path for module resolution
_project_root = str(Path(__file__).resolve().parent.parent.parent)
if _project_root not in sys.path:
    sys.path.insert(0, _project_root)

import streamlit as st
import yaml

from multiomics_kg.extraction.cluster.run_manager import RunManager
from multiomics_kg.review.review_components import (
    VERDICT_COLORS,
    render_diff_view,
    render_merge_view,
    render_quotes,
    render_review_controls,
    render_source_access,
    render_verdict,
    _verdict_badge,
)
from multiomics_kg.review.review_data import (
    compute_entry_status_color,
    export_issue_report,
    load_entry_summary,
    scan_papers_with_clusters,
)

# Project root: multiomics_kg/review/cluster_review_app.py -> project root
PROJECT_ROOT = Path(__file__).parent.parent.parent

# Paperconfig list files
PAPERCONFIG_LIST_FILES = [
    PROJECT_ROOT / "data" / "Prochlorococcus" / "papers_and_supp" / "paperconfig_files.txt",
    PROJECT_ROOT / "data" / "Synechococcus" / "papers_and_supp" / "paperconfig_files.txt",
]

# Status color mapping for sidebar
STATUS_COLOR_MAP = {
    "green": ("🟢", "All approved (current run)"),
    "light_green": ("🟡", "All approved (carried forward)"),
    "yellow": ("🟠", "Has stale reviews"),
    "red": ("🔴", "Unreviewed or rejected"),
}

# Review status icons
REVIEW_STATUS_ICONS = {
    "approve": "✅",
    "edit": "✏️",
    "reject": "❌",
    "flag-issue": "⚠️",
    "stale": "🔄",
    "unreviewed": "⬜",
}


def load_paperconfig_paths() -> list[Path]:
    """Load all paperconfig paths from the list files."""
    paths = []
    for list_file in PAPERCONFIG_LIST_FILES:
        if not list_file.exists():
            continue
        with open(list_file) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                p = PROJECT_ROOT / line
                if p.exists():
                    paths.append(p)
    return paths


@st.cache_data(ttl=30)
def cached_scan_papers(paperconfig_paths_str: list[str]) -> list[dict]:
    """Cached wrapper for scan_papers_with_clusters."""
    paths = [Path(p) for p in paperconfig_paths_str]
    return scan_papers_with_clusters(paths)


def get_paper_color(entries: list[dict]) -> str:
    """Compute worst-case color across all entries for a paper."""
    priority = {"red": 0, "yellow": 1, "light_green": 2, "green": 3}
    worst = "green"
    for e in entries:
        color = compute_entry_status_color(e["paper_dir"], e["entry_key"])
        if priority.get(color, 4) < priority.get(worst, 4):
            worst = color
    return worst


def main():
    st.set_page_config(page_title="Cluster Extraction Review", layout="wide")
    st.title("Cluster Extraction Review")

    # Load data
    pc_paths = load_paperconfig_paths()
    if not pc_paths:
        st.error("No paperconfig files found. Check paperconfig_files.txt locations.")
        return

    entries = cached_scan_papers([str(p) for p in pc_paths])
    if not entries:
        st.info("No gene_clusters entries found in any paperconfig.")
        return

    # Group entries by paper
    papers = {}
    for e in entries:
        papers.setdefault(e["paper_name"], []).append(e)

    # --- Sidebar ---
    with st.sidebar:
        st.header("Navigation")

        # Paper selector with color coding
        paper_names = sorted(papers.keys())
        paper_labels = []
        for pn in paper_names:
            color = get_paper_color(papers[pn])
            icon = STATUS_COLOR_MAP.get(color, ("⬜", ""))[0]
            paper_labels.append(f"{icon} {pn}")

        selected_paper_idx = st.selectbox(
            "Paper",
            range(len(paper_names)),
            format_func=lambda i: paper_labels[i],
            key="paper_select",
        )
        selected_paper = paper_names[selected_paper_idx]
        paper_entries = papers[selected_paper]

        # Entry selector with color coding
        entry_labels = []
        for e in paper_entries:
            color = compute_entry_status_color(e["paper_dir"], e["entry_key"])
            icon = STATUS_COLOR_MAP.get(color, ("⬜", ""))[0]
            entry_labels.append(f"{icon} {e['entry_key']} ({e['organism']})")

        selected_entry_idx = st.selectbox(
            "Analysis entry",
            range(len(paper_entries)),
            format_func=lambda i: entry_labels[i],
            key="entry_select",
        )
        selected_entry = paper_entries[selected_entry_idx]

        st.markdown("---")

        # Summary stats
        summary = load_entry_summary(selected_entry["paper_dir"], selected_entry["entry_key"])
        if summary["has_run"]:
            st.metric("Clusters", summary["cluster_count"])

            # Verdict breakdown
            st.markdown("**Verdicts:**")
            for v, count in sorted(summary["verdicts"].items()):
                bg, fg = VERDICT_COLORS.get(v, ("#9e9e9e", "white"))
                st.markdown(f"&nbsp;&nbsp;{v}: **{count}**")

            # Review status breakdown
            st.markdown("**Review status:**")
            for rs, count in sorted(summary["review_statuses"].items()):
                icon = REVIEW_STATUS_ICONS.get(rs, "")
                st.markdown(f"&nbsp;&nbsp;{icon} {rs}: **{count}**")
        else:
            st.warning("No extraction run found for this entry.")

        st.markdown("---")

        # Filters
        st.markdown("**Filters**")
        all_verdicts = sorted(set(summary.get("verdicts", {}).keys()))
        all_statuses = sorted(set(summary.get("review_statuses", {}).keys()))

        filter_verdicts = st.multiselect(
            "Filter by verdict",
            all_verdicts,
            default=[],
            key="filter_verdict",
        )
        filter_statuses = st.multiselect(
            "Filter by review status",
            all_statuses,
            default=[],
            key="filter_status",
        )

        st.markdown("---")

        # Export issue report
        if st.button("Export issue report (CSV)"):
            all_entries_for_export = [(e["paper_dir"], e["entry_key"]) for e in entries]
            report_rows = export_issue_report(all_entries_for_export)
            if report_rows:
                output = io.StringIO()
                writer = csv.DictWriter(output, fieldnames=report_rows[0].keys())
                writer.writeheader()
                writer.writerows(report_rows)
                st.download_button(
                    "Download CSV",
                    output.getvalue(),
                    file_name="cluster_issue_report.csv",
                    mime="text/csv",
                )
            else:
                st.info("No issues to report.")

        st.markdown("---")

        # Re-run Extraction section
        st.header("Re-run Extraction")
        force_mode = st.selectbox(
            "Force mode",
            ["Default (cached pre-processing)", "Force LLM", "Force all"],
            key="rerun_force_mode",
        )
        if st.button("Re-run this entry", key="sidebar_rerun"):
            from multiomics_kg.extraction.cluster.pipeline import run_pipeline

            pc_path = Path(str(selected_entry["paperconfig_path"]))
            ek = selected_entry["entry_key"]
            force_llm = force_mode == "Force LLM"
            force_all = force_mode == "Force all"
            with st.spinner(f"Re-running extraction for {ek}..."):
                run_pipeline(
                    pc_path,
                    table_key=ek,
                    force_llm=force_llm,
                    force_all=force_all,
                )
            st.success(f"Re-run complete for {ek}.")
            # Clear cached data and rerun
            cached_scan_papers.clear()
            st.rerun()

    # --- Main area ---
    if not summary.get("has_run"):
        st.warning("No extraction run available for this entry. Run the extraction pipeline first.")
        return

    # Load stage data
    cache_dir = selected_entry["paper_dir"] / ".extraction_cache"
    rm = RunManager(cache_dir, selected_entry["entry_key"])
    run_dir = rm.get_current_run()

    # Store in session state for re-run buttons in review_components
    st.session_state.paperconfig_path = str(selected_entry["paperconfig_path"])
    st.session_state.entry_key = selected_entry["entry_key"]

    if run_dir is None:
        st.error("Could not find current run directory.")
        return

    stage1 = rm.read_stage(run_dir, 1)
    stage2 = rm.read_stage(run_dir, 2)
    stage3 = rm.read_stage(run_dir, 3)
    stage4 = rm.read_stage(run_dir, 4)

    # Tabs
    tab_clusters, tab_sources = st.tabs(["Clusters", "Sources"])

    with tab_clusters:
        cluster_keys = sorted(stage2.keys(), key=lambda k: (int(k) if k.isdigit() else k))

        # Apply filters
        filtered_keys = []
        for ck in cluster_keys:
            verdict = stage3.get(ck, {}).get("verdict", "none")
            review_status = stage4.get(ck, {}).get("status", "unreviewed")

            if filter_verdicts and verdict not in filter_verdicts:
                continue
            if filter_statuses and review_status not in filter_statuses:
                continue
            filtered_keys.append(ck)

        if not filtered_keys:
            st.info("No clusters match the current filters.")
        else:
            st.markdown(f"Showing **{len(filtered_keys)}** of {len(cluster_keys)} clusters")

        for ck in filtered_keys:
            s1_cluster = stage1.get(ck, {})
            s2_cluster = stage2.get(ck, {})
            s3_cluster = stage3.get(ck, {})
            s4_cluster = stage4.get(ck, {})

            review_status = s4_cluster.get("status", "unreviewed")
            verdict = s3_cluster.get("verdict", "none")
            icon = REVIEW_STATUS_ICONS.get(review_status, "")
            cluster_id = s2_cluster.get("id", ck)

            # Auto-expand non-approved clusters
            expanded = review_status not in ("approve",)

            with st.expander(
                f"{icon} Cluster {ck}: {cluster_id} (verdict: {verdict}, review: {review_status})",
                expanded=expanded,
            ):
                # Merge view (includes synthesis result)
                render_merge_view(ck, s1_cluster, s2_cluster)

                # Validation verdict
                st.subheader("Validation")
                render_verdict(ck, s3_cluster)

                # Review controls
                st.subheader("Review")
                render_review_controls(ck, s2_cluster, s4_cluster, run_dir, rm, stage1)

                # Diff view (only if 2+ runs exist)
                if len(rm.list_runs()) >= 2:
                    with st.expander("Compare with previous run", expanded=False):
                        render_diff_view(rm, ck)

    with tab_sources:
        render_source_access(
            selected_entry["paper_dir"],
            selected_entry["paperconfig_path"],
        )


if __name__ == "__main__":
    main()
