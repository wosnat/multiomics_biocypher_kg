#!/usr/bin/env bash
# Orchestration: download + preprocess genome annotation data.
#
# Step 0 — Download genome data (NCBI GFF, Cyanorak GFF/GBK, UniProt, gene_mapping.csv)
#           calls: multiomics_kg/download/download_genome_data.py --steps 1 2 3 5
#           (eggNOG step 4 is intentionally skipped — run /eggnog-run skill separately)
#           Use --skip-cyanorak to skip Cyanorak downloads (step 2) when server is slow
# Step 1 — Build per-taxid protein annotation tables (protein_annotations.json)
#           calls: multiomics_kg/download/build_protein_annotations.py
#           Requires step 0 (UniProt data must be cached first)
# Step 2 — Build per-strain gene annotation tables (gene_annotations_merged.json)
#           calls: multiomics_kg/download/build_gene_annotations.py
#           Requires step 1 (protein_annotations.json used as UniProt source)
# Step 3 — Build per-strain gene ID mapping (gene_id_mapping.json + gene_mapping_supp.csv)
#           calls: multiomics_kg/download/build_gene_id_mapping.py
#           Reads paperconfig.yaml files to add paper-derived alt-IDs (JGI IDs, probesets, etc.)
#           Requires step 2 (gene_annotations_merged.json used as base)
#           Must run before step 4
# Step 4 — Resolve paper CSV gene IDs to locus tags
#           calls: multiomics_kg/download/resolve_paper_ids.py
#           Reads all paperconfig CSV tables, resolves name_col → locus_tag via gene_id_mapping.json
#           Writes <stem>_resolved.csv alongside each source CSV
#           Requires step 3 (gene_id_mapping.json must exist for paper-derived IDs)
#           Must run before create_knowledge_graph.py so omics_adapter uses pre-resolved files
# Step 5 — Extract eggNOG OG descriptions (og_descriptions.json)
#           calls: multiomics_kg/download/build_og_descriptions.py
#           Queries local eggnog.db for ortholog group descriptions; writes lightweight
#           cache at cache/data/eggnog/og_descriptions.json (avoids 39GB DB in Docker)
#           Requires step 2 (gene_annotations_merged.json for OG ID list)
#
# Step 6 — Build pruned KEGG data cache (kegg_data.json) + TCDB hierarchy with metabolism enrichment
#           calls: multiomics_kg/download/build_kegg_metabolism_xrefs.py
#           Walks every strain's gene_annotations_merged.json to identify gene-reachable
#           {KOs, reactions, compounds, pathways}; prunes raw KEGG to that subset; enriches
#           reactions/compounds with MNX xrefs; writes a single cache/data/kegg/kegg_data.json (~3-4 MB, indented).
#           Also downloads the 3 TCDB reference TSVs and writes cache/data/tcdb/tcdb_hierarchy.json.
#           Requires step 2 + scripts/refresh_mnx.sh
#           (MNX resolver — heavy one-time build, rerun only when MNX releases)
#           Step 6 also walks all paperconfig metabolite_assays_table entries
#           (Phase 2 metabolomics): harvests paper-measured metabolites, unions
#           "metabolomics" into evidence_sources, writes
#           cache/data/metabolomics/metabolite_id_mapping.json (consumed by step 7).
#
# Step 7 — Resolve paper metabolite names → primary IDs
#           calls: multiomics_kg/download/resolve_paper_metabolites.py
#           For each metabolite_assays_table entry, opens the source CSV and
#           writes <stem>_resolved.csv with metabolite_id + resolution_method
#           columns; mirrors step 4 for genes. Reads only metabolite_id_mapping.json
#           (no MNX resolver hit — step 6 owns that).
#           Requires step 6.
# Step 8 — Resolve paper "discusses" topic mentions → Gene / KEGG-pathway ids
#           calls: multiomics_kg/download/resolve_paper_topics.py
#           For each paper with a publication_topics/topics.json (written by the
#           /extract-discussed-topics skill — the LLM extraction step is NOT in
#           prepare_data), resolves gene + pathway mentions and writes
#           topics_resolved.json + resolution_report.txt. Deterministic; mirrors
#           steps 4/7. Papers without a topics.json are skipped.
#           Requires step 3 (gene_id_mapping.json) + step 6 (kegg_data.json).
# Step 9 — Build the central reference caches (interpro_reference.json + ncbifam_reference.json)
#           calls: multiomics_kg/download/build_interpro_reference.py
#                  multiomics_kg/download/build_ncbifam_reference.py
#           Downloads InterPro current_release entry.list + ParentChildTreeFile.txt
#           (+ interpro2go + interpro.xml.gz) and NCBI's hmm_PGAP.tsv, writing
#           cache/data/interpro/interpro_reference.json ({IPR: name/type/parent/level/...})
#           and cache/data/ncbifam/ncbifam_reference.json ({unversioned_acc: name/family_type/...}),
#           consumed by interpro_adapter / ncbifam_adapter for node names/metadata + hierarchy.
#           Global reference downloads (no per-strain data needed), but the step-2 merge
#           (build_gene_annotations.py) consumes them (interpro_reference.json already;
#           ncbifam_reference.json per spec §4), so step 9 now runs BEFORE step 2 in the
#           default order (see STEPS below). --refetch-raw re-pulls the FTP/NCBI files
#           (only on an InterPro/NCBIfam release).
# Step 10 — Build MEROPS reference cache (merops_reference.json + pfam_bridge + cleavage)
#           calls: multiomics_kg/download/build_merops_reference.py
#           Writes cache/data/merops/merops_reference.json, consumed by merops_adapter
#           for node names/clan descriptions/typing at KG-build time. No ordering
#           constraint (step 2 merge does not consume it); runs after step 2 by default
#           but can run independently via --steps 10.
#
# Logs: logs/prepare_data_step0.log … logs/prepare_data_step10.log
#       Monitor with: tail -f logs/prepare_data_step0.log
#
# Usage:
#   ./scripts/prepare_data.sh
#   ./scripts/prepare_data.sh --force
#   ./scripts/prepare_data.sh --skip-cyanorak --force   # skip slow Cyanorak server
#   ./scripts/prepare_data.sh --strains MED4 MIT9313
#   ./scripts/prepare_data.sh --steps 0 1 2 3
#   ./scripts/prepare_data.sh --steps 0 --force
#   ./scripts/prepare_data.sh --steps 1 --force         # rebuild protein_annotations.json only
#   ./scripts/prepare_data.sh --steps 2 --strains MED4 --force
#   ./scripts/prepare_data.sh --steps 3 --strains MIT9301 --force  # rebuild gene_id_mapping only
#   ./scripts/prepare_data.sh --steps 6 7 --force                  # rebuild KEGG/TCDB caches from cached raw inputs (fast iteration)
#   ./scripts/prepare_data.sh --steps 6 --refetch-raw              # also re-pull raw KEGG REST + TCDB TSVs (slow; only on upstream releases)
#   ./scripts/prepare_data.sh --steps 10                           # build MEROPS reference only
#   ./scripts/prepare_data.sh --rebuild                            # all derived steps (1-10, 9 before 2) with --force; step-0 downloads excluded

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
LOG_DIR="$PROJECT_ROOT/logs"

mkdir -p "$LOG_DIR"

# ── parse args ────────────────────────────────────────────────────────────────

FORCE=""
REFETCH_RAW=""
# 9 = central references (interpro + ncbifam); runs BEFORE the step-2 merge which consumes them (spec 2026-08-17 §4)
# 10 = MEROPS reference; no ordering constraint with other steps
STEPS="0 9 1 2 3 4 5 6 7 8 10"
STRAINS=()
SKIP_CYANORAK=0
USER_STEPS=0
REBUILD=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --force)          FORCE="--force"; shift ;;
        --refetch-raw)    REFETCH_RAW="--refetch-raw"; shift ;;
        --skip-cyanorak)  SKIP_CYANORAK=1; shift ;;
        --rebuild)        REBUILD=1; shift ;;
        --steps)
            USER_STEPS=1
            STEPS=""
            shift
            while [[ $# -gt 0 && "$1" != --* ]]; do
                STEPS="$STEPS $1"; shift
            done
            ;;
        --strains)
            shift
            while [[ $# -gt 0 && "$1" != --* ]]; do
                STRAINS+=("$1"); shift
            done
            ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

if [[ $REBUILD -eq 1 ]]; then
    if [[ $USER_STEPS -eq 1 ]]; then
        echo "--rebuild and --steps are mutually exclusive" >&2; exit 1
    fi
    # All derived steps in DEPENDENCY order (9 before 2 — step 2 consumes the
    # step-9 reference caches). Step 0 (raw downloads) deliberately excluded.
    STEPS="9 1 2 3 4 5 6 7 8 10"
    FORCE="--force"
fi

STRAINS_ARG=""
if [[ ${#STRAINS[@]} -gt 0 ]]; then
    STRAINS_ARG="--strains ${STRAINS[*]}"
fi

# ── helpers ───────────────────────────────────────────────────────────────────

run_step() {
    local step_num="$1"
    local label="$2"
    local log="$3"
    shift 3
    local cmd=("$@")

    echo ""
    echo "══════════════════════════════════════════════════════════════════════"
    echo "  Step $step_num — $label"
    echo "  Log: $log"
    echo "  Cmd: ${cmd[*]}"
    echo "══════════════════════════════════════════════════════════════════════"
    echo ""

    # RUN_STEP_APPEND=1 appends to an existing log instead of overwriting it
    # (used when a single step number runs more than one command, e.g. step 9)
    local tee_args=("$log")
    if [[ "${RUN_STEP_APPEND:-0}" -eq 1 ]]; then
        tee_args=("-a" "$log")
    fi

    # tee: show output live AND write to log
    if ! "${cmd[@]}" 2>&1 | tee "${tee_args[@]}"; then
        echo ""
        echo "ERROR: Step $step_num failed — check $log" >&2
        exit 1
    fi
    echo ""
    echo "  Step $step_num complete. Full log: $log"
}

# ── steps ─────────────────────────────────────────────────────────────────────

cd "$PROJECT_ROOT"
export PYTHONPATH="$PROJECT_ROOT${PYTHONPATH:+:$PYTHONPATH}"

echo "prepare_data.sh: steps=[${STEPS}]${STRAINS_ARG:+ strains=[${STRAINS[*]}]}${FORCE:+ (force)}${REFETCH_RAW:+ (refetch-raw)}${SKIP_CYANORAK:+ (skip-cyanorak)}"
echo "(step 1 = protein annotations, step 2 = gene annotations, step 3 = gene ID mapping, step 4 = resolve paper CSVs, step 5 = OG descriptions, step 6 = pruned KEGG + TCDB hierarchy caches, step 7 = resolve paper metabolite names, step 8 = resolve paper discuss-topics, step 9 = central references: InterPro + NCBIfam reference caches, step 10 = MEROPS reference cache)"
echo "Project root: $PROJECT_ROOT"
echo "Logs dir:     $LOG_DIR"

for step in $STEPS; do
    case "$step" in
        0)
            if [[ $SKIP_CYANORAK -eq 1 ]]; then
                DOWNLOAD_SUBSTEPS="1 3 5"
                STEP0_LABEL="Download genome data (NCBI + UniProt + gene_mapping; Cyanorak SKIPPED)"
            else
                DOWNLOAD_SUBSTEPS="1 2 3 5"
                STEP0_LABEL="Download genome data (NCBI + Cyanorak + UniProt + gene_mapping)"
            fi
            run_step 0 \
                "$STEP0_LABEL" \
                "$LOG_DIR/prepare_data_step0.log" \
                uv run python multiomics_kg/download/download_genome_data.py \
                    --steps $DOWNLOAD_SUBSTEPS \
                    $STRAINS_ARG \
                    $FORCE
            ;;
        1)
            run_step 1 \
                "Build protein annotation tables (protein_annotations.json)" \
                "$LOG_DIR/prepare_data_step1.log" \
                uv run python multiomics_kg/download/build_protein_annotations.py \
                    $STRAINS_ARG \
                    $FORCE
            ;;
        2)
            run_step 2 \
                "Build gene annotation tables (gene_annotations_merged.json)" \
                "$LOG_DIR/prepare_data_step2.log" \
                uv run python multiomics_kg/download/build_gene_annotations.py \
                    $STRAINS_ARG \
                    $FORCE
            ;;
        3)
            run_step 3 \
                "Build gene ID mapping (gene_id_mapping.json + gene_mapping_supp.csv)" \
                "$LOG_DIR/prepare_data_step3.log" \
                uv run python -m multiomics_kg.download.build_gene_id_mapping \
                    $STRAINS_ARG \
                    $FORCE
            ;;
        4)
            run_step 4 \
                "Resolve paper CSV gene IDs → locus tags (_resolved.csv)" \
                "$LOG_DIR/prepare_data_step4.log" \
                uv run python -m multiomics_kg.download.resolve_paper_ids \
                    $FORCE
            ;;
        5)
            run_step 5 \
                "Extract eggNOG OG descriptions (og_descriptions.json)" \
                "$LOG_DIR/prepare_data_step5.log" \
                uv run python -m multiomics_kg.download.build_og_descriptions \
                    $FORCE
            ;;
        6)
            run_step 6 \
                "Build pruned KEGG + TCDB hierarchy caches (kegg_data.json + tcdb_hierarchy.json)" \
                "$LOG_DIR/prepare_data_step6.log" \
                uv run python -m multiomics_kg.download.build_kegg_metabolism_xrefs \
                    $FORCE \
                    $REFETCH_RAW
            ;;
        7)
            run_step 7 \
                "Resolve paper metabolite names → primary IDs (_resolved.csv)" \
                "$LOG_DIR/prepare_data_step7.log" \
                uv run python -m multiomics_kg.download.resolve_paper_metabolites \
                    $FORCE
            ;;
        8)
            run_step 8 \
                "Resolve paper discuss-topic mentions → Gene/KEGG ids (topics_resolved.json)" \
                "$LOG_DIR/prepare_data_step8.log" \
                uv run python -m multiomics_kg.download.resolve_paper_topics \
                    $FORCE
            ;;
        9)
            run_step 9 \
                "Build InterPro reference cache (interpro_reference.json)" \
                "$LOG_DIR/prepare_data_step9.log" \
                uv run python -m multiomics_kg.download.build_interpro_reference \
                    $FORCE \
                    $REFETCH_RAW

            RUN_STEP_APPEND=1 run_step 9 \
                "Build NCBIfam reference cache (ncbifam_reference.json)" \
                "$LOG_DIR/prepare_data_step9.log" \
                uv run python -m multiomics_kg.download.build_ncbifam_reference \
                    $FORCE \
                    $REFETCH_RAW
            ;;
        10)
            run_step 10 \
                "Build MEROPS reference cache (merops_reference.json + pfam_bridge + cleavage)" \
                "$LOG_DIR/prepare_data_step10.log" \
                uv run python -m multiomics_kg.download.build_merops_reference \
                    $FORCE \
                    $REFETCH_RAW
            ;;
        *)
            echo "Unknown step: $step (valid: 0 1 2 3 4 5 6 7 8 9 10)" >&2
            exit 1
            ;;
    esac
done

echo ""
echo "══════════════════════════════════════════════════════════════════════"
echo "  All steps complete."
echo "══════════════════════════════════════════════════════════════════════"
