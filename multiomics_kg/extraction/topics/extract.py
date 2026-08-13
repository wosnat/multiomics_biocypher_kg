"""Publication "discusses" topic extraction via the OpenAI Responses API.

Stage 1 of the discuss-edges pipeline (extract -> resolve -> adapter). Feeds the
FULL paper PDF to the model and extracts the genes + KEGG pathways the paper
substantively discusses in its narrative, with verbatim identifiers and a strain
attribution so Stage 2 can resolve each mention to a strain-specific Gene node.

Usage:
    uv run python -m multiomics_kg.extraction.topics.extract              # extract all
    uv run python -m multiomics_kg.extraction.topics.extract --paper "Zinser 2009"
    uv run python -m multiomics_kg.extraction.topics.extract --force
    uv run python -m multiomics_kg.extraction.topics.extract --dry-run
    uv run python -m multiomics_kg.extraction.topics.extract --report     # no API calls

Modelled on multiomics_kg/extraction/cluster/extract.py and sharing
extraction/pdf_utils.py (upload_pdf) + the Responses API + pydantic pattern.
"""
import argparse
import logging
import os
import tempfile
import time
from datetime import datetime
from pathlib import Path
from typing import Literal

from pydantic import BaseModel

from multiomics_kg.extraction.pdf_utils import (
    count_pages,
    find_references_page,
    page_chunks,
    upload_pdf,
    write_page_range_pdf,
)
from multiomics_kg.extraction.topics.extraction_utils import (
    find_all_papers,
    has_topics,
    load_topics,
    save_topics,
)

logger = logging.getLogger(__name__)


# ── Pydantic schemas ──


class GeneMention(BaseModel):
    surface_form: str          # exact string the paper uses, e.g. "NtcA" or "PMM0552"
    gene_name: str             # standard gene symbol if any, else ""
    identifiers: list[str]     # verbatim locus tags / ORF / gene / protein ids printed
    strain: str                # one of the paper strains, or "all" / "unspecified"
    prominence: Literal["central", "peripheral"]
    evidence: str              # short direct quote + location supporting the mention


class PathwayMention(BaseModel):
    surface_form: str          # how the paper names the pathway
    kegg_pathway_id: str       # e.g. "ko00910" / "map00910" / "00910"; "" if unknown
    prominence: Literal["central", "peripheral"]
    evidence: str


class TopicsExtraction(BaseModel):
    genes: list[GeneMention]
    pathways: list[PathwayMention]
    self_assessment: Literal["high", "medium", "low"]
    assessment_notes: str
    # Triage signal: how many identifiers the paper PRINTS were left uncaptured
    # (e.g. abbreviated per-strain locus numbers in tables). Lets us later decide
    # whether better ID capture is worth a prompt refinement or is paper-specific.
    uncaptured_identifiers: Literal["none", "some", "many"]
    uncaptured_identifiers_note: str


# ── Prompt layer ──


RULES = """\
## What to extract

Extract the GENES and KEGG PATHWAYS this paper substantively DISCUSSES in its
narrative text (abstract, introduction, results, discussion) — the genes and
pathways the authors name and reason about. This is NOT the full supplementary
table: include a gene only if the paper's prose actually mentions it.

## Rules

- ONLY include genes/pathways the paper EXPLICITLY names. Do NOT infer from your
  own knowledge, and do NOT invent identifiers the paper does not print.
- identifiers (genes): capture EVERY explicit identifier the paper prints for the
  gene — locus tags (e.g. PMM0552, PMT0001, P9301_05911, SYNW1234, A9601_..., WP_...),
  ORF names, gene/protein IDs. This is the single most useful field: list them
  verbatim, exactly as printed. Empty list only when the paper gives no ID at all.
- gene_name: the standard gene symbol (e.g. "ntcA", "psbA") if the paper uses one,
  else "".
- strain attribution: set `strain` to the ONE strain from the candidate list below
  that the paper attributes this gene to. If the paper discusses it generically across
  all its strains (e.g. a conserved regulator), use "all". If you genuinely cannot
  tell, use "unspecified". A locus-tag prefix usually tells you the strain.
- prominence: "central" if the gene/pathway is a focus of the paper's findings or
  discussion; "peripheral" if mentioned only in passing.
- pathways: use KEGG pathway language. Put the KEGG pathway id in `kegg_pathway_id`
  when you know it (ko/map number); otherwise "" and rely on `surface_form`.
- evidence: a SHORT direct quote from the paper (<200 chars) plus a location hint,
  e.g. "‘ntcA was strongly induced under N stress’ (Results, p4)". REQUIRED.
- Do NOT use the characters '|' or single-quote-as-delimiter in any field.
- self_assessment: your confidence in the overall extraction (high/medium/low);
  assessment_notes: what you were unsure about (e.g. ambiguous strain attribution).
- uncaptured_identifiers: a TRIAGE signal — after extracting, judge honestly how
  many gene identifiers that the paper PRINTS were left OUT of every gene's
  `identifiers` list. This explicitly includes abbreviated / partial locus numbers
  in tables (e.g. "ntcA (0246*)" → the bare "0246"), per-strain id columns, or
  supplementary id lists you chose not to attach. Use "none", "some", or "many".
  This is NOT a request to capture them now — just report the gap.
- uncaptured_identifiers_note: one line on WHERE those uncaptured ids are
  (e.g. "Table 2 prints abbreviated per-strain locus numbers"). "" if none.
"""


def build_prompt(strains: list[str]) -> str:
    strain_block = (
        "\n".join(f"- {s}" for s in strains)
        if strains else "- (none declared in paperconfig; attribute as best you can)"
    )
    return (
        "You are extracting, from a scientific paper, the genes and KEGG pathways "
        "it discusses in its narrative, to build a literature index for a "
        "knowledge graph.\n\n"
        "## Candidate strains studied by this paper\n"
        f"{strain_block}\n\n"
        f"{RULES}"
    )


# ── LLM layer ──


def extract_topics(client, file_ids, strains, model="gpt-5.4", flex=False):
    """Run topic extraction for one paper. Returns (parsed, usage_dict)."""
    dev_msg = build_prompt(strains)
    file_inputs = [{"type": "input_file", "file_id": fid} for fid in file_ids]

    kwargs = dict(
        model=model, temperature=0,
        input=[
            {"role": "developer", "content": dev_msg},
            {"role": "user", "content": file_inputs + [
                {"type": "input_text",
                 "text": "Extract all discussed genes and KEGG pathways."},
            ]},
        ],
        text_format=TopicsExtraction,
    )
    if flex:
        kwargs["service_tier"] = "flex"

    t0 = time.time()
    resp = client.responses.parse(**kwargs)
    elapsed = time.time() - t0
    usage = {
        "input_tokens": resp.usage.input_tokens,
        "output_tokens": resp.usage.output_tokens,
        "duration_sec": elapsed,
    }
    return resp.output[0].content[0].parsed, usage


# ── Chunked extraction (split large PDFs into <=N-page chunks, skip refs, merge) ──

_CONF_RANK = {"low": 0, "medium": 1, "high": 2}
_UNC_RANK = {"none": 0, "some": 1, "many": 2}


def merge_extractions(parsed_list: list[dict]) -> dict:
    """Merge per-chunk extraction dicts into one paper-level result.

    Genes dedup by (gene_name|surface, strain): union identifiers, prominence
    central-wins. Pathways dedup by kegg id (or surface). Confidence = lowest
    chunk's; uncaptured_identifiers = highest chunk's.
    """
    genes: dict[tuple[str, str], dict] = {}
    for p in parsed_list:
        for g in p.get("genes", []):
            name = (g.get("gene_name") or g.get("surface_form") or "").strip().lower()
            strain = (g.get("strain") or "").strip().lower()
            key = (name, strain)
            ids = list(dict.fromkeys(g.get("identifiers") or []))
            if key not in genes:
                cur = dict(g)
                cur["identifiers"] = ids
                genes[key] = cur
            else:
                cur = genes[key]
                cur["identifiers"] = list(dict.fromkeys((cur.get("identifiers") or []) + ids))
                if g.get("prominence") == "central":
                    cur["prominence"] = "central"
                if not cur.get("gene_name") and g.get("gene_name"):
                    cur["gene_name"] = g["gene_name"]

    pathways: dict[str, dict] = {}
    for p in parsed_list:
        for pw in p.get("pathways", []):
            key = (pw.get("kegg_pathway_id") or "").strip().lower() \
                or "name:" + (pw.get("surface_form") or "").strip().lower()
            if key not in pathways:
                pathways[key] = dict(pw)
            else:
                cur = pathways[key]
                if pw.get("prominence") == "central":
                    cur["prominence"] = "central"
                if not cur.get("kegg_pathway_id") and pw.get("kegg_pathway_id"):
                    cur["kegg_pathway_id"] = pw["kegg_pathway_id"]

    conf = min((_CONF_RANK.get(p.get("self_assessment", "medium"), 1) for p in parsed_list), default=1)
    unc = max((_UNC_RANK.get(p.get("uncaptured_identifiers", "none"), 0) for p in parsed_list), default=0)
    notes = "; ".join(dict.fromkeys(
        n for p in parsed_list if (n := (p.get("uncaptured_identifiers_note") or "").strip())
    ))
    return {
        "genes": list(genes.values()),
        "pathways": list(pathways.values()),
        "self_assessment": {v: k for k, v in _CONF_RANK.items()}[conf],
        "uncaptured_identifiers": {v: k for k, v in _UNC_RANK.items()}[unc],
        "uncaptured_identifiers_note": notes,
    }


def extract_topics_chunked(client, pdf_path, strains, model="gpt-5.4",
                           flex=False, chunk_pages=15):
    """Extract by splitting the PDF into <=chunk_pages-page chunks (reference pages
    skipped), one request per chunk, then merging. Returns (merged, usage, n_chunks).
    """
    total = count_pages(pdf_path)
    if total <= 0:
        raise RuntimeError(f"no readable pages in {pdf_path}")
    refs = find_references_page(pdf_path)
    # Keep the page that holds the References heading (it has content above it);
    # drop everything after. If no heading found, keep all pages.
    last_page = refs if refs is not None else total - 1
    chunks = page_chunks(last_page, chunk_pages)

    parsed_list: list[dict] = []
    usage = {"input_tokens": 0, "output_tokens": 0, "duration_sec": 0.0}
    with tempfile.TemporaryDirectory() as td:
        for idx, (start, end) in enumerate(chunks):
            dest = Path(td) / f"chunk_{idx}.pdf"
            if write_page_range_pdf(pdf_path, start, end, dest) is None:
                continue
            fid = upload_pdf(client, dest)
            try:
                parsed, u = extract_topics(client, [fid], strains, model=model, flex=flex)
                parsed_list.append(parsed.model_dump())
                usage["input_tokens"] += u["input_tokens"]
                usage["output_tokens"] += u["output_tokens"]
                usage["duration_sec"] += u["duration_sec"]
            except Exception as e:
                logger.warning("    chunk %d (pp %d-%d) failed: %s", idx, start, end, e)
            finally:
                try:
                    client.files.delete(fid)
                except Exception:
                    pass
                time.sleep(2)

    if not parsed_list:
        raise RuntimeError("all chunks failed")
    return merge_extractions(parsed_list), usage, len(chunks)


# ── Report layer ──


def generate_report(papers) -> str:
    """Diff-friendly markdown report from existing topics JSONs (no API calls)."""
    lines = ["# Publication Topics Extraction Report\n"]
    for paper_dir, pub, _pdf, _strains in sorted(papers, key=lambda e: e[1].get("papername", "")):
        data = load_topics(paper_dir)
        if not data:
            continue
        genes = data.get("genes", [])
        pathways = data.get("pathways", [])
        meta = data.get("metadata", {})
        lines.append(f"## {pub.get('papername', paper_dir.name)}")
        unc = meta.get("uncaptured_identifiers", "?")
        unc_note = meta.get("uncaptured_identifiers_note", "")
        unc_str = f"uncaptured_ids={unc}" + (f" ({unc_note})" if unc not in ("none", "?") and unc_note else "")
        lines.append(
            f"_{len(genes)} genes, {len(pathways)} pathways; "
            f"confidence={meta.get('self_assessment', '?')}; {unc_str}_\n"
        )
        for g in sorted(genes, key=lambda x: (x.get("prominence") != "central",
                                              x.get("gene_name", ""))):
            ids = ", ".join(g.get("identifiers", [])) or "—"
            lines.append(
                f"- **{g.get('gene_name') or g.get('surface_form')}** "
                f"[{g.get('prominence')}] strain={g.get('strain')} ids=[{ids}]"
            )
        if pathways:
            lines.append("\n  Pathways:")
            for p in pathways:
                pid = p.get("kegg_pathway_id") or "—"
                lines.append(f"  - {p.get('surface_form')} ({pid}) [{p.get('prominence')}]")
        lines.append("")
    return "\n".join(lines)


# ── CLI ──


def main():
    from dotenv import load_dotenv
    load_dotenv()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    parser = argparse.ArgumentParser(description="Publication topics extraction pipeline")
    parser.add_argument("--paper", help="Filter by paper name (partial match)")
    parser.add_argument("--model", default=os.environ.get("TOPICS_EXTRACTION_MODEL", "gpt-5.4"))
    parser.add_argument("--flex", action="store_true", help="Flex processing (cheaper)")
    parser.add_argument("--force", action="store_true", help="Overwrite existing")
    parser.add_argument("--dry-run", action="store_true", help="Show what would run")
    parser.add_argument("--report", action="store_true", help="Generate report (no API calls)")
    parser.add_argument("--chunk-pages", type=int, default=15,
                        help="Max pages per request chunk (reference pages skipped)")
    args = parser.parse_args()

    papers = find_all_papers()
    if args.paper:
        papers = [p for p in papers
                  if args.paper.lower() in p[1].get("papername", "").lower()]

    if args.report:
        report = generate_report(papers)
        report_path = Path("data/publication_topics_report.md")
        report_path.write_text(report)
        logger.info("Report written to %s", report_path)
        return

    if not papers:
        logger.warning("No matching papers found")
        return

    if args.dry_run:
        print(f"\nWould process {len(papers)} papers:\n")
        for paper_dir, pub, pdf_path, strains in papers:
            status = "EXISTS" if has_topics(paper_dir) else "NEW"
            print(f"  {pub.get('papername', paper_dir.name)}: "
                  f"{len(strains)} strain(s) [{status}] <- {pdf_path.name}")
        return

    from openai import OpenAI
    client = OpenAI()
    total_in = total_out = total_genes = total_pathways = 0

    for paper_dir, pub, pdf_path, strains in papers:
        paper = pub.get("papername", paper_dir.name)
        if has_topics(paper_dir) and not args.force:
            logger.info("%s: exists, skipping", paper)
            continue

        logger.info("Paper: %s (%d strain(s))", paper, len(strains))
        try:
            merged, usage, n_chunks = extract_topics_chunked(
                client, pdf_path, strains, model=args.model, flex=args.flex,
                chunk_pages=args.chunk_pages,
            )
            genes = merged["genes"]
            pathways = merged["pathways"]
            metadata = {
                "paper": paper,
                "doi": pub.get("doi"),
                "strains": strains,
                "model": args.model,
                "n_chunks": n_chunks,
                "chunk_pages": args.chunk_pages,
                "self_assessment": merged["self_assessment"],
                "uncaptured_identifiers": merged["uncaptured_identifiers"],
                "uncaptured_identifiers_note": merged["uncaptured_identifiers_note"],
                "extracted_at": datetime.now().isoformat(),
                "input_tokens": usage["input_tokens"],
                "output_tokens": usage["output_tokens"],
            }
            save_topics(paper_dir, metadata, genes, pathways)
            total_in += usage["input_tokens"]
            total_out += usage["output_tokens"]
            total_genes += len(genes)
            total_pathways += len(pathways)
            logger.info("    %d chunks → %d genes, %d pathways [%s] | uncaptured_ids=%s",
                        n_chunks, len(genes), len(pathways), merged["self_assessment"],
                        merged["uncaptured_identifiers"])
        except Exception as e:
            logger.error("  %s: FAILED — %s", paper, e)
        finally:
            time.sleep(1)

    logger.info("Done: %d genes + %d pathways, %d in + %d out tokens",
                total_genes, total_pathways, total_in, total_out)


if __name__ == "__main__":
    main()
