"""Interactive RAG query experiments for cluster extraction.

Usage:
    uv run python scripts/rag_experiment.py "data/.../tolonen 2006" --query "MED4 cluster 1"
    uv run python scripts/rag_experiment.py "data/.../tolonen 2006" --queries-file queries.txt
    uv run python scripts/rag_experiment.py "data/.../tolonen 2006" --interactive

Caches PDF text, chunks, and embeddings to {paper_dir}/.extraction_cache/shared/
so subsequent queries are instant (no API calls).
"""

import argparse
import json
import logging
import sys
from pathlib import Path

from dotenv import load_dotenv
load_dotenv()

import numpy as np

# Ensure project root on path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from multiomics_kg.extraction.pdf_utils import extract_pdf_text, collect_pdf_files
from multiomics_kg.extraction.rag import chunk_text, embed_texts, retrieve_top_k

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)


def get_or_build_cache(paper_dir: Path) -> tuple[list[str], list[list[float]]]:
    """Load cached chunks+embeddings, or build and cache them."""
    cache_dir = paper_dir / ".extraction_cache" / "shared"
    chunks_path = cache_dir / "chunks.json"
    embeddings_path = cache_dir / "embeddings.npy"
    pdf_text_path = cache_dir / "pdf_text.json"

    # Try loading from cache
    if chunks_path.exists() and embeddings_path.exists():
        logger.info("Loading cached chunks and embeddings from %s", cache_dir)
        with open(chunks_path) as f:
            chunks = json.load(f)
        embeddings = np.load(embeddings_path).tolist()
        logger.info("  %d chunks, %d embeddings loaded", len(chunks), len(embeddings))
        return chunks, embeddings

    # Build from scratch
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Extract PDF text
    if pdf_text_path.exists():
        logger.info("Loading cached PDF text")
        with open(pdf_text_path) as f:
            pdf_texts = json.load(f)
    else:
        logger.info("Extracting PDF text...")
        main_pdf = None
        pc_path = paper_dir / "paperconfig.yaml"
        if pc_path.exists():
            import yaml
            with open(pc_path) as f:
                config = yaml.safe_load(f)
            pdf_rel = config.get("publication", {}).get("papermainpdf", "")
            if pdf_rel:
                project_root = Path(__file__).resolve().parent.parent
                candidate = project_root / pdf_rel
                if candidate.exists():
                    main_pdf = candidate

        pdf_texts = {}
        for pdf_path in collect_pdf_files(paper_dir, main_pdf):
            text = extract_pdf_text(pdf_path)
            if text:
                pdf_texts[str(pdf_path.name)] = text
                logger.info("  %s: %d chars", pdf_path.name, len(text))

        with open(pdf_text_path, "w") as f:
            json.dump(pdf_texts, f)

    # Chunk
    combined = "\n\n".join(pdf_texts.values())
    chunks = chunk_text(combined)
    logger.info("Chunked into %d chunks from %d chars", len(chunks), len(combined))

    with open(chunks_path, "w") as f:
        json.dump(chunks, f)

    # Embed
    logger.info("Embedding %d chunks (API call)...", len(chunks))
    embeddings = embed_texts(chunks)
    logger.info("  Got %d embeddings", len(embeddings))

    np.save(embeddings_path, np.array(embeddings))

    return chunks, embeddings


def run_query(query: str, chunks: list[str], embeddings: list[list[float]],
              top_k: int = 8) -> list[tuple[str, float]]:
    """Run a single RAG query and return results."""
    return retrieve_top_k(query, chunks, embeddings, top_k=top_k)


def print_results(query: str, results: list[tuple[str, float]]):
    """Pretty-print RAG results."""
    print(f"\n{'='*80}")
    print(f"QUERY: {query}")
    print(f"{'='*80}")
    for i, (text, score) in enumerate(results):
        print(f"\n--- [{i+1}] relevance={score:.3f} ---")
        print(text[:500])
        if len(text) > 500:
            print(f"  ... ({len(text)} chars total)")


def interactive_mode(chunks, embeddings, top_k=8):
    """Interactive query loop."""
    print(f"\nInteractive RAG query mode ({len(chunks)} chunks indexed)")
    print("Type a query, or 'quit' to exit.\n")
    while True:
        try:
            query = input("query> ").strip()
        except (EOFError, KeyboardInterrupt):
            break
        if not query or query.lower() in ("quit", "exit", "q"):
            break
        results = run_query(query, chunks, embeddings, top_k=top_k)
        print_results(query, results)


def main():
    parser = argparse.ArgumentParser(description="RAG query experiments")
    parser.add_argument("paper_dir", type=Path, help="Path to paper directory")
    parser.add_argument("--query", "-q", help="Single query to run")
    parser.add_argument("--interactive", "-i", action="store_true",
                        help="Interactive query mode")
    parser.add_argument("--top-k", "-k", type=int, default=8,
                        help="Number of results (default: 8)")
    parser.add_argument("--batch", nargs="+",
                        help="Multiple queries to run")
    args = parser.parse_args()

    paper_dir = args.paper_dir.resolve()
    if not paper_dir.exists():
        print(f"Paper directory not found: {paper_dir}")
        sys.exit(1)

    chunks, embeddings = get_or_build_cache(paper_dir)

    if args.query:
        results = run_query(args.query, chunks, embeddings, top_k=args.top_k)
        print_results(args.query, results)

    if args.batch:
        for q in args.batch:
            results = run_query(q, chunks, embeddings, top_k=args.top_k)
            print_results(q, results)

    if args.interactive:
        interactive_mode(chunks, embeddings, top_k=args.top_k)

    if not args.query and not args.batch and not args.interactive:
        # Default: run a few test queries
        test_queries = [
            "MED4 cluster 1",
            'MED4 "cluster 1"',
            "cluster 1 nitrogen transport upregulated",
            "cluster 1 urtA cynA transport",
        ]
        for q in test_queries:
            results = run_query(q, chunks, embeddings, top_k=args.top_k)
            print_results(q, results)


if __name__ == "__main__":
    main()
