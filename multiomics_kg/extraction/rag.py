# multiomics_kg/extraction/rag.py
"""Shared RAG utilities: text chunking, embedding, cosine retrieval."""
import logging
import re
from typing import Optional

import numpy as np

logger = logging.getLogger(__name__)

try:
    import openai as _openai
except Exception:
    _openai = None


def chunk_text(text: str, chunk_size: int = 400, overlap: int = 100) -> list[str]:
    """Split text into overlapping chunks at sentence boundaries."""
    sentences = re.split(r'(?<=[.!?])\s+', text)
    chunks = []
    current = ""
    for sent in sentences:
        if len(current) + len(sent) > chunk_size and current:
            chunks.append(current.strip())
            words = current.split()
            keep = max(1, overlap // 5)
            overlap_text = " ".join(words[-keep:]) if len(words) > keep else ""
            current = overlap_text + " " + sent
        else:
            current = (current + " " + sent).strip()
    if current.strip():
        chunks.append(current.strip())
    return chunks


def embed_texts(texts: list[str],
                model: str = "text-embedding-3-small",
                ) -> list[list[float]]:
    """Embed a list of texts using OpenAI embeddings API.

    Handles batching for large lists (API limit ~2048 per call).
    """
    if _openai is None or not texts:
        return []
    client = _openai.OpenAI()
    all_embeddings = []
    batch_size = 2000
    for i in range(0, len(texts), batch_size):
        batch = texts[i:i + batch_size]
        response = client.embeddings.create(model=model, input=batch)
        all_embeddings.extend(item.embedding for item in response.data)
    return all_embeddings


def retrieve_top_k(query: str,
                   chunks: list[str],
                   chunk_embeddings: list[list[float]],
                   top_k: int = 8,
                   model: str = "text-embedding-3-small",
                   ) -> list[tuple[str, float]]:
    """Embed a query and retrieve top-k most similar chunks.

    Returns list of (chunk_text, relevance_score) sorted by score desc.
    """
    if not chunk_embeddings or not chunks:
        return []
    # Embed query
    query_emb = embed_texts([query], model=model)
    if not query_emb:
        return []
    q = np.array(query_emb[0])
    matrix = np.array(chunk_embeddings)
    sims = matrix @ q / (np.linalg.norm(matrix, axis=1) * np.linalg.norm(q) + 1e-10)
    top_indices = np.argsort(sims)[-top_k:][::-1]
    return [(chunks[i], float(sims[i])) for i in top_indices]
