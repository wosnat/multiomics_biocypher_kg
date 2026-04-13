"""OpenAI Responses API client for timepoint extraction.

Mirrors the cluster-extraction pattern. Attaches the paper PDF + any
extras.additional_pdfs as multimodal inputs. Returns a parsed dict plus
token/model metadata.
"""
from __future__ import annotations

import json
import logging
import os
from pathlib import Path

from openai import OpenAI

logger = logging.getLogger(__name__)

DEFAULT_MODEL = "gpt-4.1-mini"


class OpenAIResponsesClient:
    """Minimal wrapper around OpenAI Responses API for structured-JSON extraction."""

    def __init__(self, api_key: str | None = None, model: str = DEFAULT_MODEL):
        self._client = OpenAI(api_key=api_key or os.environ.get("OPENAI_API_KEY"))
        self._model = model

    def call(
        self,
        prompt: str,
        pdf_paths: list[Path],
        model: str | None = None,
    ) -> tuple[dict, dict]:
        """Call the LLM with the prompt + attached PDFs; parse JSON response."""
        model = model or self._model

        # Upload PDFs that exist. Missing PDFs are a warning, not a hard fail.
        file_inputs: list[dict] = []
        uploaded_ids: list[str] = []
        try:
            for p in pdf_paths:
                if not p.exists():
                    logger.warning("PDF not found, skipping: %s", p)
                    continue
                with open(p, "rb") as fh:
                    uploaded = self._client.files.create(file=fh, purpose="user_data")
                uploaded_ids.append(uploaded.id)
                file_inputs.append({"type": "input_file", "file_id": uploaded.id})

            content = file_inputs + [{"type": "input_text", "text": prompt}]

            response = self._client.responses.create(
                model=model,
                input=[{"role": "user", "content": content}],
                text={"format": {"type": "json_object"}},
            )

            raw_text = response.output_text
            try:
                payload = json.loads(raw_text)
            except json.JSONDecodeError as e:
                raise ValueError(f"LLM returned unparseable JSON: {e}\n{raw_text[:500]}")

            meta = {
                "input_tokens": response.usage.input_tokens if response.usage else 0,
                "output_tokens": response.usage.output_tokens if response.usage else 0,
                "model": model,
            }
            return payload, meta
        finally:
            for fid in uploaded_ids:
                try:
                    self._client.files.delete(fid)
                except Exception:
                    logger.warning("Failed to delete uploaded file %s", fid)
