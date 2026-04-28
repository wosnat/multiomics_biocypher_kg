"""Resolver accessor API.

Phase 1.1A skeleton: every function raises NotImplementedError. Real
implementations land in Phase 1.1B once the actual MNX file shapes are
confirmed via the audit.

API contract (consumed by download/build_gene_annotations.py via the
transforms framework, by the Spec 1.2 scaffold builder, and by Phase 2
paper-measurement extraction):

    open_resolver(path: Path | None = None) -> sqlite3.Connection
    resolve_metabolite(value: str, conn) -> tuple[str | None, str]
    resolve_reaction (value: str, conn) -> tuple[str | None, str]
"""
from __future__ import annotations

import sqlite3
from pathlib import Path


def open_resolver(path: Path | None = None) -> sqlite3.Connection:
    raise NotImplementedError("Phase 1.1B")


def resolve_metabolite(value: str, conn: sqlite3.Connection) -> tuple[str | None, str]:
    raise NotImplementedError("Phase 1.1B")


def resolve_reaction(value: str, conn: sqlite3.Connection) -> tuple[str | None, str]:
    raise NotImplementedError("Phase 1.1B")
