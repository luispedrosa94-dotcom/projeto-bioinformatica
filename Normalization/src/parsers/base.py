"""
Base parser interface.

Each tool gets its own parser module (upimapi.py, eggnog.py, ...) that
implements the `Parser` protocol. The orchestrator runs each parser, gets
back a list of AnnotationRecord objects, and concatenates them all into
the final long-format annotations table.
"""
from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Iterable

from ..schema import AnnotationRecord, ProteinRecord


class BaseParser(ABC):
    """Abstract parser. Subclasses parse one tool's output."""

    #: short name used in logs and config (e.g. 'upimapi')
    tool_name: str = ""

    def __init__(self, raw_data_root: Path):
        self.raw_data_root = Path(raw_data_root)
        self.log = logging.getLogger(f"parser.{self.tool_name}")

    @abstractmethod
    def parse(self) -> list[AnnotationRecord]:
        """Read the tool's output files and return canonical AnnotationRecord objects."""
        raise NotImplementedError

    def discover_proteins(self) -> list[ProteinRecord]:
        """Optional: return ProteinRecord objects for proteins seen by this parser.

        The default implementation returns an empty list; the orchestrator will
        deduplicate proteins across all parsers anyway. Override only if a parser
        has unique knowledge about a protein (e.g. ColabFold knowing the length).
        """
        return []


def filter_to_scope(
    records: Iterable[AnnotationRecord],
    scope: set[str] | None,
) -> list[AnnotationRecord]:
    """Filter records to only those whose accession is in scope. None = keep all."""
    if scope is None:
        return list(records)
    return [r for r in records if r.uniprot_accession in scope]
