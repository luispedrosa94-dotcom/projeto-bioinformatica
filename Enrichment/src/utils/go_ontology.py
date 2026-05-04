"""
Local GO ontology parser.

Downloads go-basic.obo once from the GO consortium and caches it locally.
Parses it to build a GO term → aspect (BP/MF/CC) map — no API calls needed.

File: http://purl.obolibrary.org/obo/go/go-basic.obo
Size: ~32 MB, downloaded once, cached permanently.
"""
from __future__ import annotations

import logging
import re
from pathlib import Path

import requests

log = logging.getLogger(__name__)

GO_OBO_URL   = "http://purl.obolibrary.org/obo/go/go-basic.obo"
DEFAULT_CACHE = Path(__file__).resolve().parents[3] / "outputs" / "go-basic.obo"

_NAMESPACE_MAP = {
    "biological_process": "BP",
    "molecular_function": "MF",
    "cellular_component": "CC",
}


def load_go_aspects(cache_path: Path = DEFAULT_CACHE) -> dict[str, str]:
    """
    Returns a dict mapping GO term ID → aspect short code (BP/MF/CC).

    Downloads go-basic.obo if not cached. Subsequent calls use the cache.
    """
    obo_text = _get_obo(cache_path)
    return _parse_obo(obo_text)


def _get_obo(cache_path: Path) -> str:
    if cache_path.exists():
        log.info("GO ontology: loading from cache %s", cache_path)
        return cache_path.read_text(encoding="utf-8", errors="replace")

    log.info("GO ontology: downloading from %s ...", GO_OBO_URL)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    resp = requests.get(GO_OBO_URL, timeout=120, stream=True)
    resp.raise_for_status()

    content = b""
    for chunk in resp.iter_content(chunk_size=1024 * 64):
        content += chunk

    text = content.decode("utf-8", errors="replace")
    cache_path.write_text(text, encoding="utf-8")
    log.info("GO ontology: saved to %s (%d KB)", cache_path, len(content) // 1024)
    return text


def _parse_obo(text: str) -> dict[str, str]:
    """Parse OBO format and return {go_id: aspect_short}."""
    aspects: dict[str, str] = {}
    current_id = None
    current_ns  = None
    in_term     = False

    for line in text.splitlines():
        line = line.strip()
        if line == "[Term]":
            in_term     = True
            current_id  = None
            current_ns  = None
        elif line == "" and in_term:
            # End of a [Term] block
            if current_id and current_ns:
                short = _NAMESPACE_MAP.get(current_ns)
                if short:
                    aspects[current_id] = short
            in_term    = False
            current_id = None
            current_ns = None
        elif in_term:
            if line.startswith("id: GO:"):
                current_id = line[4:].strip()
            elif line.startswith("namespace: "):
                current_ns = line[len("namespace: "):].strip()
            elif line.startswith("alt_id: GO:"):
                # Also map alternative IDs to the same namespace
                alt_id = line[7:].strip()
                if current_ns:
                    short = _NAMESPACE_MAP.get(current_ns)
                    if short:
                        aspects[alt_id] = short

    log.info("GO ontology: parsed %d terms", len(aspects))
    return aspects
