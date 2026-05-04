"""
Checkpoint utilities — save and load progress so the pipeline can resume
after interruption without redoing completed work.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path

log = logging.getLogger(__name__)


def save(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False)
    tmp.replace(path)
    log.debug("Checkpoint saved → %s", path)


def load(path: Path) -> dict | None:
    if path.exists():
        try:
            with open(path, encoding="utf-8") as f:
                data = json.load(f)
            log.info("Checkpoint loaded from %s", path)
            return data
        except Exception as e:
            log.warning("Checkpoint corrupted (%s) — ignoring", e)
    return None
