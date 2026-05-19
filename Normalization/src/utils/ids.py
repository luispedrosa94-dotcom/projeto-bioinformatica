"""
UniProt ID parsing utilities.

Across the tools in this project, protein identifiers appear in three variants:

    sp|A0B9K2|AATA_METTP                                          (UPIMAPI, reCOGnizer, eggNOG, DeepFRI, DeepGO2, CLEAN)
    sp_A0B9K2_AATA_METTP                                          (Foldseek query without suffix; ColabFold filenames)
    sp_A0B9K2_AATA_METTP_unrelaxed_rank_001_alphafold2_model_...  (Foldseek query with model suffix)

All carry a UniProt accession we can extract with one regex. We treat the
accession as the canonical key for the rest of the pipeline.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Optional

# UniProt accession regex (official spec: https://www.uniprot.org/help/accession_numbers).
# 6 chars: [OPQ][0-9][A-Z0-9]{3}[0-9]   OR   [A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}
# We use a permissive pattern that also accepts the 10-char form.
_ACCESSION_RE = re.compile(
    r"\b("
    r"[OPQ][0-9][A-Z0-9]{3}[0-9]"
    r"|"
    r"[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}"
    r")\b"
)

# Canonical pipe-form: sp|ACCESSION|ENTRY_NAME or tr|ACCESSION|ENTRY_NAME
_PIPE_RE = re.compile(r"^(sp|tr)\|([A-Z0-9]+)\|([A-Z0-9_]+)$", re.IGNORECASE)
# Underscore form (Foldseek/ColabFold): sp_ACCESSION_ENTRY_NAME[_extra...]
_UNDERSCORE_RE = re.compile(r"^(sp|tr)_([A-Z0-9]+)_([A-Z0-9]+(?:_[A-Z0-9]+)?)", re.IGNORECASE)


@dataclass(frozen=True)
class ParsedID:
    accession: str        # canonical UniProt accession
    db_source: str        # 'sp' | 'tr' | 'unknown'
    entry_name: Optional[str]   # e.g. AATA_METTP
    original: str         # the original string we received


def parse_protein_id(raw_id: str) -> ParsedID:
    """Parse any of the ID variants used in this project.

    Strategy:
      1. Try canonical pipe form (most tools use this).
      2. Try underscore form (Foldseek, ColabFold).
      3. Fall back to extracting the first valid UniProt accession.
      4. If nothing matches, return the raw string as accession with db_source='unknown'.
    """
    s = raw_id.strip()

    m = _PIPE_RE.match(s)
    if m:
        return ParsedID(
            accession=m.group(2).upper(),
            db_source=m.group(1).lower(),
            entry_name=m.group(3).upper(),
            original=raw_id,
        )

    m = _UNDERSCORE_RE.match(s)
    if m:
        # Reconstruct entry_name from groups 3 onward up to a known suffix marker
        entry = m.group(3).upper()
        return ParsedID(
            accession=m.group(2).upper(),
            db_source=m.group(1).lower(),
            entry_name=entry,
            original=raw_id,
        )

    # Last resort: any UniProt-shaped accession anywhere in the string
    m = _ACCESSION_RE.search(s.upper())
    if m:
        return ParsedID(
            accession=m.group(1),
            db_source="unknown",
            entry_name=None,
            original=raw_id,
        )

    return ParsedID(
        accession=s.upper(),
        db_source="unknown",
        entry_name=None,
        original=raw_id,
    )
