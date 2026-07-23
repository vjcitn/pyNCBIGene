"""BiocFileCache helpers for live and frozen parquet caching."""

import re
from typing import Optional

from ._state import get_cache, taxid_column


def _live_key(gres: str, taxid: int) -> str:
    return f"{gres}_taxid{taxid}.parquet"


def _frozen_key(gres: str, taxid: int, tag: str) -> str:
    return f"{gres}_taxid{taxid}_frozen_{tag}.parquet"


def _bfc_query(pattern: str):
    """Return cache entries whose rname matches pattern."""
    cache = get_cache()
    try:
        # pyBiocFileCache >= 0.2 API
        return cache.query(query=pattern)
    except TypeError:
        return cache.query(pattern)


def _bfc_remove(rids):
    cache = get_cache()
    for rid in (rids if hasattr(rids, "__iter__") and not isinstance(rids, str) else [rids]):
        try:
            cache.remove(rid)
        except Exception:
            pass


def _lookup_cache_path(gres: str, taxid: int) -> Optional[str]:
    key = _live_key(gres, taxid)
    hits = _bfc_query(key)
    if hits is None or len(hits) == 0:
        return None
    # exclude frozen entries
    hits = [h for h in hits if "_frozen_" not in h.get("rname", "")]
    if not hits:
        return None
    return hits[0].get("rpath") or hits[0].get("fpath")


def _lookup_frozen_path(gres: str, taxid: int, tag: str) -> str:
    key = _frozen_key(gres, taxid, tag)
    hits = _bfc_query(key)
    if not hits or len(hits) == 0:
        raise KeyError(
            f"No frozen snapshot '{tag}' for resource '{gres}' taxid {taxid}. "
            f"Call freeze_taxon_cache({taxid}, tag='{tag}') first."
        )
    return hits[0].get("rpath") or hits[0].get("fpath")
