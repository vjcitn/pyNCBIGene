"""Taxon-level caching and reproducibility freezing."""

import re
import shutil
import tempfile
from pathlib import Path
from typing import Optional

import pandas as pd

from ._state import get_cache, get_connection, taxid_column
from ._cache import _live_key, _frozen_key, _bfc_query, _bfc_remove
from .remote import available_ncbi_parquet, open_ncbi_gene


def cache_by_taxon(
    taxid: int,
    resources: Optional[list] = None,
    force: bool = False,
    verbose: bool = True,
) -> dict:
    """Filter and cache remote parquet resources for a single taxon locally.

    After calling this function, :func:`open_ncbi_gene` routes to the local
    file for the specified taxon with no API change.

    Parameters
    ----------
    taxid : int
        NCBI taxonomy ID (e.g. ``9606`` for human, ``10090`` for mouse).
    resources : list or None
        Resource names (without ``.parquet``) to cache. Defaults to all.
    force : bool
        If ``True``, re-download and overwrite existing cache entries.
    verbose : bool
        Print progress messages.

    Returns
    -------
    dict
        ``{resource_name: local_path}`` for each cached resource.
    """
    if resources is None:
        resources = [r.replace(".parquet", "") for r in available_ncbi_parquet()]

    cache = get_cache()
    con = get_connection()
    paths = {}

    for gres in resources:
        key = _live_key(gres, taxid)
        hits = _bfc_query(key)
        live_hits = [h for h in (hits or []) if "_frozen_" not in h.get("rname", "")]

        if live_hits and not force:
            if verbose:
                print(f"  {key}: already cached")
            paths[gres] = live_hits[0].get("rpath") or live_hits[0].get("fpath")
            continue

        if verbose:
            print(f"  {gres}: filtering from OSN ...")

        # create or refresh the remote view
        open_ncbi_gene(gres)
        vname = f"v_{re.sub(r'[^A-Za-z0-9]', '_', gres)}"
        tcol = taxid_column(gres)

        tmp = tempfile.mktemp(suffix=".parquet")
        con.execute(
            f"COPY (SELECT * FROM {vname} WHERE \"{tcol}\" = {taxid}) "
            f"TO '{tmp}' (FORMAT parquet, COMPRESSION zstd, COMPRESSION_LEVEL 15)"
        )

        if live_hits and force:
            _bfc_remove([h["rid"] for h in live_hits])

        try:
            result = cache.add(rname=key, fpath=tmp, rtype="local")
        except TypeError:
            result = cache.add(key, tmp)
        Path(tmp).unlink(missing_ok=True)

        local_path = result if isinstance(result, str) else result.get("rpath", "")
        paths[gres] = local_path
        if verbose:
            print(f"  {key}: done")

    return paths


def taxon_cache_info(taxid: Optional[int] = None) -> pd.DataFrame:
    """List taxon-specific parquets stored in the local cache.

    Parameters
    ----------
    taxid : int or None
        When ``None``, lists all cached taxon parquets.

    Returns
    -------
    pd.DataFrame
        With columns ``rname``, ``rpath``, and ``create_time`` (or equivalent).
    """
    pattern = f"_taxid{taxid}" if taxid is not None else "_taxid"
    hits = _bfc_query(pattern) or []
    rows = []
    for h in hits:
        rows.append({
            "rname": h.get("rname", ""),
            "rpath": h.get("rpath") or h.get("fpath", ""),
            "create_time": h.get("create_time") or h.get("access_time", ""),
        })
    return pd.DataFrame(rows)


def clear_taxon_cache(taxid: Optional[int] = None) -> list:
    """Remove live (non-frozen) taxon parquets from the local cache.

    Frozen snapshots are never removed by this function.

    Parameters
    ----------
    taxid : int or None
        When ``None``, removes live cache for all taxa.

    Returns
    -------
    list
        Names of removed cache entries.
    """
    pattern = f"_taxid{taxid}" if taxid is not None else "_taxid"
    hits = _bfc_query(pattern) or []
    live = [h for h in hits if "_frozen_" not in h.get("rname", "")]
    if not live:
        print("No live cached entries found.")
        return []
    _bfc_remove([h["rid"] for h in live])
    names = [h["rname"] for h in live]
    print(f"Removed {len(names)} live cached parquet(s).")
    return names


def freeze_taxon_cache(
    taxid: int,
    tag: str,
    resources: Optional[list] = None,
    force: bool = False,
) -> dict:
    """Freeze a snapshot of taxon-cached parquets for reproducibility.

    Makes a physical copy of each live cached parquet under ``tag``.  Frozen
    copies survive :func:`clear_taxon_cache`.  Use
    ``open_ncbi_gene(resource, taxid, freeze_tag=tag)`` to query a snapshot.

    Parameters
    ----------
    taxid : int
        NCBI taxonomy ID.
    tag : str
        Identifying label; must not already exist unless ``force=True``.
    resources : list or None
        Resources to freeze. Defaults to all currently live-cached resources
        for this taxon.
    force : bool
        If ``True``, overwrite an existing frozen snapshot with the same tag.

    Returns
    -------
    dict
        ``{resource_name: frozen_local_path}``
    """
    if not tag:
        raise ValueError("tag must be a non-empty string")

    cache = get_cache()

    # reject duplicate tag
    existing_tag = _bfc_query(f"_frozen_{tag}.parquet") or []
    if existing_tag and not force:
        raise ValueError(
            f"Tag '{tag}' already exists in cache. Use force=True to overwrite."
        )

    # default: whatever is live-cached for this taxid
    if resources is None:
        pattern = f"_taxid{taxid}.parquet"
        all_live = [
            h for h in (_bfc_query(pattern) or [])
            if "_frozen_" not in h.get("rname", "")
        ]
        if not all_live:
            raise RuntimeError(
                f"No live cache for taxid {taxid}. "
                f"Call cache_by_taxon({taxid}) first."
            )
        resources = [
            re.sub(rf"_taxid{taxid}\.parquet$", "", h["rname"])
            for h in all_live
        ]

    paths = {}
    for gres in resources:
        live_key = _live_key(gres, taxid)
        live = [
            h for h in (_bfc_query(live_key) or [])
            if "_frozen_" not in h.get("rname", "")
        ]
        if not live:
            raise RuntimeError(
                f"No live cache for '{gres}' taxid {taxid}. "
                f"Call cache_by_taxon({taxid}) first."
            )

        frozen_key = _frozen_key(gres, taxid, tag)
        old = _bfc_query(frozen_key) or []
        if old and force:
            _bfc_remove([h["rid"] for h in old])

        src = live[0].get("rpath") or live[0].get("fpath", "")
        tmp = tempfile.mktemp(suffix=".parquet")
        shutil.copy2(src, tmp)

        try:
            result = cache.add(rname=frozen_key, fpath=tmp, rtype="local")
        except TypeError:
            result = cache.add(frozen_key, tmp)
        Path(tmp).unlink(missing_ok=True)

        local_path = result if isinstance(result, str) else result.get("rpath", "")
        paths[gres] = local_path
        print(f"  frozen: {frozen_key}")

    return paths
