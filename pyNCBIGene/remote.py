"""Lazy remote queries over NCBI Gene parquet in the OSN bucket."""

import re
import urllib.request
from typing import Optional

import pandas as pd

from ._state import (
    OSN_BASE,
    OSN_LISTING_URL,
    OSN_PROVENANCE_URL,
    get_connection,
    taxid_column,
)

# session-level cache of available resource names (strings, not user input)
_resource_names: list[str] | None = None


def available_ncbi_parquet() -> list[str]:
    """List parquet resources currently in the OSN bucket."""
    with urllib.request.urlopen(OSN_LISTING_URL) as resp:
        xml = resp.read().decode()
    keys = re.findall(r"<Key>([^<]+\.parquet)</Key>", xml)
    if not keys:
        raise RuntimeError(
            "Bucket listing returned no parquet files -- check the OSN URL."
        )
    return sorted(re.sub(r".*/", "", k) for k in keys)


def _resource_name_list() -> list[str]:
    global _resource_names
    if _resource_names is None:
        _resource_names = [r.replace(".parquet", "") for r in available_ncbi_parquet()]
    return _resource_names


def _validate_resource(gres: str) -> None:
    if gres not in _resource_name_list():
        raise ValueError(
            f"'{gres}' is not an available resource. "
            "Call available_ncbi_parquet() to see options."
        )


def _validate_column(col: str, gres: str) -> None:
    """Raise ValueError if col is not a valid column name for gres."""
    valid = ncbi_gene_fields(gres)["column_name"].tolist()
    if col not in valid:
        raise ValueError(
            f"Column '{col}' not found in '{gres}'. "
            f"Available: {valid}"
        )


def _safe_vname(gres: str) -> str:
    """Return the duckdb VIEW name for a resource (always derived, never user-supplied)."""
    return "v_" + re.sub(r"[^A-Za-z0-9]", "_", gres)


def _ensure_view(gres: str) -> str:
    """Ensure the remote VIEW for gres exists; return its name."""
    _validate_resource(gres)
    con = get_connection()
    url = OSN_BASE.format(gres)
    vname = _safe_vname(gres)
    con.execute(
        f"CREATE OR REPLACE VIEW {vname} AS "
        f"SELECT * FROM read_parquet('{url}')"
    )
    return vname


def _cached_parquet_path(gres: str, taxid: Optional[int]) -> Optional[str]:
    if taxid is None:
        return None
    try:
        from ._cache import _lookup_cache_path
        return _lookup_cache_path(gres, taxid)
    except Exception:
        return None


def _frozen_parquet_path(gres: str, taxid: Optional[int], tag: str) -> str:
    from ._cache import _lookup_frozen_path
    return _lookup_frozen_path(gres, taxid, tag)


def _resolve_source(
    gres: str,
    taxid: Optional[int],
    freeze_tag: Optional[str],
) -> tuple[str, bool]:
    """Return (view_name, taxid_already_applied) for a resource request.

    Creates the appropriate duckdb VIEW if it does not yet exist.

    * frozen snapshot  -> VIEW over local parquet, taxid pre-filtered
    * live local cache -> VIEW over local parquet, taxid pre-filtered
    * remote OSN       -> VIEW over remote parquet, taxid NOT applied
    """
    con = get_connection()

    if freeze_tag is not None:
        local = _frozen_parquet_path(gres, taxid, freeze_tag)
        vname = "v_frozen_" + re.sub(r"[^A-Za-z0-9]", "_", gres) + f"_{taxid}_{freeze_tag}"
        con.execute(
            f"CREATE OR REPLACE VIEW {vname} AS "
            f"SELECT * FROM read_parquet('{local}')"
        )
        return vname, True

    local = _cached_parquet_path(gres, taxid)
    if local is not None:
        vname = "v_local_" + re.sub(r"[^A-Za-z0-9]", "_", gres) + f"_{taxid}"
        con.execute(
            f"CREATE OR REPLACE VIEW {vname} AS "
            f"SELECT * FROM read_parquet('{local}')"
        )
        return vname, True

    return _ensure_view(gres), False


def open_ncbi_gene(
    resource: str = "gene_info",
    taxid: Optional[int] = None,
    freeze_tag: Optional[str] = None,
) -> "duckdb.DuckDBPyRelation":
    """Open a lazy duckdb relation over a remote NCBI Gene parquet resource.

    Parameters
    ----------
    resource : str
        Resource name with or without ``.parquet`` suffix.
    taxid : int or None
        When provided, filters to this NCBI taxonomy ID and drops the taxid
        column from the result.
    freeze_tag : str or None
        When provided, opens a frozen snapshot.  Raises ``KeyError`` if not
        found.

    Returns
    -------
    duckdb.DuckDBPyRelation
        Lazy relation -- chain ``.filter()``, ``.select()``, ``.df()`` etc.
    """
    gres = resource.replace(".parquet", "")
    con = get_connection()
    tcol = taxid_column(gres)

    vname, taxid_applied = _resolve_source(gres, taxid, freeze_tag)
    all_cols = con.table(vname).columns
    cols = [c for c in all_cols if c != tcol]
    quoted = ", ".join(f'"{c}"' for c in cols)

    if taxid is not None and not taxid_applied:
        return con.sql(
            f'SELECT {quoted} FROM {vname} WHERE "{tcol}" = {int(taxid)}'
        )
    if taxid is not None:
        return con.sql(f"SELECT {quoted} FROM {vname}")
    return con.table(vname)


def ncbi_gene_fields(resource: str = "gene_info") -> pd.DataFrame:
    """Return column names and types for a resource.

    Parameters
    ----------
    resource : str
        Resource name with or without ``.parquet`` suffix.

    Returns
    -------
    pd.DataFrame
        Data frame with columns ``column_name`` and ``column_type``.
    """
    gres = resource.replace(".parquet", "")
    vname = _ensure_view(gres)
    result = get_connection().execute(f"DESCRIBE {vname}").df()
    return result[["column_name", "column_type"]]


def ncbi_parquet_info() -> pd.DataFrame:
    """Retrieve size, upload date, and NCBI source date for each bucket resource."""
    with urllib.request.urlopen(OSN_LISTING_URL) as resp:
        xml = resp.read().decode()

    blocks = xml.split("<Contents>")[1:]
    rows = []
    for block in blocks:
        key = re.search(r"<Key>([^<]+)</Key>", block)
        size = re.search(r"<Size>([^<]+)</Size>", block)
        lm = re.search(r"<LastModified>([^<]+)</LastModified>", block)
        if key and key.group(1).endswith(".parquet"):
            rows.append({
                "resource": re.sub(r".*/", "", key.group(1)),
                "size_bytes": int(size.group(1)) if size else None,
                "bucket_modified": lm.group(1) if lm else None,
            })

    info = pd.DataFrame(rows).sort_values("resource").reset_index(drop=True)

    try:
        with urllib.request.urlopen(OSN_PROVENANCE_URL) as resp:
            json_text = resp.read().decode()
        pat = r'"([^"]+\.parquet)"\s*:\s*\{[^}]*"ncbi_last_modified"\s*:\s*"([^"]*)"'
        prov = {m.group(1): m.group(2) for m in re.finditer(pat, json_text)}
        info["ncbi_last_modified"] = info["resource"].map(prov)
    except Exception:
        info["ncbi_last_modified"] = None

    return info
