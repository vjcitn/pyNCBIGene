"""Module-level duckdb connection and BiocFileCache singleton."""

import os
from pathlib import Path

import duckdb

try:
    from pyBiocFileCache import BiocFileCache as _BFC
    _HAS_BFC = True
except ImportError:
    _HAS_BFC = False

_con: duckdb.DuckDBPyConnection | None = None
_cache = None

# Resources whose taxid column is not "#tax_id"
TAXID_COL: dict[str, str] = {
    "gene_refseq_uniprotkb_collab": "NCBI_tax_id",
}

OSN_BASE = (
    "https://mghp.osn.xsede.org/bir190004-bucket01/BiocParquetNCBI/{}.parquet"
)
OSN_LISTING_URL = (
    "https://mghp.osn.xsede.org/bir190004-bucket01?prefix=BiocParquetNCBI/"
)
OSN_PROVENANCE_URL = (
    "https://mghp.osn.xsede.org/bir190004-bucket01/"
    "BiocParquetNCBI/provenance.json"
)


def taxid_column(gres: str) -> str:
    return TAXID_COL.get(gres, "#tax_id")


def get_connection() -> duckdb.DuckDBPyConnection:
    global _con
    if _con is None or _con.is_closed():
        ext_dir = Path.home() / ".cache" / "pyNCBIGene" / "duckdb_extensions"
        ext_dir.mkdir(parents=True, exist_ok=True)
        _con = duckdb.connect(config={"extension_directory": str(ext_dir)})
        _con.execute("INSTALL httpfs; LOAD httpfs;")
    return _con


def get_cache():
    global _cache
    if _cache is None:
        if not _HAS_BFC:
            raise ImportError(
                "pyBiocFileCache is required for caching. "
                "Install with: pip install pybiocfilecache"
            )
        cache_dir = Path.home() / ".cache" / "pyNCBIGene" / "BiocFileCache"
        cache_dir.mkdir(parents=True, exist_ok=True)
        _cache = _BFC(str(cache_dir))
    return _cache


def set_cache(cache):
    """Override the package-level cache (useful for testing)."""
    global _cache
    _cache = cache
