"""Module-level duckdb connection and BiocFileCache singleton.

The connection and cache are process-level singletons, not thread-safe beyond
the lock protecting initial creation.  For multi-threaded use, create
per-thread connections via ``duckdb.connect()`` directly.
"""

import threading
from pathlib import Path

import duckdb

try:
    from pyBiocFileCache import BiocFileCache as _BFC
    _HAS_BFC = True
except ImportError:
    _HAS_BFC = False

_con: duckdb.DuckDBPyConnection | None = None
_con_lock = threading.Lock()
_cache = None
_cache_lock = threading.Lock()

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


def _connection_alive() -> bool:
    try:
        _con.execute("SELECT 1")
        return True
    except Exception:
        return False


def get_connection() -> duckdb.DuckDBPyConnection:
    """Return the process-level duckdb connection, creating it if needed.

    Thread-safe for initial creation.  Long-lived sessions should be aware
    that all callers share this connection.
    """
    global _con
    with _con_lock:
        if _con is None or not _connection_alive():
            ext_dir = Path.home() / ".cache" / "pyNCBIGene" / "duckdb_extensions"
            ext_dir.mkdir(parents=True, exist_ok=True)
            _con = duckdb.connect(config={"extension_directory": str(ext_dir)})
            _con.execute("INSTALL httpfs; LOAD httpfs;")
    return _con


def get_cache():
    """Return the process-level BiocFileCache, creating it if needed."""
    global _cache
    with _cache_lock:
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
    """Override the process-level cache (useful for testing)."""
    global _cache
    _cache = cache
