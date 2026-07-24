"""pyNCBIGene: access NCBI Gene parquet resources in the OSN bucket via duckdb.

Mirrors the functionality of the Bioconductor RNCBIGene R package.
"""

from .remote import (
    available_ncbi_parquet,
    ncbi_gene_fields,
    ncbi_parquet_info,
    open_ncbi_gene,
)
from .mapids import map_ids_ng
from .joins import join_ncbi_gene
from .taxon_cache import (
    cache_by_taxon,
    cached_ncbi_resources,
    clear_taxon_cache,
    freeze_taxon_cache,
    taxon_cache_info,
)

__version__ = "0.1.1"
__all__ = [
    "available_ncbi_parquet",
    "ncbi_gene_fields",
    "ncbi_parquet_info",
    "open_ncbi_gene",
    "map_ids_ng",
    "join_ncbi_gene",
    "cache_by_taxon",
    "cached_ncbi_resources",
    "clear_taxon_cache",
    "freeze_taxon_cache",
    "taxon_cache_info",
]
