# pyNCBIGene

Python access to NCBI Gene annotation stored as Apache Parquet in an NSF Open
Storage Network (OSN) bucket.  Mirrors the Bioconductor
[RNCBIGene](https://github.com/vjcitn/RNCBIGene) R package.

All queries run through a persistent [duckdb](https://duckdb.org) connection
that reads remote parquet directly over HTTPS, pushing filters and column
selections to the parquet layer before any data lands in Python.

## Installation

```
pip install pyNCBIGene
```

Or from source:

```
pip install -e /path/to/pyNCBIGene
```

Optional caching requires [pyBiocFileCache](https://github.com/biocpy/pyBiocFileCache):

```
pip install pybiocfilecache
```

## Quick start

```python
from pyNCBIGene import open_ncbi_gene, map_ids_ng, join_ncbi_gene

# list available resources
from pyNCBIGene import available_ncbi_parquet
available_ncbi_parquet()

# lazy query -- no data fetched yet
tbl = open_ncbi_gene("gene_info", taxid=9606)
tbl.filter('"Symbol" IN (\'TP53\', \'BRCA1\')') \
   .select('"Symbol", "GeneID", "map_location"') \
   .df()

# identifier mapping (push-down to duckdb)
map_ids_ng(["ORMDL3", "TP53", "XyZZY"], keytype="Symbol",
           column="GeneID", taxid=9606)
# {"ORMDL3": 94103, "TP53": 7157, "XyZZY": None}

# Ensembl to Symbol (JOIN in duckdb)
map_ids_ng(["ENSG00000073605", "ENSG00000141510"],
           keytype="Ensembl", column="Symbol", taxid=9606)

# join a local DataFrame to a remote resource
import pandas as pd
local = pd.DataFrame({"GeneID": [94103, 7157, 672]})
join_ncbi_gene(local, by="GeneID", resource="gene2go", taxid=9606).df()
```

## Field discovery

```python
from pyNCBIGene import ncbi_gene_fields
ncbi_gene_fields("gene_info")
ncbi_gene_fields("gene2go")
```

## Non-human organisms

All eight resources cover all organisms in NCBI Gene; just change `taxid`:

```python
# mouse Lilrb4a -- Ensembl gene and transcript IDs
map_ids_ng(["Lilrb4a"], keytype="Symbol", column="Ensembl", taxid=10090)

info = open_ncbi_gene("gene_info", taxid=10090)
g2e  = open_ncbi_gene("gene2ensembl", taxid=10090)
# join in duckdb before collecting
```

## Local caching

```python
from pyNCBIGene import cache_by_taxon, taxon_cache_info, clear_taxon_cache

cache_by_taxon(9606)    # one-time filter + local write (slow for large resources)
taxon_cache_info(9606)  # list cached entries
clear_taxon_cache(9606) # remove live cache, route back to OSN
```

After caching, `open_ncbi_gene("gene_info", taxid=9606)` automatically routes
to the local file with no API change.

## Reproducibility snapshots

```python
from pyNCBIGene import freeze_taxon_cache

freeze_taxon_cache(9606, tag="paper_2026_07")  # physical copy of live cache

# query a frozen snapshot
open_ncbi_gene("gene_info", taxid=9606, freeze_tag="paper_2026_07") \
    .filter('"Symbol" = \'TP53\'').df()
```

Duplicate tags are rejected unless `force=True`.  Frozen entries survive
`clear_taxon_cache()`.

## Bucket metadata and provenance

```python
from pyNCBIGene import ncbi_parquet_info
ncbi_parquet_info()   # size, upload date, NCBI source date per resource
```
