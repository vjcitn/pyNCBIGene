# API Changes

Changes to the public API of pyNCBIGene.  When a function is added or its
signature changes, record it here and open a corresponding issue in
[RNCBIGene](https://github.com/vjcitn/RNCBIGene) to keep the R package
in sync.

---

## 0.1.x (current)

### Added
- `open_ncbi_gene(resource, taxid, freeze_tag)` -- lazy `DuckDBPyRelation`
  over any OSN bucket resource; `taxid` drops the taxid column; `freeze_tag`
  opens a frozen snapshot
- `ncbi_gene_fields(resource)` -- column names and types via `DESCRIBE`
- `available_ncbi_parquet()` -- live S3 listing of bucket resources
- `ncbi_parquet_info()` -- sizes, upload dates, NCBI source dates
- `map_ids_ng(keys, keytype, column, taxid, freeze_tag)` -- push-down
  identifier mapping; keytypes: Symbol, GeneID, Ensembl; returns dict with
  None for unrecognised keys
- `join_ncbi_gene(local_df, by, resource, taxid, how, freeze_tag)` -- join
  local DataFrame to remote resource; validates columns before SQL; returns
  lazy relation
- `cache_by_taxon(taxid, resources, force, verbose)` -- filter remote parquet
  to taxon and store in BiocFileCache; subsequent `open_ncbi_gene` calls route
  to local file transparently
- `taxon_cache_info(taxid)` -- list cached entries
- `clear_taxon_cache(taxid)` -- remove live cache entries (frozen entries
  are protected)
- `freeze_taxon_cache(taxid, tag, resources, force)` -- physical snapshot of
  live cache for reproducibility; duplicate tags rejected without `force=True`

### Behaviour notes
- `open_ncbi_gene()` with `taxid` always drops the taxid column from the
  result
- `gene_refseq_uniprotkb_collab` uses `NCBI_tax_id` as its taxid column;
  all other resources use `#tax_id`
- User-supplied key values are passed as `?` parameters to `con.execute()`;
  column names are validated against schema before SQL construction
