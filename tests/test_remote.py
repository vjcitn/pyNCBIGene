"""Tests for remote query functions -- skipped when offline."""

import pytest
import pandas as pd

try:
    import urllib.request
    urllib.request.urlopen("http://example.com", timeout=3)
    ONLINE = True
except Exception:
    ONLINE = False

skip_offline = pytest.mark.skipif(not ONLINE, reason="no internet connection")

from pyNCBIGene import (
    available_ncbi_parquet,
    ncbi_gene_fields,
    open_ncbi_gene,
    map_ids_ng,
    join_ncbi_gene,
)


@skip_offline
def test_available_ncbi_parquet_returns_list():
    res = available_ncbi_parquet()
    assert isinstance(res, list)
    assert len(res) > 0
    assert all(r.endswith(".parquet") for r in res)


@skip_offline
def test_available_ncbi_parquet_includes_core_resources():
    res = available_ncbi_parquet()
    assert "gene_info.parquet" in res
    assert "gene2go.parquet" in res
    assert "gene2ensembl.parquet" in res


@skip_offline
def test_open_ncbi_gene_rejects_unknown_resource():
    with pytest.raises(ValueError, match="not an available resource"):
        open_ncbi_gene("gene2banana")


@skip_offline
def test_open_ncbi_gene_returns_lazy_relation():
    import duckdb
    rel = open_ncbi_gene("gene_info", taxid=9606)
    assert hasattr(rel, "df"), "Expected a duckdb relation with .df() method"


@skip_offline
def test_open_ncbi_gene_taxid_drops_tax_id_column():
    rel = open_ncbi_gene("gene_info", taxid=9606)
    df = rel.limit(1).df()
    assert "#tax_id" not in df.columns


@skip_offline
def test_open_ncbi_gene_no_taxid_retains_tax_id_column():
    rel = open_ncbi_gene("gene_info")
    df = rel.limit(1).df()
    assert "#tax_id" in df.columns


@skip_offline
def test_ncbi_gene_fields_returns_dataframe():
    fields = ncbi_gene_fields("gene_info")
    assert isinstance(fields, pd.DataFrame)
    assert "column_name" in fields.columns
    assert "column_type" in fields.columns
    assert "Symbol" in fields["column_name"].values
    assert "GeneID" in fields["column_name"].values


@skip_offline
def test_ncbi_gene_fields_rejects_unknown_resource():
    with pytest.raises(ValueError, match="not an available resource"):
        ncbi_gene_fields("gene2banana")


@skip_offline
def test_ncbi_gene_fields_differ_across_resources():
    gi = ncbi_gene_fields("gene_info")["column_name"].tolist()
    g2g = ncbi_gene_fields("gene2go")["column_name"].tolist()
    assert "Symbol" in gi
    assert "Symbol" not in g2g
    assert "GO_ID" in g2g


@skip_offline
def test_map_ids_ng_symbol_to_geneid():
    result = map_ids_ng(["ORMDL3", "TP53", "XyZZY"], keytype="Symbol",
                        column="GeneID", taxid=9606)
    assert result["ORMDL3"] == 94103
    assert result["TP53"] == 7157
    assert result["XyZZY"] is None


@skip_offline
def test_map_ids_ng_rejects_unsupported_keytype():
    with pytest.raises(ValueError, match="keytype not supported"):
        map_ids_ng(["TP53"], keytype="MIM", column="GeneID", taxid=9606)


@skip_offline
def test_map_ids_ng_ensembl_to_symbol():
    result = map_ids_ng(
        ["ENSG00000073605", "ENSG00000141510"],
        keytype="Ensembl", column="Symbol", taxid=9606
    )
    assert result["ENSG00000073605"] == "GSDMB"
    assert result["ENSG00000141510"] == "TP53"


@skip_offline
def test_join_ncbi_gene_left_join():
    local = pd.DataFrame({"Symbol": ["ORMDL3", "TP53", "XyZZY"]})
    result = join_ncbi_gene(local, by="Symbol", resource="gene_info",
                            taxid=9606).df()
    assert "GeneID" in result.columns
    assert len(result) >= 2  # XyZZY has no match, left join keeps the row


@skip_offline
def test_join_ncbi_gene_rejects_bad_resource():
    local = pd.DataFrame({"GeneID": [94103]})
    with pytest.raises(ValueError, match="not an available resource"):
        join_ncbi_gene(local, by="GeneID", resource="gene2banana")


@skip_offline
def test_join_ncbi_gene_rejects_missing_local_column():
    local = pd.DataFrame({"Symbol": ["TP53"]})
    with pytest.raises(ValueError, match="not in local_df"):
        join_ncbi_gene(local, by="GeneID", resource="gene_info", taxid=9606)


@skip_offline
def test_join_ncbi_gene_rejects_missing_remote_column():
    local = pd.DataFrame({"Banana": ["TP53"]})
    with pytest.raises(ValueError, match="not in 'gene_info'"):
        join_ncbi_gene(local, by="Banana", resource="gene_info", taxid=9606)
