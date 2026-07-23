"""Validate pyNCBIGene against the shared golden-query file in the OSN bucket."""

import json
import urllib.request

import pytest

OSN_GOLDEN_URL = (
    "https://mghp.osn.xsede.org/bir190004-bucket01/"
    "BiocParquetNCBI/test_queries.json"
)

try:
    urllib.request.urlopen("http://example.com", timeout=3)
    ONLINE = True
except Exception:
    ONLINE = False

skip_offline = pytest.mark.skipif(not ONLINE, reason="no internet connection")


def _fetch_queries() -> list:
    try:
        with urllib.request.urlopen(OSN_GOLDEN_URL, timeout=10) as resp:
            data = json.loads(resp.read().decode())
        return data["queries"]
    except Exception as exc:
        pytest.skip(f"could not fetch golden query file: {exc}")


@skip_offline
def test_golden_queries():
    """Run every query in test_queries.json and compare to expected output."""
    from pyNCBIGene import available_ncbi_parquet, map_ids_ng

    queries = _fetch_queries()

    for q in queries:
        qid = q["id"]
        fn = q["function"]
        args = q["args"]
        expected = q["expected"]

        if fn == "available_ncbi_parquet":
            result = available_ncbi_parquet()
            assert len(result) >= expected["min_count"], (
                f"golden[{qid}]: expected >= {expected['min_count']} resources, "
                f"got {len(result)}"
            )

        elif fn == "map_ids_ng":
            result = map_ids_ng(
                keys=args["keys"],
                keytype=args["keytype"],
                column=args["column"],
                taxid=args["taxid"],
            )
            for key, exp_val in expected.items():
                got = result.get(key)
                # JSON null -> Python None; GeneIDs come back as int from duckdb
                if exp_val is None:
                    assert got is None, (
                        f"golden[{qid}] key={key}: expected None, got {got!r}"
                    )
                elif isinstance(exp_val, int):
                    assert got == exp_val or int(got) == exp_val, (
                        f"golden[{qid}] key={key}: expected {exp_val}, got {got!r}"
                    )
                else:
                    assert got == exp_val, (
                        f"golden[{qid}] key={key}: expected {exp_val!r}, got {got!r}"
                    )

        else:
            pytest.skip(f"golden[{qid}]: unknown function '{fn}'")
