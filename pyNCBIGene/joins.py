"""join_ncbi_gene: join a local DataFrame to a remote NCBI Gene resource."""

from typing import Optional, Union

import pandas as pd

from .remote import _ensure_view, _validate_column, ncbi_gene_fields
from ._state import get_connection, taxid_column


def join_ncbi_gene(
    local_df: pd.DataFrame,
    by: Union[str, list],
    resource: str = "gene_info",
    taxid: Optional[int] = None,
    how: str = "left",
    freeze_tag: Optional[str] = None,
) -> "duckdb.DuckDBPyRelation":
    """Join a local DataFrame against a remote NCBI Gene parquet resource.

    The local DataFrame is registered as an in-memory duckdb table; the remote
    resource is referenced by its existing VIEW name -- no remote data is
    materialised into Python before the join executes.  Column names in ``by``
    are validated against both the local DataFrame and the remote schema before
    any SQL is constructed.

    Parameters
    ----------
    local_df : pd.DataFrame
        Local data to join.
    by : str or list of str
        Column(s) to join on.
    resource : str
        Resource name with or without ``.parquet`` suffix.
    taxid : int or None
        When provided, pre-filters the remote resource to this taxonomy ID.
    how : str
        ``"left"`` (default) or ``"inner"``.
    freeze_tag : str or None
        When set, uses a frozen snapshot.

    Returns
    -------
    duckdb.DuckDBPyRelation
        Lazy join result -- call ``.df()`` to collect.
    """
    if how not in ("left", "inner"):
        raise ValueError("how must be 'left' or 'inner'")

    by_cols = [by] if isinstance(by, str) else list(by)
    gres = resource.replace(".parquet", "")

    # validate join columns against local df
    missing_local = [c for c in by_cols if c not in local_df.columns]
    if missing_local:
        raise ValueError(f"Column(s) not in local_df: {missing_local}")

    # validate join columns against remote schema (also ensures VIEW exists)
    remote_cols = ncbi_gene_fields(resource)["column_name"].tolist()
    missing_remote = [c for c in by_cols if c not in remote_cols]
    if missing_remote:
        raise ValueError(
            f"Column(s) not in '{gres}': {missing_remote}. "
            f"Available: {remote_cols}"
        )

    con = get_connection()

    # use frozen or live local parquet if available, else the remote VIEW
    if freeze_tag is not None:
        from .remote import _frozen_parquet_path
        local_path = _frozen_parquet_path(gres, taxid, freeze_tag)
        vname = "v_frozen_" + gres.replace("-", "_") + f"_{taxid}_{freeze_tag}"
        con.execute(
            f"CREATE OR REPLACE VIEW {vname} AS "
            f"SELECT * FROM read_parquet('{local_path}')"
        )
    else:
        from .remote import _cached_parquet_path
        local_path = _cached_parquet_path(gres, taxid)
        if local_path is not None:
            vname = "v_local_" + gres.replace("-", "_") + f"_{taxid}"
            con.execute(
                f"CREATE OR REPLACE VIEW {vname} AS "
                f"SELECT * FROM read_parquet('{local_path}')"
            )
        else:
            vname = _ensure_view(gres)

    # register the local DataFrame -- zero-copy, stays in-process
    tmp = f"_local_{abs(id(local_df))}"
    con.register(tmp, local_df)

    # join column names are validated above; taxid is int (safe to format)
    on_clause = " AND ".join(f'lhs."{c}" = rhs."{c}"' for c in by_cols)
    tcol = taxid_column(gres)

    if taxid is not None and local_path is None:
        # remote VIEW: add taxid WHERE clause (int, not user string)
        where = f'WHERE rhs."{tcol}" = {int(taxid)}'
    else:
        where = ""

    return con.sql(
        f"SELECT * FROM {tmp} lhs "
        f"{how.upper()} JOIN {vname} rhs ON {on_clause} {where}"
    )
