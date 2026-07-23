"""join_ncbi_gene: join a local DataFrame to a remote NCBI Gene resource."""

from typing import Optional, Union

import pandas as pd

from .remote import ncbi_gene_fields, open_ncbi_gene
from ._state import get_connection


def join_ncbi_gene(
    local_df: pd.DataFrame,
    by: Union[str, list],
    resource: str = "gene_info",
    taxid: Optional[int] = None,
    how: str = "left",
    freeze_tag: Optional[str] = None,
) -> "duckdb.DuckDBPyRelation":
    """Join a local DataFrame against a remote NCBI Gene parquet resource.

    Both sides of the join execute in duckdb.  Call ``.df()`` on the result
    to collect into a pandas DataFrame.

    Parameters
    ----------
    local_df : pd.DataFrame
        Local data to join.
    by : str or list of str
        Column(s) to join on; must exist in both ``local_df`` and the remote
        resource.
    resource : str
        Resource name with or without ``.parquet`` suffix.
    taxid : int or None
        When provided, pre-filters the remote resource to this taxonomy ID.
    how : str
        ``"left"`` (default) or ``"inner"``.
    freeze_tag : str or None
        When set, uses a frozen snapshot (see :func:`freeze_taxon_cache`).

    Returns
    -------
    duckdb.DuckDBPyRelation
        Lazy join result -- call ``.df()`` to collect.
    """
    if how not in ("left", "inner"):
        raise ValueError("how must be 'left' or 'inner'")

    by_cols = [by] if isinstance(by, str) else list(by)
    missing_local = [c for c in by_cols if c not in local_df.columns]
    if missing_local:
        raise ValueError(f"Column(s) not in local_df: {missing_local}")

    remote_cols = ncbi_gene_fields(resource)["column_name"].tolist()
    missing_remote = [c for c in by_cols if c not in remote_cols]
    if missing_remote:
        raise ValueError(
            f"Column(s) not in '{resource}': {missing_remote}. "
            f"Available: {remote_cols}"
        )

    con = get_connection()
    tmp = f"_local_join_{abs(id(local_df))}"
    con.register(tmp, local_df)

    remote_rel = open_ncbi_gene(resource, taxid, freeze_tag=freeze_tag)
    # materialise remote into a temp view so we can join by name
    remote_tmp = f"_remote_join_{abs(id(local_df))}"
    con.register(remote_tmp, remote_rel.df())

    on_clause = " AND ".join(
        f'lhs."{c}" = rhs."{c}"' for c in by_cols
    )
    return con.sql(
        f'SELECT * FROM {tmp} lhs '
        f'{how.upper()} JOIN {remote_tmp} rhs ON {on_clause}'
    )
