"""map_ids_ng: push-down identifier mapping analogous to AnnotationDbi mapIds."""

from typing import Optional

from ._state import get_connection, taxid_column
from .remote import _ensure_view, _validate_column, ncbi_gene_fields


def map_ids_ng(
    keys: list,
    keytype: str = "Symbol",
    column: str = "GeneID",
    taxid: int = 9606,
    freeze_tag: Optional[str] = None,
) -> dict:
    """Map identifiers using NCBI Gene parquet resources.

    All filtering uses parameterized duckdb queries -- user-supplied key values
    are never interpolated into SQL strings.  Column names (keytype, column) are
    validated against the resource schema before use.

    Parameters
    ----------
    keys : list
        Identifiers to translate.
    keytype : str
        One of ``"Symbol"``, ``"GeneID"``, or ``"Ensembl"``.
    column : str
        Output annotation column.
    taxid : int
        NCBI taxonomy ID.
    freeze_tag : str or None
        When set, uses a frozen snapshot.

    Returns
    -------
    dict
        ``{key: value}`` mapping; value is ``None`` for unrecognised keys.
    """
    supported = ("Symbol", "GeneID", "Ensembl")
    if keytype not in supported:
        raise ValueError(f"{keytype} keytype not supported. Use one of {supported}.")

    con = get_connection()
    keys = list(keys)
    placeholders = ", ".join(["?"] * len(keys))
    # taxid is always an int -- safe to format directly
    taxid = int(taxid)

    # validate column names against known schema (not user data, but still checked)
    if keytype != "Ensembl":
        _validate_column(keytype, "gene_info")
    if column not in ("Ensembl",):
        _validate_column(column, "gene_info")

    # ensure required VIEWs exist
    gi_vname  = _ensure_view("gene_info")
    g2e_vname = _ensure_view("gene2ensembl") if (keytype == "Ensembl" or column == "Ensembl") else None
    gi_tcol   = taxid_column("gene_info")
    g2e_tcol  = taxid_column("gene2ensembl")

    if keytype == "Ensembl":
        if column == "GeneID":
            # gene2ensembl only -- single parameterized query
            rows = con.execute(
                f'SELECT "Ensembl_gene_identifier", "GeneID" '
                f'FROM {g2e_vname} '
                f'WHERE "{g2e_tcol}" = {taxid} '
                f'AND "Ensembl_gene_identifier" IN ({placeholders})',
                keys,
            ).df()
            result = dict(zip(rows["Ensembl_gene_identifier"], rows["GeneID"]))
        else:
            # JOIN gene2ensembl x gene_info entirely in duckdb, parameterized
            rows = con.execute(
                f'SELECT g."Ensembl_gene_identifier", i."{column}" '
                f'FROM {g2e_vname} g '
                f'JOIN {gi_vname} i ON g."GeneID" = i."GeneID" '
                f'WHERE g."{g2e_tcol}" = {taxid} '
                f'AND i."{gi_tcol}" = {taxid} '
                f'AND g."Ensembl_gene_identifier" IN ({placeholders})',
                keys,
            ).df()
            result = dict(zip(rows["Ensembl_gene_identifier"], rows[column]))

    elif column == "Ensembl":
        # gene_info filter then join gene2ensembl -- parameterized
        rows = con.execute(
            f'SELECT i."{keytype}", g."Ensembl_gene_identifier" '
            f'FROM {gi_vname} i '
            f'LEFT JOIN {g2e_vname} g ON i."GeneID" = g."GeneID" '
            f'  AND g."{g2e_tcol}" = {taxid} '
            f'WHERE i."{gi_tcol}" = {taxid} '
            f'AND i."{keytype}" IN ({placeholders})',
            keys,
        ).df()
        result = dict(zip(rows[keytype], rows["Ensembl_gene_identifier"]))

    else:
        # pure gene_info query -- parameterized
        rows = con.execute(
            f'SELECT "{keytype}", "{column}" '
            f'FROM {gi_vname} '
            f'WHERE "{gi_tcol}" = {taxid} '
            f'AND "{keytype}" IN ({placeholders})',
            keys,
        ).df()
        # keep first match per key
        result = {}
        for _, row in rows.iterrows():
            k = row[keytype]
            if k not in result:
                result[k] = row[column]

    return {k: result.get(k) for k in keys}
