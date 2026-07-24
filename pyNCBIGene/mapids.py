"""map_ids_ng: push-down identifier mapping analogous to AnnotationDbi mapIds."""

from typing import Optional

from ._state import get_connection, taxid_column
from .remote import _resolve_source, _validate_column


def map_ids_ng(
    keys: list,
    keytype: str = "Symbol",
    column: str = "GeneID",
    taxid: int = 9606,
    freeze_tag: Optional[str] = None,
) -> dict:
    """Map identifiers using NCBI Gene parquet resources.

    Uses the local BiocFileCache when available, falling back to the remote
    OSN bucket.  Works offline for cached taxa.  All filtering uses
    parameterized duckdb queries.

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
    taxid = int(taxid)

    need_g2e = keytype == "Ensembl" or column == "Ensembl"

    # Resolve views via cache-first logic (no network when data is local).
    gi_vname,  gi_taxid_applied  = _resolve_source("gene_info",    taxid, freeze_tag)
    g2e_vname, g2e_taxid_applied = _resolve_source("gene2ensembl", taxid, freeze_tag) if need_g2e else (None, False)

    gi_tcol  = taxid_column("gene_info")
    g2e_tcol = taxid_column("gene2ensembl")

    # Build WHERE clauses for taxid -- omit if already applied by local view.
    gi_where  = "" if gi_taxid_applied  else f'AND "{gi_tcol}" = {taxid}'
    g2e_where = "" if g2e_taxid_applied else f'AND "{g2e_tcol}" = {taxid}'

    # Validate identifier column names (offline-safe -- uses ncbi_gene_fields
    # which also checks local cache first).
    if keytype != "Ensembl":
        _validate_column(keytype, "gene_info")
    if column not in ("Ensembl",):
        _validate_column(column, "gene_info")

    if keytype == "Ensembl":
        if column == "GeneID":
            rows = con.execute(
                f'SELECT "Ensembl_gene_identifier", "GeneID" '
                f'FROM {g2e_vname} '
                f'WHERE "Ensembl_gene_identifier" IN ({placeholders}) '
                f'{g2e_where}',
                keys,
            ).df()
            result = dict(zip(rows["Ensembl_gene_identifier"], rows["GeneID"]))
        else:
            rows = con.execute(
                f'SELECT g."Ensembl_gene_identifier", i."{column}" '
                f'FROM {g2e_vname} g '
                f'JOIN {gi_vname} i ON g."GeneID" = i."GeneID" '
                f'WHERE g."Ensembl_gene_identifier" IN ({placeholders}) '
                f'{g2e_where} {gi_where}',
                keys,
            ).df()
            result = dict(zip(rows["Ensembl_gene_identifier"], rows[column]))

    elif column == "Ensembl":
        rows = con.execute(
            f'SELECT i."{keytype}", g."Ensembl_gene_identifier" '
            f'FROM {gi_vname} i '
            f'LEFT JOIN {g2e_vname} g ON i."GeneID" = g."GeneID" '
            f'{g2e_where.replace("AND", "AND g.", 1) if g2e_where else ""} '
            f'WHERE i."{keytype}" IN ({placeholders}) '
            f'{gi_where}',
            keys,
        ).df()
        result = dict(zip(rows[keytype], rows["Ensembl_gene_identifier"]))

    else:
        rows = con.execute(
            f'SELECT "{keytype}", "{column}" '
            f'FROM {gi_vname} '
            f'WHERE "{keytype}" IN ({placeholders}) '
            f'{gi_where}',
            keys,
        ).df()
        result = {}
        for _, row in rows.iterrows():
            k = row[keytype]
            if k not in result:
                result[k] = row[column]

    return {k: result.get(k) for k in keys}
