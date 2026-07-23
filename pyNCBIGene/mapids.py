"""map_ids_ng: push-down identifier mapping analogous to AnnotationDbi mapIds."""

from typing import Optional

from .remote import open_ncbi_gene


def map_ids_ng(
    keys: list,
    keytype: str = "Symbol",
    column: str = "GeneID",
    taxid: int = 9606,
    freeze_tag: Optional[str] = None,
) -> dict:
    """Map identifiers using NCBI Gene parquet resources.

    All filtering is pushed to duckdb before a single ``.df()`` retrieves only
    the matched rows.  For Ensembl keytype or column, a JOIN between
    ``gene2ensembl`` and ``gene_info`` runs entirely in duckdb.

    Parameters
    ----------
    keys : list
        Identifiers to translate.
    keytype : str
        One of ``"Symbol"``, ``"GeneID"``, or ``"Ensembl"``.
    column : str
        Output annotation column (e.g. ``"GeneID"``, ``"map_location"``,
        ``"Ensembl"``).
    taxid : int
        NCBI taxonomy ID.
    freeze_tag : str or None
        When set, uses a frozen snapshot (see :func:`freeze_taxon_cache`).

    Returns
    -------
    dict
        ``{key: value}`` mapping; value is ``None`` for unrecognised keys.
    """
    supported = ("Symbol", "GeneID", "Ensembl")
    if keytype not in supported:
        raise ValueError(f"{keytype} keytype not supported. Use one of {supported}.")

    keys_sql = ", ".join(f"'{k}'" for k in keys)

    def _info(cols):
        return (
            open_ncbi_gene("gene_info", taxid, freeze_tag=freeze_tag)
            .select(cols)
        )

    def _g2e():
        return open_ncbi_gene("gene2ensembl", taxid, freeze_tag=freeze_tag)

    if keytype == "Ensembl":
        g2e = _g2e()
        g2e_f = g2e.filter(
            f'"Ensembl_gene_identifier" IN ({keys_sql})'
        ).select('"GeneID", "Ensembl_gene_identifier"')

        if column == "GeneID":
            df = g2e_f.df().rename(columns={"Ensembl_gene_identifier": "Ensembl"})
        else:
            info = _info(f'"GeneID", "{column}"')
            from ._state import get_connection
            con = get_connection()
            con.register("_g2e_tmp", g2e_f.df())
            con.register("_info_tmp", info.df())
            df = con.sql(
                f'SELECT g."Ensembl_gene_identifier" AS "Ensembl", i."{column}" '
                'FROM _g2e_tmp g LEFT JOIN _info_tmp i USING ("GeneID")'
            ).df()
        result = dict(zip(df["Ensembl"], df[column if column != "Ensembl" else "Ensembl"]))

    elif column == "Ensembl":
        info = _info(f'"GeneID", "{keytype}"').filter(
            f'"{keytype}" IN ({keys_sql})'
        )
        g2e = _g2e().select('"GeneID", "Ensembl_gene_identifier"')
        from ._state import get_connection
        con = get_connection()
        con.register("_info_tmp", info.df())
        con.register("_g2e_tmp", g2e.df())
        df = con.sql(
            f'SELECT i."{keytype}", g."Ensembl_gene_identifier" AS "Ensembl" '
            'FROM _info_tmp i LEFT JOIN _g2e_tmp g USING ("GeneID")'
        ).df()
        result = dict(zip(df[keytype], df["Ensembl"]))

    else:
        df = (
            open_ncbi_gene("gene_info", taxid, freeze_tag=freeze_tag)
            .filter(f'"{keytype}" IN ({keys_sql})')
            .select(f'"{keytype}", "{column}"')
            .df()
        )
        # keep first match per key
        result = {}
        for _, row in df.iterrows():
            k = row[keytype]
            if k not in result:
                result[k] = row[column]

    return {k: result.get(k) for k in keys}
