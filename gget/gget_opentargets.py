import requests
import json
import time
import pandas as pd

from .constants import OPENTARGETS_GRAPHQL_API
from .utils import set_up_logger, wrap_cols_func, graphql_query, json_list_to_df

logger = set_up_logger()


def opentargets(
    ensembl_id: str,
    resource: str = "diseases",
    limit: int | None = None,
    verbose: bool = True,
    wrap_text: bool = False,
) -> pd.DataFrame:
    """
    Query OpenTargets for data associated with a given Ensembl gene ID.

    Args:

    - ensembl_id    Ensembl gene ID to be queried (str), e.g. "ENSG00000169194".
    - resource      Defines type of information to be returned.
                    "diseases": Returns diseases associated with the gene (default).
    - limit         Limit the number of results returned (default: No limit).
    - verbose       Print progress messages (default: True).
    - wrap_text     If True, displays data frame with wrapped text for easy reading. Default: False.

    Returns requested information in DataFrame format.
    """

    # check if resource argument is valid
    resources = {"diseases": _opentargets_disease}
    if resource not in resources:
        raise ValueError(
            f"'resource' argument specified as {resource}. Expected one of: {', '.join(resources)}"
        )

    return resources[resource](
        ensembl_id, limit=limit, verbose=verbose, wrap_text=wrap_text
    )


def _opentargets_disease(
    ensembl_id: str,
    limit: int | None = None,
    verbose: bool = True,
    wrap_text: bool = False,
) -> pd.DataFrame:
    query_string = """
    query target($ensemblId: String!, $pagination: Pagination) {
        target(ensemblId: $ensemblId) {
            associatedDiseases(page: $pagination){
                count
                rows{
                    score
                    disease{
                        id
                        name
                        description
                    }
                }
            }
        }
    }
    """

    if limit is None:
        # when limit is None, we don't fetch any results
        pagination = {"index": 0, "size": 0}
    else:
        pagination = {"index": 0, "size": limit}
    variables = {"ensemblId": ensembl_id, "pagination": pagination}

    results = graphql_query(OPENTARGETS_GRAPHQL_API, query_string, variables)
    target: dict[str, ...] = results["data"]["target"]
    if target is None:
        raise ValueError(
            f"No data found for Ensembl ID: {ensembl_id}. Please double-check the ID and try again."
        )
    data: dict[str, ...] = target["associatedDiseases"]

    total_count: int = data["count"]
    rows: list[dict[str, ...]] = data["rows"]

    if verbose:
        explanation = ""
        if limit is None:
            explanation = " (Querying count, will fetch all results next.)"
        logger.info(
            f"Retrieved {len(rows)}/{total_count} associated diseases." + explanation
        )

    if limit is None:
        # wait 1 second as a courtesy
        time.sleep(1)
        variables["pagination"] = {"index": 0, "size": total_count}

        new_results = graphql_query(OPENTARGETS_GRAPHQL_API, query_string, variables)
        new_data: dict[str, ...] = new_results["data"]["target"]["associatedDiseases"]
        new_rows: list[dict[str, ...]] = new_data["rows"]
        # we re-fetched the original 1, so we need to replace them
        rows = new_rows
        if verbose:
            logger.info(f"Retrieved {len(rows)}/{total_count} associated diseases.")

    df = json_list_to_df(
        rows,
        [
            ("id", "disease.id"),
            ("name", "disease.name"),
            ("description", "disease.description"),
            ("score", "score"),
        ],
    )

    if wrap_text:
        df_wrapped = df.copy()
        wrap_cols_func(df_wrapped, ["description"])

    return df
