import requests
import json
import time
import pandas as pd

from .constants import OPENTARGETS_GRAPHQL_API
from .utils import set_up_logger, wrap_cols_func

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
    r = requests.post(
        OPENTARGETS_GRAPHQL_API, json={"query": query_string, "variables": variables}
    )
    if r.status_code != 200:
        raise RuntimeError(
            f"The OpenTargets server responded with status code: {r.status_code}. "
            "Please double-check arguments and try again.\n"
        )

    results: dict[str, ...] = json.loads(r.text)
    data: dict[str, ...] = results["data"]["target"]["associatedDiseases"]

    total_count: int = data["count"]
    rows: list[dict[str, ...]] = data["rows"]

    if verbose:
        explanation = ""
        if limit is None:
            explanation = " (Querying count, will fetch all results next.)"
        logger.info(
            f"Retrieved {len(rows)}/{total_count} associated diseases." + explanation
        )

    if len(rows) < total_count:
        if limit is None:
            # wait 1 second as a courtesy
            time.sleep(1)
            variables["pagination"] = {"index": 0, "size": total_count}

            r = requests.post(
                OPENTARGETS_GRAPHQL_API,
                json={"query": query_string, "variables": variables},
            )
            if r.status_code != 200:
                raise RuntimeError(
                    f"The OpenTargets server responded with status code: {r.status_code}. "
                    "Please double-check arguments and try again.\n"
                )

            new_results: dict[str, ...] = json.loads(r.text)
            new_data: dict[str, ...] = new_results["data"]["target"][
                "associatedDiseases"
            ]
            new_rows: list[dict[str, ...]] = new_data["rows"]
            # we re-fetched the original 1, so we need to replace them
            rows = new_rows
            if verbose:
                logger.info(f"Retrieved {len(rows)}/{total_count} associated diseases.")
        else:
            page_length = len(rows)
            page_index = 1
            # While fewer rows have been returned than requested, keep requesting more
            while len(rows) < total_count and len(rows) < limit:
                # wait 1 second as a courtesy
                time.sleep(1)
                variables["pagination"] = {"index": page_index, "size": page_length}

                r = requests.post(
                    OPENTARGETS_GRAPHQL_API,
                    json={"query": query_string, "variables": variables},
                )
                if r.status_code != 200:
                    raise RuntimeError(
                        f"The OpenTargets server responded with status code: {r.status_code}. "
                        "Please double-check arguments and try again.\n"
                    )

                new_results: dict[str, ...] = json.loads(r.text)
                new_data: dict[str, ...] = new_results["data"]["target"][
                    "associatedDiseases"
                ]
                new_rows: list[dict[str, ...]] = new_data["rows"]
                rows.extend(new_rows)

                if verbose:
                    logger.info(
                        f"Retrieved {len(rows)}/{total_count} associated diseases."
                    )
                page_index += 1

    ids: list[str] = []
    names: list[str] = []
    descriptions: list[str] = []
    scores: list[float] = []

    for row in rows:
        score: float = row["score"]
        disease: dict[str, ...] = row["disease"]
        ids.append(disease["id"])
        names.append(disease["name"])
        descriptions.append(disease["description"])
        scores.append(score)

    df = pd.DataFrame(
        data={"id": ids, "name": names, "description": descriptions, "score": scores}
    )

    if wrap_text:
        df_wrapped = df.copy()
        wrap_cols_func(df_wrapped, ["description"])

    return df
