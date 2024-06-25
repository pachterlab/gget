import time
import typing
from typing import Literal

import pandas as pd

from .constants import OPENTARGETS_GRAPHQL_API
from .utils import set_up_logger, wrap_cols_func, graphql_query, json_list_to_df

logger = set_up_logger()


def opentargets(
    ensembl_id: str,
    resource: Literal["diseases", "drugs", "tractability"] = "diseases",
    limit: int | None = None,
    verbose: bool = True,
    wrap_text: bool = False,
) -> pd.DataFrame:
    """
    Query OpenTargets for data associated with a given Ensembl gene ID.

    Args:

    - ensembl_id    Ensembl gene ID to be queried (str), e.g. "ENSG00000169194".
    - resource      Defines type of information to be returned.
                    "diseases":     Returns diseases associated with the gene (default).
                    "drugs":        Returns drugs associated with the gene.
                    "tractability": Returns tractability data for the gene.
    - limit         Limit the number of results returned (default: No limit).
    - verbose       Print progress messages (default: True).
    - wrap_text     If True, displays data frame with wrapped text for easy reading. Default: False.

    Returns requested information in DataFrame format.
    """

    # check if resource argument is valid
    if resource not in _RESOURCES:
        raise ValueError(
            f"'resource' argument specified as {resource}. Expected one of: {', '.join(_RESOURCES)}"
        )

    return _RESOURCES[resource](
        ensembl_id, limit=limit, verbose=verbose, wrap_text=wrap_text
    )


def _limit_pagination() -> tuple[str, str, callable]:
    def f(limit: int | None):
        if limit is None:
            limit = 1
        return {"index": 0, "size": limit}
    return "page", "Pagination", f


def _limit_size() -> tuple[str, str, callable]:
    def f(limit: int | None):
        if limit is None:
            limit = 1
        return limit
    return "size", "Int", f


def _limit_not_supported() -> tuple[None, None, callable]:
    def f(limit: int):
        if limit is not None:
            raise ValueError("Limit is not supported for this resource.")
        return 1

    return None, None, f


def _tractability_converter(row: dict[str, ...]):
    _modality_map = {
        "SM": "Small molecule",
        "AB": "Antibody",
        "PR": "PROTAC",
        "OC": "Other",
    }
    row["modality"] = _modality_map.get(row["modality"], row["modality"])


def _make_query_fun(
    top_level_key: str,
    row_query: str,
    human_readable_tlk: str,
    df_schema: list[tuple[str, str]],
    wrap_columns: list[str],
    limit_func: callable,
    filter_: typing.Callable[[dict[str, ...]], bool] = lambda x: True,
    converter: typing.Callable[[dict[str, ...]], dict[str, ...]] = lambda x: None,
) -> callable:
    """
    Make a query function for OpenTargets API.

    Args:

    - top_level_key         Top level key in the GraphQL query response, e.g. "associatedDiseases".
    - row_query             Query string for the row data, e.g. "score disease{id name description}".
    - human_readable_tlk    Human readable version of the top level key, e.g. "associated diseases".
    - df_schema             Schema for the DataFrame, e.g. [("id", "disease.id"), ("name", "disease.name")].
    - wrap_columns          Columns to wrap text for easy reading, e.g. ["description"].
    - limit_func            Function to set the limit in the pagination object.
    - filter_               Function to filter the raw data (return False to skip value).
                            Note: applied after limit but before converter.
    - converter             Function to optionally manipulate the raw data IN PLACE before it is converted to a DataFrame.

    Returns a function that queries the OpenTargets API and returns the data in DataFrame format.
    """

    limit_key, limit_type, limit_func = limit_func()

    def fun(
        ensembl_id: str,
        limit: int | None = None,
        verbose: bool = True,
        wrap_text: bool = False,
    ) -> pd.DataFrame:
        if limit_key is None:
            query_string = """
            query target($ensemblId: String!) {
                target(ensemblId: $ensemblId) {
                    <TOP_LEVEL_KEY>{
                        <ROW_QUERY>
                    }
                }
            }
            """.replace(
                "<TOP_LEVEL_KEY>", top_level_key
            ).replace(
                "<ROW_QUERY>", row_query
            )
        else:
            query_string = """
            query target($ensemblId: String!, $pagination: <LIMIT_TYPE>) {
                target(ensemblId: $ensemblId) {
                    <TOP_LEVEL_KEY>(<LIMIT_KEY>: $pagination){
                        count
                        rows{
                            <ROW_QUERY>
                        }
                    }
                }
            }
            """.replace(
                "<TOP_LEVEL_KEY>", top_level_key
            ).replace(
                "<LIMIT_TYPE>", limit_type
            ).replace(
                "<LIMIT_KEY>", limit_key
            ).replace(
                "<ROW_QUERY>", row_query
            )

        pagination = limit_func(limit)
        variables = {"ensemblId": ensembl_id, "pagination": pagination}

        if limit_key is None:
            del variables["pagination"]

        results = graphql_query(OPENTARGETS_GRAPHQL_API, query_string, variables)
        target: dict[str, ...] = results["data"]["target"]
        if target is None:
            raise ValueError(
                f"No data found for Ensembl ID: {ensembl_id}. Please double-check the ID and try again."
            )
        data: dict[str, ...] = target[top_level_key]

        if limit_key is None:
            # noinspection PyTypeChecker
            rows: list[dict[str, ...]] = data
            total_count = len(data)
        else:
            total_count: int = data["count"]
            rows: list[dict[str, ...]] = data["rows"]

        if verbose:
            explanation = ""
            if limit is None and limit_key is not None:
                explanation = " (Querying count, will fetch all results next.)"
            logger.info(
                f"Retrieved {len(rows)}/{total_count} {human_readable_tlk}.{explanation}"
            )

        if limit is None and limit_key is not None:
            # wait 1 second as a courtesy
            time.sleep(1)
            variables["pagination"] = limit_func(total_count)

            new_results = graphql_query(
                OPENTARGETS_GRAPHQL_API, query_string, variables
            )
            new_data: dict[str, ...] = new_results["data"]["target"][top_level_key]
            new_rows: list[dict[str, ...]] = new_data["rows"]
            # we re-fetched the original 1, so we need to replace them
            rows = new_rows
            if verbose:
                logger.info(
                    f"Retrieved {len(rows)}/{total_count} {human_readable_tlk}."
                )

        if limit is not None:
            rows = rows[:limit]

        rows = [row for row in rows if filter_(row)]
        for row in rows:
            converter(row)

        df = json_list_to_df(
            rows,
            df_schema,
        )

        if wrap_text:
            df_wrapped = df.copy()
            wrap_cols_func(df_wrapped, wrap_columns)

        return df

    return fun


_RESOURCES = {
    "diseases": _make_query_fun(
        "associatedDiseases",
        """
        score
        disease{
            id
            name
            description
        }""",
        "associated diseases",
        [
            ("id", "disease.id"),
            ("name", "disease.name"),
            ("description", "disease.description"),
            ("score", "score"),
        ],
        ["description"],
        _limit_pagination,
    ),
    "drugs": _make_query_fun(
        "knownDrugs",
        """
        # Basic drug data
        drugId
        prefName
        drugType
        mechanismOfAction
        # Disease data
        disease{
            id
            name
        }
        # Clinical trial data
        phase
        status
        ctIds
        # More drug data
        drug{
            synonyms
            tradeNames
            description
            isApproved
        }
        """,
        "known drugs",
        [
            # Drug data
            ("id", "drugId"),
            ("name", "prefName"),
            ("type", "drugType"),
            ("action_mechanism", "mechanismOfAction"),
            ("description", "drug.description"),
            ("synonyms", "drug.synonyms"),
            ("trade_names", "drug.tradeNames"),
            # Disease data
            ("disease_id", "disease.id"),
            ("disease_name", "disease.name"),
            # Trial data
            ("trial_phase", "phase"),
            ("trial_status", "status"),
            ("trial_ids", "ctIds"),
            ("approved", "drug.isApproved"),
        ],
        ["description", "synonyms", "trade_names", "trial_ids"],
        _limit_size,
    ),
    "tractability": _make_query_fun(
        "tractability",
        """
        label
        modality
        value
        """,
        "tractability states",
        [
            ("label", "label"),
            ("modality", "modality"),
        ],
        [],
        _limit_not_supported,
        filter_=lambda x: x["value"],
        converter=_tractability_converter,
    )
}

OPENTARGETS_RESOURCES = list(_RESOURCES.keys())
