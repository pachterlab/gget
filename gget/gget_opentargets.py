import json as json_
import pandas as pd

from .constants import OPENTARGETS_GRAPHQL_API
from .utils import set_up_logger, wrap_cols_func, graphql_query, json_list_to_df

logger = set_up_logger()


def opentargets(
    ensembl_id,
    resource="diseases",
    limit=None,
    verbose=True,
    wrap_text=False,
    filters=None,
    filter_mode="and",
    json=False,
):
    """
    Query OpenTargets for data associated with a given Ensembl gene ID.

    Args:

    - ensembl_id    Ensembl gene ID to be queried (str), e.g. "ENSG00000169194".
    - resource      Defines type of information to be returned.
                    "diseases":         Returns diseases associated with the gene (default).
                    "drugs":            Returns drugs associated with the gene.
                    "tractability":     Returns tractability data for the gene.
                    "pharmacogenetics": Returns pharmacogenetics data for the gene.
                    "expression":       Returns gene expression data (by tissues, organs, and anatomical systems).
                    "depmap":           Returns DepMap gene-disease effect data for the gene.
                    "interactions":     Returns protein-protein interactions for the gene.
    - limit         Limit the number of results returned (default: No limit).
                    Note: Not compatible with the 'tractability' and 'depmap' resources.
    - verbose       Print progress messages (default: True).
    - wrap_text     If True, displays data frame with wrapped text for easy reading. Default: False.
    - filters       Filters to apply to the data. Supported filters by resource:
                    "diseases": None
                    "drugs": disease_id (e.g. "EFO_0000274")
                    "tractability": None
                    "pharmacogenetics": drug_id (e.g. "CHEMBL535")
                    "expression": tissue_id (e.g. "UBERON_0002245"), anatomical_system (e.g. "nervous system"), organ (e.g. "brain")
                    "depmap": tissue_id (e.g. "UBERON_0002245")
                    "interactions": protein_a_id (e.g. "ENSP00000304915"), protein_b_id (e.g. "ENSP00000379111"), gene_b_id (e.g. "ENSG00000077238")
    - filter_mode   For resources that support multiple types of filters, this argument specifies how to combine them.
    - json          If True, returns results in JSON format instead of as a Data Frame. Default: False.


    Returns requested information in DataFrame format.
    """

    # check if resource argument is valid
    if resource not in _RESOURCES:
        raise ValueError(
            f"'resource' argument specified as {resource}. Expected one of: {', '.join(_RESOURCES)}"
        )

    # Wrap everything into a list
    if filters is not None:
        filters = {k: v if isinstance(v, list) else [v] for k, v in filters.items()}

    df: pd.DataFrame = _RESOURCES[resource](
        ensembl_id,
        limit=limit,
        verbose=verbose,
        wrap_text=wrap_text,
        filters=filters,
        filter_mode=filter_mode,
    )

    if json:
        return json_.loads(df.to_json(orient="records", force_ascii=False))
    else:
        return df


def _limit_pagination():
    """
    Limit is expressed as (page: {"index": 0, "size": limit}).
    """

    def f(limit, is_rows_based_query):
        if limit is None:
            # special case because `None` is used to probe the total count
            if is_rows_based_query:
                limit = 1
            else:
                return None
        return {"index": 0, "size": limit}

    return "page", "Pagination", f


def _limit_size():
    """
    Limit is expressed as (size: limit).
    """

    def f(limit, is_rows_based_query):
        # special case because `None` is used to probe the total count
        if limit is None and is_rows_based_query:
            limit = 1
        return limit

    return "size", "Int", f


def _limit_not_supported():
    """
    Limit is not supported for this resource (it has no GraphQL-support, and it is meaningless).
    """

    def f(limit, _is_rows_based_query):
        if limit is not None:
            raise ValueError("Limit is not supported for this resource.")
        return None

    return None, None, f


def _limit_deferred():
    """
    Limit is handled after fetching the data (it is not supported by the GraphQL query, but does have meaning).
    """

    def f(_limit, _is_rows_based_query):
        return None

    return None, None, f


def _tractability_converter(row):
    _modality_map = {
        "SM": "Small molecule",
        "AB": "Antibody",
        "PR": "PROTAC",
        "OC": "Other",
    }
    row["modality"] = _modality_map.get(row["modality"], row["modality"])


def _pharmacogenetics_converter(row):
    # need to modify the drugs field to parse it into a dataframe
    drugs = row["drugs"]
    drugs_df = json_list_to_df(drugs, [("id", "drugId"), ("name", "drugFromSource")])
    row["drugs"] = drugs_df


def _mk_list_converter(f):
    def fun(rows):
        for row in rows:
            f(row)

    return fun


def _flatten_depmap(json_entries):
    old_entries = json_entries.copy()
    json_entries.clear()

    for entry in old_entries:
        for screen in entry["screens"]:
            new_entry = entry.copy()
            del new_entry["screens"]
            new_entry["screen"] = screen
            json_entries.append(new_entry)


def _mk_filter_applicator(id_key):
    """
    Make a filter applicator function based on the provided mapping.
    :param id_key: Mapping of filter key to the key in the row data. (e.g. {"disease_id": "disease.id"})
    """

    def f(row, mode, filters):
        for filter_id, filter_values in filters.items():
            split_key = id_key[filter_id].split(".")
            actual_value = row
            for k in split_key:
                if actual_value is None:
                    break
                if type(actual_value) is list:
                    actual_value = [v[k] for v in actual_value]
                else:
                    actual_value = actual_value[k]

            if mode == "and":
                if type(actual_value) is list:
                    if not any(
                        v in filter_values for v in actual_value
                    ):  # if none match, return False
                        return False
                elif actual_value not in filter_values:
                    return False
            else:
                if type(actual_value) is list:
                    if any(
                        v in filter_values for v in actual_value
                    ):  # if any match, return True
                        return True
                elif actual_value in filter_values:
                    return True

        return mode == "and"  # fallthrough return for "and" is True, for "or" is False

    return set(id_key.keys()), f


def _make_query_fun(
    top_level_key,
    inner_query,
    human_readable_tlk,
    df_schema,
    wrap_columns,
    limit_func,
    is_rows_based_query=True,
    sorter=None,
    sort_reverse=False,
    filter_=lambda x: True,
    converter=lambda x: None,
    user_filter=None,
):
    """
    Make a query function for OpenTargets API.

    Args:

    - top_level_key         Top level key in the GraphQL query response, e.g. "associatedDiseases".
    - inner_query           Query string for the row data, e.g. "score disease{id name description}".
    - human_readable_tlk    Human readable version of the top level key, e.g. "associated diseases".
    - df_schema             Schema for the DataFrame, e.g. [("id", "disease.id"), ("name", "disease.name")].
    - wrap_columns          Columns to wrap text for easy reading, e.g. ["description"].
    - limit_func            Function to convert a limit into a pagination variable.
    - is_rows_based_query   If True, the query is wrapped inside `count rows{query}`.
    - sorter                Function to sort the raw data before it is limited, filtered, or converted.
                            Note: you may want to use `_limit_deferred`, otherwise the limit will be
                            applied through GraphQL before sorting.
    - sort_reverse          If True, sort in reverse order.
    - filter_               Function to filter the raw data (return False to skip value).
                            Note: applied after limit but before converter.
    - converter             Function to optionally manipulate the raw data IN PLACE before it is converted to a DataFrame.
    - user_filter           Tuple of a set of valid filter keys and a function to apply filters to the raw data.
                            Run before limit, `filter_`, and `converter`. Only called if filters are provided.

    Returns a function that queries the OpenTargets API and returns the data in DataFrame format.
    """

    limit_key, limit_type, limit_func = limit_func()

    def fun(
        ensembl_id,
        limit=None,
        verbose=True,
        wrap_text=False,
        filters=None,
        filter_mode="and",
    ):
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
                "<ROW_QUERY>", inner_query
            )
        else:
            if is_rows_based_query:
                query_tmp = """
                count
                rows{
                    <INNER_QUERY>
                }""".replace(
                    "<INNER_QUERY>", inner_query
                )
            else:
                query_tmp = inner_query
            query_string = (
                """
            query target($ensemblId: String!, $pagination: <LIMIT_TYPE>) {
                target(ensemblId: $ensemblId) {
                    <TOP_LEVEL_KEY>(<LIMIT_KEY>: $pagination){
                        <QUERY>
                    }
                }
            }
            """.replace(
                    "<TOP_LEVEL_KEY>", top_level_key
                )
                .replace("<LIMIT_TYPE>", limit_type)
                .replace("<LIMIT_KEY>", limit_key)
                .replace("<QUERY>", query_tmp)
            )

        actual_limit = limit

        if filters is not None and len(filters) > 0:
            if user_filter is None:
                raise ValueError("Filters are not supported for this resource.")
            valid_keys, filter_func = user_filter
            invalid_keys = set(filters.keys()) - valid_keys
            if len(invalid_keys) > 0:
                invalid_keys = [f"'{k}'" for k in invalid_keys]
                valid_keys = [f"'{k}'" for k in valid_keys]
                raise ValueError(
                    f"The following filter keys are invalid for this resource: {', '.join(invalid_keys)}. Valid keys are: {', '.join(valid_keys)}"
                )
            # we have to fetch all data to be able to apply limit in a sane way when filters are specified
            limit = None

        pagination = limit_func(limit, is_rows_based_query)
        variables = {"ensemblId": ensembl_id, "pagination": pagination}

        if limit_key is None:
            del variables["pagination"]

        results = graphql_query(OPENTARGETS_GRAPHQL_API, query_string, variables)
        target = results["data"]["target"]
        if target is None:
            raise ValueError(
                f"No data found for Ensembl ID: {ensembl_id}. Please double-check the ID and try again."
            )
        data = target[top_level_key]

        if is_rows_based_query:
            total_count = data["count"]
            rows = data["rows"]
        else:
            # noinspection PyTypeChecker
            rows = data
            total_count = len(data)

        if verbose:
            explanation = ""
            if limit is None and limit_key is not None and is_rows_based_query:
                explanation = " (Querying count, will fetch all results next.)"
            logger.info(
                f"Retrieved {len(rows)}/{total_count} {human_readable_tlk}.{explanation}"
            )

        if limit is None and limit_key is not None and is_rows_based_query:
            variables["pagination"] = limit_func(total_count, is_rows_based_query)

            new_results = graphql_query(
                OPENTARGETS_GRAPHQL_API, query_string, variables
            )
            new_data = new_results["data"]["target"][top_level_key]
            new_rows = new_data["rows"]
            # we re-fetched the original 1, so we need to replace them
            rows = new_rows
            if verbose:
                logger.info(
                    f"Retrieved {len(rows)}/{total_count} {human_readable_tlk}."
                )

        if sorter is not None:
            rows.sort(key=sorter, reverse=sort_reverse)

        if filters is not None and len(filters) > 0:
            rows = [row for row in rows if user_filter[1](row, filter_mode, filters)]

        if actual_limit is not None:
            rows = rows[:actual_limit]

        rows = [row for row in rows if filter_(row)]
        converter(rows)

        if actual_limit is not None:  # just in case the converter changed the length
            rows = rows[:actual_limit]

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
        user_filter=_mk_filter_applicator({"disease_id": "disease.id"}),
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
        converter=_mk_list_converter(_tractability_converter),
        is_rows_based_query=False,
    ),
    "pharmacogenetics": _make_query_fun(
        "pharmacogenomics",
        """
        variantRsId
        
        genotypeId
        genotype
        
        variantFunctionalConsequence{
            id
            label
        }
        
        drugs{
            drugId
            drugFromSource
        }
        phenotypeText
        genotypeAnnotationText
        
        pgxCategory
        isDirectTarget
        
        evidenceLevel
        
        datasourceId
        literature
        """,
        "pharmacogenetic responses",
        [
            ("rs_id", "variantRsId"),
            ("genotype_id", "genotypeId"),
            ("genotype", "genotype"),
            ("variant_consequence_id", "variantFunctionalConsequence.id"),
            ("variant_consequence_label", "variantFunctionalConsequence.label"),
            ("drugs", "drugs"),  # this is processed into a DataFrame by the converter
            ("phenotype", "phenotypeText"),
            ("genotype_annotation", "genotypeAnnotationText"),
            ("response_category", "pgxCategory"),
            ("direct_target", "isDirectTarget"),
            ("evidence_level", "evidenceLevel"),
            ("source", "datasourceId"),
            ("literature", "literature"),
        ],
        ["phenotype", "genotype_annotation"],
        _limit_pagination,
        is_rows_based_query=False,
        converter=_mk_list_converter(_pharmacogenetics_converter),
        user_filter=_mk_filter_applicator({"drug_id": "drugs.drugId"}),
    ),
    "expression": _make_query_fun(
        "expressions",
        """
        tissue{
            id
            label
            anatomicalSystems
            organs
        }
        rna{
            zscore
            value
            unit
            level
        }
        """,
        "baseline expressions",
        [
            ("tissue_id", "tissue.id"),
            ("tissue_name", "tissue.label"),
            ("rna_zscore", "rna.zscore"),
            ("rna_value", "rna.value"),
            ("rna_unit", "rna.unit"),
            ("rna_level", "rna.level"),
            ("anatomical_systems", "tissue.anatomicalSystems"),
            ("organs", "tissue.organs"),
        ],
        ["tissue_name", "anatomical_systems", "organs"],
        _limit_deferred,
        is_rows_based_query=False,
        sorter=lambda x: (x["rna"]["value"], x["rna"]["zscore"]),
        sort_reverse=True,
        user_filter=_mk_filter_applicator(
            {
                "tissue_id": "tissue.id",
                "anatomical_system": "tissue.anatomicalSystems",
                "organ": "tissue.organs",
            }
        ),
    ),
    "depmap": _make_query_fun(
        "depMapEssentiality",
        """
        tissueId
        tissueName
        screens{
            depmapId
            expression
            geneEffect
            cellLineName
            diseaseCellLineId
            diseaseFromSource
            mutation
        }
        """,
        "DepMap entries",
        [
            ("depmap_id", "screen.depmapId"),
            ("expression", "screen.expression"),
            ("effect", "screen.geneEffect"),
            ("tissue_id", "tissueId"),
            ("tissue_name", "tissueName"),
            ("cell_line_name", "screen.cellLineName"),
            ("disease_cell_line_id", "screen.diseaseCellLineId"),
            ("disease_name", "screen.diseaseFromSource"),
            ("mutation", "screen.mutation"),
        ],
        [],
        _limit_not_supported,
        is_rows_based_query=False,
        converter=_flatten_depmap,
        user_filter=_mk_filter_applicator({"tissue_id": "tissueId"}),
    ),
    "interactions": _make_query_fun(
        "interactions",
        """
        score
        count
        sourceDatabase
        
        intA
        targetA{
            id
            approvedSymbol
        }
        intABiologicalRole
        speciesA{
            taxonId
        }
        
        intB
        targetB{
            id
            approvedSymbol
        }
        intBBiologicalRole
        speciesB{
            taxonId
        }
        """,
        "interactions",
        [
            ("evidence_score", "score"),
            ("evidence_count", "count"),
            ("source_db", "sourceDatabase"),
            ("protein_a_id", "intA"),
            ("gene_a_id", "targetA.id"),
            ("gene_a_symbol", "targetA.approvedSymbol"),
            ("role_a", "intABiologicalRole"),
            ("taxon_a", "speciesA.taxonId"),
            ("protein_b_id", "intB"),
            ("gene_b_id", "targetB.id"),
            ("gene_b_symbol", "targetB.approvedSymbol"),
            ("role_b", "intBBiologicalRole"),
            ("taxon_b", "speciesB.taxonId"),
        ],
        [],
        _limit_pagination,
        user_filter=_mk_filter_applicator(
            {"protein_a_id": "intA", "protein_b_id": "intB", "gene_b_id": "targetB.id"}
        ),
    ),
}

OPENTARGETS_RESOURCES = list(_RESOURCES.keys())
