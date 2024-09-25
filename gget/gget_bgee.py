import pandas as pd
import json as json_

import requests

from .utils import set_up_logger, json_list_to_df

logger = set_up_logger()


def _bgee_species(gene_id: str, verbose=True):
    """
    Get species ID from Bgee
    :param gene_id: Ensembl gene ID
    :param verbose: log progress
    :return: species ID
    """

    if verbose:
        logger.info(f"Getting species ID for gene {gene_id} from Bgee")

    response = requests.get(
        "https://bgee.org/api/",
        params={
            "display_type": "json",
            "page": "gene",
            "action": "general_info",
            "gene_id": gene_id,
        },
    )

    if not response.ok:
        raise RuntimeError(
            f"Bgee API request returned with error code: {response.status_code}. "
            "Please double-check the arguments and try again.\n"
        )

    genes_data = response.json()["data"]["genes"]
    assert len(genes_data) == 1
    gene_data = genes_data[0]

    species: int = gene_data["species"]["genomeSpeciesId"]
    return species


def _bgee_orthologs(gene_id, json=False, verbose=True):
    """
    Get orthologs for a gene from Bgee

    Args:

    :param gene_id: Ensembl gene ID
    :param json:    return JSON instead of DataFrame
    :param verbose: log progress

    Returns requested information as a DataFrame or JSON
    """

    # must first obtain species
    species = _bgee_species(gene_id, verbose=verbose)

    if verbose:
        logger.info(f"Getting orthologs for gene {gene_id} from Bgee")

    # then obtain homologs
    response = requests.get(
        f"https://bgee.org/api/",
        params={
            "display_type": "json",
            "page": "gene",
            "action": "homologs",
            "gene_id": gene_id,
            "species_id": species,
        },
    )

    if not response.ok:
        raise RuntimeError(
            f"Bgee API request returned with error code: {response.status_code}. "
            "Please double-check the arguments and try again.\n"
        )

    homologs_data = response.json()["data"]["orthologsByTaxon"]
    homologs_data = sum([v["genes"] for v in homologs_data], [])

    df = json_list_to_df(
        homologs_data,
        [
            ("gene_id", "geneId"),
            ("gene_name", "name"),
            ("species_id", "species.id"),
            ("genus", "species.genus"),
            ("species", "species.speciesName"),
        ],
    )

    if json:
        return json_.loads(df.to_json(orient="records", force_ascii=False))
    else:
        return df


def _bgee_expression(gene_id, json=False, verbose=True):
    """
    Get expression data from Bgee

    Args:

    :param gene_id: Ensembl gene ID
    :param json:    return JSON instead of DataFrame
    :param verbose: log progress

    Returns requested information as a DataFrame or JSON
    """
    # must first obtain species
    species = _bgee_species(gene_id, verbose=verbose)

    if verbose:
        logger.info(f"Getting expression data for gene {gene_id} from Bgee")

    # then obtain expression data
    response = requests.get(
        "https://bgee.org/api/",
        params={
            "display_type": "json",
            "page": "gene",
            "action": "expression",
            "gene_id": gene_id,
            "species_id": species,
            "cond_param": ["anat_entity", "cell_type"],
            "data_type": "all",
        },
    )

    if not response.ok:
        raise RuntimeError(
            f"Bgee API request returned with error code: {response.status_code}. "
            "Please double-check the arguments and try again.\n"
        )

    expression_data = response.json()["data"]["calls"]

    df = json_list_to_df(
        expression_data,
        [
            ("anat_entity_id", "condition.anatEntity.id"),
            ("anat_entity_name", "condition.anatEntity.name"),
            ("score", "expressionScore.expressionScore"),
            ("score_confidence", "expressionScore.expressionScoreConfidence"),
            ("expression_state", "expressionState"),
        ],
    )
    df["score"] = df["score"].astype(float)

    if json:
        return json_.loads(df.to_json(orient="records", force_ascii=False))
    else:
        return df


# noinspection PyShadowingBuiltins
def bgee(
    gene_id,
    type,
    json=False,
    verbose=True,
):
    """
    Get orthologs/expression data for a gene from Bgee (https://www.bgee.org/).

    Args:
    type        type of data to retrieve (expression or orthologs)
    gene_id     Ensembl gene ID
    json        return JSON instead of DataFrame
    verbose     log progress

    Returns requested information as a DataFrame or JSON.
    """
    if type == "expression":
        return _bgee_expression(gene_id, json=json, verbose=verbose)
    elif type == "orthologs":
        return _bgee_orthologs(gene_id, json=json, verbose=verbose)
    else:
        raise ValueError(f"Unknown type: {type}")
