import pandas as pd

import requests

from .utils import set_up_logger, json_list_to_df

logger = set_up_logger()


def _bgee_species(gene_id: str) -> int:
    """
    Get species ID from Bgee
    :param gene_id: Ensembl gene ID
    :return: species ID
    """
    response = requests.get("https://bgee.org/api/", params={
        "display_type": "json",
        "page": "gene",
        "action": "general_info",
        "gene_id": gene_id,
    })

    if not response.ok:
        raise RuntimeError(
            f"Bgee API request returned with error code: {response.status_code}. "
            "Please double-check the arguments and try again.\n"
        )

    genes_data: list[dict[str, ...]] = response.json()['data']['genes']
    assert len(genes_data) == 1
    gene_data: dict[str, ...] = genes_data[0]

    species: int = gene_data['species']['genomeSpeciesId']
    return species


def bgee_orthologs(gene_id: str) -> pd.DataFrame:
    """
    Search for gene in Bgee
    :param gene_id: Ensembl gene ID
    :return: DataFrame with orthologs
    """

    # must first obtain species
    species = _bgee_species(gene_id)

    # then obtain homologs
    response = requests.get(f"https://bgee.org/api/", params={
        "display_type": "json",
        "page": "gene",
        "action": "homologs",
        "gene_id": gene_id,
        "species_id": species,
    })

    if not response.ok:
        raise RuntimeError(
            f"Bgee API request returned with error code: {response.status_code}. "
            "Please double-check the arguments and try again.\n"
        )

    homologs_data: list[dict[str, ...]] = response.json()['data']['orthologsByTaxon']
    homologs_data: list[dict[str, ...]] = sum([v['genes'] for v in homologs_data], [])

    return json_list_to_df(homologs_data, [
        ("gene_id", "geneId"),
        ("gene_name", "name"),
        ("species_id", "species.id"),
        ("genus", "species.genus"),
        ("species", "species.speciesName"),
    ])


def bgee_expression(gene_id: str) -> pd.DataFrame:
    """
    Get expression data from Bgee
    :param gene_id: Ensembl gene ID
    :return: DataFrame with expression data
    """
    species = _bgee_species(gene_id)

    response = requests.get("https://bgee.org/api/", params={
        "display_type": "json",
        "page": "gene",
        "action": "expression",
        "gene_id": gene_id,
        "species_id": species,
        "cond_param": ["anat_entity", "cell_type"],
        "data_type": "all",
    })

    if not response.ok:
        raise RuntimeError(
            f"Bgee API request returned with error code: {response.status_code}. "
            "Please double-check the arguments and try again.\n"
        )

    expression_data: list[dict[str, ...]] = response.json()['data']['calls']

    df =  json_list_to_df(expression_data, [
        ("anat_entity_id", "condition.anatEntity.id"),
        ("anat_entity_name", "condition.anatEntity.name"),
        ("score", "expressionScore.expressionScore"),
        ("score_confidence", "expressionScore.expressionScoreConfidence"),
        ("expression_state", "expressionState"),
    ])
    df["score"] = df["score"].astype(float)
    return df
