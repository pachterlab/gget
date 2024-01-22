import requests
import pandas as pd
import json as json_package
import logging

# Add and format time stamp in logging messages
logging.basicConfig(
    format="%(asctime)s %(levelname)s %(message)s",
    level=logging.INFO,
    datefmt="%c",
)
# Mute numexpr threads info
logging.getLogger("numexpr").setLevel(logging.WARNING)

# Constants
from .constants import COSMIC_GET_URL


def cosmic(
    searchterm, entity="mutations", limit=100, save=False, verbose=True, json=False
):
    """
    Search for genes, mutations, etc associated with cancers using the COSMIC
    (Catalogue Of Somatic Mutations In Cancer) database
    (https://cancer.sanger.ac.uk/cosmic).

    Args:
    - searchterm    (str) Search term, which can be a mutation, or gene name (or Ensembl ID), or sample, etc.
                    as defined using the 'entity' argument. Example: 'EGFR'.
    - entity        (str) Defines the type of the supplied search term. One of the following:
                    'mutations' (default), 'genes', 'cancer', 'tumour site', 'studies', 'pubmed', or 'samples'.
    - limit         (int) Number of hits to return. Default: 100
    - json          (True/False) If True, returns results in json format instead of data frame. Default: False
    - save          (True/False) whether to save the results in the local directory. Default: False
    - verbose       (True/False) whether to print progress information. Default: True

    Returns a data frame with the requested results.
    """

    if verbose:
        logging.info("NOTE: Licence fees apply for the commercial use of COSMIC.")

    # Check if 'entity' argument is valid
    sps = [
        "mutations",
        "pubmed",
        "genes",
        "studies",
        "samples",
        "cancer",
        "tumour site",
    ]
    if entity not in sps:
        raise ValueError(
            f"'entity' argument specified as {entity}. Expected one of: {', '.join(sps)}"
        )

    # Translate categories to match COSMIC data table IDs
    if entity == "cancer":
        entity = "disease"

    if entity == "tumour site":
        entity = "tumour"

    if verbose:
        logging.info(
            f"Fetching the {limit} most correlated to {searchterm} from COSMIC."
        )

    r = requests.get(url=COSMIC_GET_URL + entity + "?q=" + searchterm + "&export=json")

    # Check if the request returned an error (e.g. gene not found)
    if not r.ok:
        raise RuntimeError(
            f"COSMIC API request returned error {r.status_code}. "
            "Please double-check the arguments and try again.\n"
        )

    if r.text == "\n":
        logging.warning(
            f"searchterm = '{searchterm}' did not return any results with entity = '{entity}'. "
            "Please double-check the arguments and try again.\n"
        )
        return None

    data = r.text.split("\n")
    dicts = {}
    counter = 1
    if entity == "mutations":
        dicts = {"Gene": [], "Syntax": [], "Alternate IDs": [], "Canonical": []}
        for i in data:
            if len(i) > 2:
                parsing_mutations = i.split("\t")
                dicts["Gene"].append(parsing_mutations[0])
                dicts["Syntax"].append(parsing_mutations[1])
                dicts["Alternate IDs"].append(
                    parsing_mutations[2].replace('" ', "").replace('"', "")
                )
                dicts["Canonical"].append(parsing_mutations[3])
                counter = counter + 1
                if limit < counter:
                    break

    elif entity == "pubmed":
        dicts = {"Pubmed": [], "Paper title": [], "Author": []}
        for i in data:
            if len(i) > 2:
                parsing_mutations = i.split("\t")
                dicts["Pubmed"].append(parsing_mutations[0])
                dicts["Paper title"].append(
                    parsing_mutations[1].replace('" ', "").replace('"', "").capitalize()
                )
                dicts["Author"].append(
                    parsing_mutations[2].replace('" ', "").replace('"', "")
                )
                counter = counter + 1
                if limit < counter:
                    break

    elif entity == "genes":
        dicts = {
            "Gene": [],
            "Alternate IDs": [],
            "Tested samples": [],
            "Simple Mutations": [],
            "Fusions": [],
            "Coding Mutations": [],
        }
        for i in data:
            if len(i) > 2:
                parsing_mutations = i.split("\t")
                dicts["Gene"].append(parsing_mutations[0])
                dicts["Alternate IDs"].append(
                    parsing_mutations[1].replace('" ', "").replace('"', "")
                )
                dicts["Tested samples"].append(parsing_mutations[2])
                dicts["Simple Mutations"].append(parsing_mutations[3])
                dicts["Fusions"].append(parsing_mutations[4])
                dicts["Coding Mutations"].append(parsing_mutations[5])
                counter = counter + 1
                if limit < counter:
                    break

    elif entity == "samples":
        dicts = {
            "Sample Name": [],
            "Sites & Histologies": [],
            "Analysed Genes": [],
            "Mutations": [],
            "Fusions": [],
            "Structual variants": [],
        }
        for i in data:
            if len(i) > 2:
                parsing_mutations = i.split("\t")
                dicts["Sample Name"].append(parsing_mutations[0])
                dicts["Sites & Histologies"].append(
                    parsing_mutations[1].replace(":", ", ")
                )
                dicts["Analysed Genes"].append(parsing_mutations[2])
                dicts["Mutations"].append(parsing_mutations[3])
                dicts["Fusions"].append(parsing_mutations[4])
                dicts["Structual variants"].append(parsing_mutations[5])
                counter = counter + 1
                if limit < counter:
                    break

    elif entity == "studies":
        dicts = {
            "Study Id": [],
            "Project Code": [],
            "Description": [],
        }
        for i in data:
            if len(i) > 2:
                parsing_mutations = i.split("\t")
                dicts["Study Id"].append(parsing_mutations[0])
                dicts["Project Code"].append(parsing_mutations[1])
                dicts["Description"].append(parsing_mutations[2])
                counter = counter + 1
                if limit < counter:
                    break

    elif entity == "disease":
        dicts = {
            "COSMIC classification": [],
            "Paper description": [],
            "Tested samples": [],
            "Mutations": [],
        }
        for i in data:
            if len(i) > 2:
                parsing_mutations = i.split("\t")
                dicts["COSMIC classification"].append(
                    parsing_mutations[0].replace('" ', "").replace('"', "")
                )
                dicts["Paper description"].append(
                    parsing_mutations[1].replace('" ', "").replace('"', "")
                )
                dicts["Tested samples"].append(parsing_mutations[2])
                dicts["Mutations"].append(parsing_mutations[3])
                counter = counter + 1
                if limit < counter:
                    break

    elif entity == "tumour":
        dicts = {
            "Primary Site": [],
            "Tested sample": [],
            "Analyzed genes": [],
            "Mutations": [],
            "Fusions": [],
            "Structural variants": [],
        }
        for i in data:
            if len(i) > 2:
                parsing_mutations = i.split("\t")
                dicts["Primary Site"].append(parsing_mutations[0])
                dicts["Tested sample"].append(parsing_mutations[1])
                dicts["Analyzed genes"].append(parsing_mutations[2])
                dicts["Mutations"].append(parsing_mutations[3])
                dicts["Fusions"].append(parsing_mutations[4])
                dicts["Structural variants"].append(parsing_mutations[5])
                counter = counter + 1
                if limit < counter:
                    break

    corr_df = pd.DataFrame(dicts)

    if json:
        results_dict = json_package.loads(corr_df.to_json(orient="records"))
        if save:
            with open(
                f"gget_cosmic_{entity}_{searchterm}.json", "w", encoding="utf-8"
            ) as f:
                json_package.dump(results_dict, f, ensure_ascii=False, indent=4)

        return results_dict

    else:
        if save:
            corr_df.to_csv(f"gget_cosmic_{entity}_{searchterm}.csv", index=False)

        return corr_df
