import numpy as np
import pandas as pd
import json as json_package
import requests
from bs4 import BeautifulSoup

# import json
import logging

# Add and format time stamp in logging messages
logging.basicConfig(
    format="%(asctime)s %(levelname)s %(message)s",
    level=logging.INFO,
    datefmt="%c",
)
# Custom functions
from .utils import rest_query, get_uniprot_info, wrap_cols_func

# Constants
from .constants import ENSEMBL_REST_API, UNIPROT_REST_API, NCBI_URL

## gget info
def info(ens_ids, expand=False, wrap_text=False, json=False, verbose=True, save=False):
    """
    Fetch gene and transcript metadata using Ensembl IDs.

    Args:
    - ens_ids       One or more Ensembl IDs to look up (string or list of strings).
    - expand        True/False wether to expand returned information (only for genes and transcripts). Default: False.
                    For genes, this adds transcript isoform information.
                    For transcripts, this adds translation and exon information.
    - wrap_text     If True, displays data frame with wrapped text for easy reading. Default: False.
    - json          If True, returns results in json/dictionary format instead of data frame. Default: False.
    - verbose       True/False whether to print progress information. Default True.
    - save          True/False wether to save csv with query results in current working directory. Default: False.

    Returns a data frame containing the requested information about the Ensembl IDs.
    """
    # Define Ensembl REST API server
    server = ENSEMBL_REST_API
    # Define type of returned content from REST
    content_type = "application/json"

    ## Clean up Ensembl IDs
    # If single Ensembl ID passed as string, convert to list
    if type(ens_ids) == str:
        ens_ids = [ens_ids]
    # Remove Ensembl ID version if passed
    ens_ids_clean = []
    temp = 0
    for ensembl_ID in ens_ids:
        ens_ids_clean.append(ensembl_ID.split(".")[0])

        if "." in ensembl_ID and temp == 0:
            if verbose is True:
                logging.info(
                    "We noticed that you may have passed a version number with your Ensembl ID.\n"
                    "Please note that gget info will always return information linked to the latest Ensembl ID version (see 'ensembl_id')."
                )
            temp = +1

    # Initiate dictionary to save results for all IDs in
    master_dict = {}

    # Query REST APIs from https://rest.ensembl.org/
    for ensembl_ID in ens_ids_clean:
        # Create dict to save query results
        results_dict = {ensembl_ID: {}}

        ## lookup/id/ query: Find the species and database for a single identifier
        # Define the REST query
        if expand:
            query = "lookup/id/" + ensembl_ID + "?" + "expand=1"
        else:
            query = "lookup/id/" + ensembl_ID + "?"
        # Submit query
        try:
            df_temp = rest_query(server, query, content_type)

            try:
                # Add Ensembl ID with latest version number to df_temp
                ensembl_id_dict = {
                    "ensembl_id": str(df_temp["id"]) + "." + str(df_temp["version"])
                }
                df_temp.update(ensembl_id_dict)

            except KeyError:
                # Just add Ensembl ID if no version found
                ensembl_id_dict = {"ensembl_id": str(df_temp["id"])}
                df_temp.update(ensembl_id_dict)

        # If query returns in an error:
        except RuntimeError:
            # Try submitting query without expand (expand does not work for exons and translation IDs)
            try:
                query = "lookup/id/" + ensembl_ID + "?"
                df_temp = rest_query(server, query, content_type)
                # Add Ensembl ID with latest version number to df_temp
                ensembl_id_dict = {
                    "ensembl_id": str(df_temp["id"]) + "." + str(df_temp["version"])
                }
                df_temp.update(ensembl_id_dict)

            # Log error if this also did not work
            except RuntimeError:
                logging.error(
                    f"Ensembl ID '{ensembl_ID}' not found. "
                    "Please double-check spelling/arguments and try again."
                )
                continue

        # Add results to master dict
        results_dict[ensembl_ID].update(df_temp)
        master_dict.update(results_dict)

    # Return None if none of the Ensembl IDs were found
    if len(master_dict) == 0:
        return None

    ## Build data frame from dictionary
    df = pd.DataFrame.from_dict(master_dict)

    # Rename indeces
    df = df.rename(
        index={
            "description": "ensembl_description",
            "Parent": "parent_gene",
            "display_name": "ensembl_gene_name",
        }
    )

    # UniProt version 2022_02 requires the ID to be passed including the version number for homo sapiens
    # https://www.uniprot.org/news/2022/05/25/release
    # Collect all IDs with their latest version number
    uniprot_ens_ids = []
    for id_ in ens_ids_clean:
        uniprot_ens_ids.append(master_dict[id_]["ensembl_id"])

    ## For genes and transcripts, get gene names and descriptions from UniProt
    df_temp = pd.DataFrame()
    for ens_id, uniprot_ens_id, id_type in zip(
        ens_ids_clean, uniprot_ens_ids, df.loc["object_type"].values
    ):
        if id_type == "Gene" or id_type == "Transcript":

            # Check if this is a wrombase ID:
            if ens_id.startswith("WB"):
                if id_type == "Gene":
                    df_uniprot = get_uniprot_info(
                        UNIPROT_REST_API, ens_id, id_type="WB_Gene", verbose=verbose
                    )

                else:
                    df_uniprot = get_uniprot_info(
                        UNIPROT_REST_API,
                        ens_id,
                        id_type="WB_Transcript",
                        verbose=verbose,
                    )

            # Check if this is a flybase ID:
            elif ens_id.startswith("FB"):
                df_uniprot = get_uniprot_info(
                    UNIPROT_REST_API, ens_id, id_type="Flybase", verbose=verbose
                )

            # Check if this ID requires a version number
            elif master_dict[ens_id]["species"] == "homo_sapiens":
                df_uniprot = get_uniprot_info(
                    UNIPROT_REST_API, uniprot_ens_id, id_type=id_type, verbose=verbose
                )

            else:
                df_uniprot = get_uniprot_info(
                    UNIPROT_REST_API, ens_id, id_type=id_type, verbose=verbose
                )

            if not isinstance(df_uniprot, type(None)):
                # If two different UniProt IDs for a single query ID are returned, they should be merged into one column
                # So len(df_uniprot) should always be 1
                if len(df_uniprot) > 1:
                    df_uniprot = df_uniprot.iloc[[0]]
                    # This should not be necessary
                    if verbose is True:
                        logging.warning(
                            f"More than one UniProt match was found for Ensembl ID {ens_id}. Only the first match and its associated information will be returned."
                        )

            else:
                if verbose is True:
                    logging.warning(
                        f"No UniProt entry was found for Ensembl ID {ens_id}."
                    )

            ## Get NCBI gene ID and description (for genes only)
            url = NCBI_URL + f"/gene/?term={ens_id}"
            html = requests.get(url)

            # Raise error if status code not "OK" Response
            if html.status_code != 200:
                raise RuntimeError(
                    f"NCBI returned error status code {html.status_code}. Please try again."
                )

            ## Web scrape NCBI website for gene ID, synonyms and description
            soup = BeautifulSoup(html.text, "html.parser")

            # Check if NCBI gene ID is available
            try:
                ncbi_gene_id = soup.find("input", {"id": "gene-id-value"}).get("value")
            except:
                ncbi_gene_id = np.nan

            # Check if NCBI description is available
            try:
                ncbi_description = (
                    soup.find("div", class_="section", id="summaryDiv")
                    .find("dt", text="Summary")
                    .find_next_sibling("dd")
                    .text
                )
            except:
                ncbi_description = np.nan

            # Check if NCBI synonyms are available
            try:
                ncbi_synonyms = (
                    soup.find("div", class_="section", id="summaryDiv")
                    .find("dt", text="Also known as")
                    .find_next_sibling("dd")
                    .text
                )
                # Split NCBI synonyms
                ncbi_synonyms = ncbi_synonyms.split("; ")
            except:
                ncbi_synonyms = None

            # If both NCBI and UniProt synonyms available,
            # final synonyms list will be combined a set of both lists
            if ncbi_synonyms is not None and not isinstance(df_uniprot, type(None)):
                # Collect and flatten UniProt synonyms
                uni_synonyms = df_uniprot["uni_synonyms"].values[0]
                synonyms = list(set().union(uni_synonyms, ncbi_synonyms))
            # Add only UniProt synonyms if NCBI syns not available
            elif ncbi_synonyms is None and not isinstance(df_uniprot, type(None)):
                synonyms = df_uniprot["uni_synonyms"].values[0]
            else:
                synonyms = np.nan

            # Sort synonyms alphabetically (is sortable)
            try:
                synonyms = sorted(synonyms)
            except:
                None

            # Save NCBI info to data frame
            df_ncbi = pd.DataFrame(
                {
                    "ncbi_gene_id": [ncbi_gene_id],
                    "ncbi_description": [ncbi_description],
                    "synonyms": [synonyms],
                },
            )

            # Transpost NCBI df and add Ensembl ID as column name
            df_ncbi = df_ncbi.T
            df_ncbi.columns = [ens_id]

            ## Add NCBI and UniProt info to data frame
            if not isinstance(df_uniprot, type(None)):
                # Transpose UniProt data frame and add Ensembl ID as column name
                df_uniprot = df_uniprot.T
                df_uniprot.columns = [ens_id]

                # Combine Ensembl and NCBI info
                df_uni_ncbi = pd.concat([df_uniprot, df_ncbi])

                # Append NCBI and UniProt info to df_temp
                df_temp = pd.concat([df_temp, df_uni_ncbi], axis=1)

            else:
                # Add only NCBI info to df_temp
                df_temp = pd.concat([df_temp, df_ncbi], axis=1)

    # Append UniProt and NCBI info to df
    df = pd.concat([df, df_temp])

    # Reindex df (this also drops all unmentioned indeces)
    df_final = df.reindex(
        [
            "ensembl_id",
            "uniprot_id",
            "ncbi_gene_id",
            "species",
            "assembly_name",
            "primary_gene_name",
            "ensembl_gene_name",
            "synonyms",
            "parent_gene",
            "protein_names",
            "ensembl_description",
            "uniprot_description",
            "ncbi_description",
            "object_type",
            "biotype",
            "canonical_transcript",
            "seq_region_name",
            "strand",
            "start",
            "end",
        ]
    )

    if expand:
        ens_ids = []
        # Dictionary to save clean info
        data = {
            "all_transcripts": [],
            "transcript_biotypes": [],
            "transcript_names": [],
            "all_exons": [],
            "exon_starts": [],
            "exon_ends": [],
            "all_translations": [],
            "translation_starts": [],
            "translation_ends": [],
        }

        # Collect the transcripts info for each ID in separate lists
        for ens_id in df.columns:
            # Record Ensembl ID
            ens_ids.append(ens_id)

            # Clean up transcript info, if available
            all_transcripts = []
            transcript_biotypes = []
            transcript_names = []
            try:
                for trans_dict in df[ens_id]["Transcript"]:
                    try:
                        all_transcripts.append(trans_dict["id"])
                    except:
                        all_transcripts.append(np.NaN)
                    try:
                        transcript_names.append(trans_dict["display_name"])
                    except:
                        transcript_names.append(np.NaN)
                    try:
                        transcript_biotypes.append(trans_dict["biotype"])
                    except:
                        transcript_biotypes.append(np.NaN)

                data["all_transcripts"].append(all_transcripts)
                data["transcript_biotypes"].append(transcript_biotypes)
                data["transcript_names"].append(transcript_names)

            except:
                data["all_transcripts"].append(np.NaN)
                data["transcript_biotypes"].append(np.NaN)
                data["transcript_names"].append(np.NaN)

            # Clean up exon info, if available
            all_exons = []
            exon_starts = []
            exon_ends = []
            try:
                for exon_dict in df[ens_id]["Exon"]:
                    try:
                        all_exons.append(exon_dict["id"])
                    except:
                        all_exons.append(np.NaN)
                    try:
                        exon_starts.append(exon_dict["start"])
                    except:
                        exon_starts.append(np.NaN)
                    try:
                        exon_ends.append(exon_dict["end"])
                    except:
                        exon_ends.append(np.NaN)

                data["all_exons"].append(all_exons)
                data["exon_starts"].append(exon_starts)
                data["exon_ends"].append(exon_ends)

            except:
                data["all_exons"].append(np.NaN)
                data["exon_starts"].append(np.NaN)
                data["exon_ends"].append(np.NaN)

            # Clean up translation info, if available
            all_translations = []
            translation_starts = []
            translation_ends = []
            try:
                for transl_dict in df[ens_id]["Exon"]:
                    try:
                        all_translations.append(transl_dict["id"])
                    except:
                        all_translations.append(np.NaN)
                    try:
                        translation_starts.append(transl_dict["start"])
                    except:
                        translation_starts.append(np.NaN)
                    try:
                        translation_ends.append(transl_dict["end"])
                    except:
                        translation_ends.append(np.NaN)

                data["all_translations"].append(all_translations)
                data["translation_starts"].append(translation_starts)
                data["translation_ends"].append(translation_ends)

            except:
                data["all_translations"].append(np.NaN)
                data["translation_starts"].append(np.NaN)
                data["translation_ends"].append(np.NaN)

        # Append cleaned up info to df_final
        df_final = pd.concat(
            [df_final, pd.DataFrame.from_dict(data, orient="index", columns=ens_ids)]
        )

    ## Transpose data frame so each row corresponds to one Ensembl ID
    df_final = df_final.T

    # # Add Ensembl ID column from index
    # df_final.insert(0, "ensembl_id", df_final.index)

    if wrap_text:
        df_wrapped = df_final.copy()
        wrap_cols_func(df_wrapped, ["uniprot_description", "ensembl_description"])

    if json:
        results_dict = json_package.loads(df_final.to_json(orient="index"))
        if save:
            with open("gget_info_results.json", "w", encoding="utf-8") as f:
                json_package.dump(results_dict, f, ensure_ascii=False, indent=4)

        return results_dict

    else:
        if save:
            df_final.to_csv("gget_info_results.csv", index=False)

        return df_final
