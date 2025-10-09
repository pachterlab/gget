from bs4 import BeautifulSoup
import requests

# from requests.adapters import HTTPAdapter, Retry
# import time
import re
import os
import sys
import uuid
import pandas as pd
import numpy as np
from IPython.display import display, HTML
import logging

# from datetime import datetime

# Mute numexpr threads info
logging.getLogger("numexpr").setLevel(logging.WARNING)

from .constants import (
    ENSEMBL_FTP_URL,
    ENSEMBL_FTP_URL_NV,
    ENS_TO_PDB_API,
    COSMIC_RELEASE_URL,
)


def set_up_logger():
    logging_level_name = os.getenv("GGET_LOGLEVEL", "INFO")
    logging_level = logging.getLevelName(logging_level_name)
    if type(logging_level) != int:  # unknown log level
        logging_level = logging.INFO
    logger = logging.getLogger(__name__)
    logger.setLevel(logging_level)

    if not logger.hasHandlers():
        formatter = logging.Formatter(
            "%(asctime)s - %(levelname)s - %(message)s", "%H:%M:%S"
        )

        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

        # log_dir = "logs"
        # if not os.path.exists(log_dir):
        #     os.makedirs(log_dir)

        # log_file = os.path.join(
        #     log_dir, f"log_file_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
        # )

        # file_handler = logging.FileHandler(log_file)
        # file_handler.setFormatter(formatter)
        # logger.addHandler(file_handler)

    return logger


logger = set_up_logger()


def flatten(xss):
    """
    Function to flatten a list of lists.
    """
    return [x for xs in xss for x in xs]


def print_sys(s):
    """system print
    Args:
        s (str): the string to print
    """
    print(s, flush = True, file = sys.stderr)


def get_latest_cosmic():
    html = requests.get(COSMIC_RELEASE_URL)
    if html.status_code != 200:
        raise RuntimeError(
            f"The COSMIC server returned error status code {html.status_code}. Please try again."
        )
    soup = BeautifulSoup(html.text, "html.parser")

    return int(soup.find("div", class_="news").get("id").split("v")[-1])


def check_file_for_error_message(filepath, filename, download_path):
    with open(filepath, "r", encoding="utf-8") as file:
        content = file.read().strip()

    # Define common error indicators
    error_keywords = [
        "Internal Server Error",
        "Bad Gateway",
        "Service Unavailable",
        "502 Bad Gateway",
        "500 Internal Server Error",
    ]

    # Check if the file contains an error message and raise ValueError + print error message if so
    if any(keyword in content for keyword in error_keywords):
        raise ValueError(
            f"""
            The {filename} downloaded from {download_path} 
            contains an error message instead of valid data.\n
            Error message:\n{content}\n
            Please try again. If the problem persists, please report it here: https://github.com/pachterlab/gget/issues/new?template=issue_report.yml
            """
        )


def read_fasta(fasta):
    """
    Args:
    - fasta     (str) Path to fasta file.

    Returns titles and seqs from fasta file as two list objects.
    """
    titles = []
    seqs = []
    title_last = False
    new_seq = False
    with open(fasta) as fasta_file:
        for i, line in enumerate(fasta_file):
            if i == 0 and line[0] != ">":
                raise ValueError("Expected FASTA file to start with a '>' character. ")
            elif line[0] == ">":
                if title_last:
                    raise ValueError(
                        "FASTA file contains two lines starting with '>' in a row -> missing sequence line. "
                    )

                if new_seq:
                    # Append last recorded sequence
                    seqs.append(new_seq)

                # Append title line to titles list
                titles.append(line.strip().replace(">", ""))
                title_last = True

            else:
                if title_last:
                    # Start recording new sequence if title was last
                    new_seq = line.strip()
                else:
                    new_seq = new_seq + line.strip()
                title_last = False

        # Append last sequence
        seqs.append(new_seq)

    return titles, seqs


def n_colors(nucleotide):
    """
    Returns a string format to print the nucleotide
    with its appropriate background color according to the Clustal Colour Scheme.
    """

    # Raw python background colors
    # References:
    # https://stackabuse.com/how-to-print-colored-text-in-python/
    # https://www.ditig.com/256-colors-cheat-sheet
    raw_colors = {
        "black": 0,
        "red": 9,
        "green": 10,
        "yellow": 11,
        "blue": 12,
        "white": 15,
    }

    # Define nulceotide colors dict
    n_colors_dict = {
        "blue": ["C"],
        "red": ["A"],
        "green": ["T"],
        "yellow": ["G"],
    }

    # Define background color based on which nucleotide was passed
    bkg_color = None
    letter_color = None
    for color in n_colors_dict:
        if nucleotide in n_colors_dict[color]:
            bkg_color = raw_colors[color]
            letter_color = color

    # If the nucleotide does not fall into the defined color categories,
    # make it white (e.g. "-")
    if bkg_color == None:
        bkg_color = raw_colors["white"]

    if letter_color is not None and letter_color in ["blue", "red"]:
        # Define textcolor as white for darker colors
        textcolor = raw_colors["white"]
    else:
        # Define textcolor as black
        textcolor = raw_colors["black"]

    return f"\033[38;5;{textcolor}m\033[48;5;{bkg_color}m{nucleotide}\033[0;0m"


def aa_colors(amino_acid):
    """
    Returns a string format to print the amino acid
    with its appropriate background color according to the Clustal Colour Scheme:
    http://www.jalview.org/help/html/colourSchemes/clustal.html
    """

    # Raw python background colors
    # References:
    # https://stackabuse.com/how-to-print-colored-text-in-python/
    # https://www.ditig.com/256-colors-cheat-sheet
    raw_colors = {
        "black": 0,
        "red": 9,
        "green": 10,
        "yellow": 11,
        "blue": 12,
        "cyan": 14,
        "magenta": 5,
        "pink": 13,
        "orange": 1,  # This is maroon because the system colors don't have orange
        "white": 15,
    }

    # Define amino acid colors dict according to http://www.jalview.org/help/html/colourSchemes/clustal.html
    aa_colors_dict = {
        "blue": ["A", "I", "L", "M", "F", "W", "V", "C"],  # hydrophobic AAs
        "red": ["K", "R"],  # positive charge
        "magenta": ["E", "D"],  # negative charge
        "green": ["N", "Q", "S", "T"],  # polar
        "pink": ["C"],  # cysteines
        "orange": ["G"],  # glycines
        "yellow": ["P"],  # prolines
        "cyan": ["H", "Y"],  # aromatic
    }

    # Define background color based on which amino acid was passed
    bkg_color = None
    letter_color = None
    for color in aa_colors_dict:
        if amino_acid in aa_colors_dict[color]:
            bkg_color = raw_colors[color]
            letter_color = color

    # If the amino acid does not fall into the defined color categories,
    # make it white (e.g. "-")
    if bkg_color == None:
        bkg_color = raw_colors["white"]

    if letter_color is not None and letter_color in [
        "blue",
        "red",
        "magenta",
        "orange",
    ]:
        # Define textcolor as white for darker colors
        textcolor = raw_colors["white"]
    else:
        # Define textcolor as black
        textcolor = raw_colors["black"]

    return f"\033[38;5;{textcolor}m\033[48;5;{bkg_color}m{amino_acid}\033[0;0m"


def get_uniprot_seqs(server, ensembl_ids):
    """
    Retrieve UniProt sequences based on Ensemsbl, WormBase or FlyBase identifiers.

    Args:
    - server        Link to UniProt REST API server.
    - ensembl_ids   One or more Ensembl, WormBase or FlyBase IDs (string or list of strings).

    Returns data frame with UniProt ID, gene name, organism, sequence, sequence length, and query ID.
    """

    # If a single UniProt ID is passed as string, convert to list
    if type(ensembl_ids) == str:
        ensembl_ids = [ensembl_ids]

    # Initiate data frame so empty df will be returned if no matches are found
    master_df = pd.DataFrame()

    for id_ in ensembl_ids:
        # API documentation: https://www.uniprot.org/help/api_queries
        # Submit server request
        r = requests.get(server + id_ + "+AND+reviewed:true")
        if not r.ok:
            logger.error(
                f"UniProt server request returned with error status code: {r.status_code}. Please double-check arguments or try again later."
            )
        # Convert to json
        json = r.json()

        # If no reviewed results were found, try again for unreviewed results
        if not len(json["results"]) > 0:
            # Submit server request
            r = requests.get(server + id_)
            if not r.ok:
                logger.error(
                    f"UniProt server request returned with error status code: {r.status_code}. Please double-check arguments or try again later."
                )
            # Convert to json
            json = r.json()

            # Warn user if unreviewed results were found
            if len(json["results"]) > 0:
                logger.warning(
                    f"No reviewed UniProt results were found for ID {id_}. Returning all unreviewed results."
                )

        if len(json["results"]) > 0:
            # Convert results to data frame
            df = pd.json_normalize(json["results"])

            # Remove non-relevant columns
            df = df[
                [
                    "primaryAccession",
                    "organism.scientificName",
                    "sequence.value",
                    "sequence.length",
                ]
            ]

            # Rename columns
            df.columns = [
                "uniprot_id",
                "organism",
                "sequence",
                "sequence_length",
            ]

            # Add gene name and query columns
            gene_names = []
            for i in np.arange(len(json["results"])):
                try:
                    gene_names.append(
                        json["results"][i]["genes"][0]["geneName"]["value"]
                    )
                except:
                    gene_names.append(np.nan)
            df["gene_name"] = gene_names
            df["query"] = id_

            # Append results for this ID to master data frame
            master_df = pd.concat([master_df, df], axis=0)

        else:
            # If no results were found, warn user and do nothing -> returns empty df
            logger.warning(f"No UniProt sequences were found for ID {id_}.")

    return master_df


def get_uniprot_info(server, ensembl_id, verbose=True):
    """
    Retrieve UniProt synonyms and description based on Ensemsbl identifiers.

    Args:
    - server          Link to UniProt REST API server.
    - ensembl_id      Ensembl, WormBase or FlyBase ID (str).
    - verbose         True/False to print logging messages.

    Returns data frame with UniProt ID, gene name, organism, sequence, sequence length, and query ID.
    """
    # API documentation: https://www.uniprot.org/help/api_queries
    # Submit server request for reviewed entries
    r = requests.get(server + ensembl_id + "+AND+reviewed:true")
    if not r.ok:
        logger.error(
            f"UniProt server request returned with error status code: {r.status_code}. Please double-check arguments or try again later."
        )
    # Convert to json
    json = r.json()

    # If no reviewed entries were found, try again for unreviewed entries
    if not len(json["results"]) > 0:
        # Submit server request
        r = requests.get(server + ensembl_id)
        if not r.ok:
            logger.error(
                f"UniProt server request returned with error status code: {r.status_code}. Please double-check arguments or try again later."
            )
        # Convert to json
        json = r.json()

        # Warn user if unreviewed results were found
        if len(json["results"]) > 0:
            if verbose is True:
                logger.warning(
                    f"No reviewed UniProt results were found for ID {ensembl_id}. Returning all unreviewed results."
                )

    if len(json["results"]) > 0:
        # Convert results to data frame
        df = pd.json_normalize(json["results"])

        # Remove non-relevant columns
        df = df[
            [
                "primaryAccession",
            ]
        ]
        # Rename column
        df.columns = [
            "uniprot_id",
        ]

        # Get primary gene name for each result
        gene_names = []
        for i in np.arange(len(json["results"])):
            try:
                gene_names.append(json["results"][i]["genes"][0]["geneName"]["value"])
            except:
                gene_names.append(np.nan)
        df["primary_gene_name"] = gene_names

        # Get synonyms for each result
        uni_synonyms = []
        for i in np.arange(len(json["results"])):
            uni_syn_temp = []
            try:
                for syn in json["results"][i]["genes"][0]["synonyms"]:
                    uni_syn_temp.append(syn["value"])
            except:
                uni_syn_temp.append(np.nan)
            uni_synonyms.append(uni_syn_temp)
        df["uni_synonyms"] = uni_synonyms

        # Get protein names for each result
        protein_names = []
        for i in np.arange(len(json["results"])):
            try:
                protein_names.append(
                    json["results"][i]["proteinDescription"]["recommendedName"][
                        "fullName"
                    ]["value"]
                )
            except:
                protein_names.append(np.nan)
        df["protein_names"] = protein_names

        # Get descriptions for each result
        descriptions = []
        for i in np.arange(len(json["results"])):
            des_temp = []
            try:
                for text in json["results"][i]["comments"]:
                    if text["commentType"] == "FUNCTION":
                        des_temp.append(text["texts"][0]["value"])
                # Keep only unique descriptions
                des_temp = np.unique(np.array(des_temp))
                # Append all descriptions to a single string object
                des_temp = " ".join(des_temp)
            except:
                des_temp.append(np.nan)

            descriptions.append(des_temp)
        df["uniprot_description"] = descriptions

        # Get subcellular localisations for each result
        subcel_locs_final = []
        for i in np.arange(len(json["results"])):
            subcel_locs = []
            try:
                for comment_idx in np.arange(len(json["results"][i]["comments"])):
                    comment_json = json["results"][i]["comments"][comment_idx]
                    if comment_json["commentType"] == "SUBCELLULAR LOCATION":
                        for location_dict in comment_json["subcellularLocations"]:
                            subcel_locs.append(location_dict["location"]["value"])
            except:
                pass
            subcel_locs_final.append(subcel_locs)

        if any(subcel_locs_final):
            df["subcellular_localisation"] = subcel_locs_final
        else:
            # No subcellular localisation data will return as nan
            nan_list = np.empty(len(subcel_locs_final))
            nan_list[:] = np.nan
            df["subcellular_localisation"] = nan_list

        # Add query colunm
        df["query"] = ensembl_id

        # Return set of all results if more than one UniProt entry was found for this Ensembl ID
        if len(df) > 1:
            final_df = pd.DataFrame()
            for column in df.columns:
                if column == "uni_synonyms" or column == "subcellular_localisation":
                    # Flatten synonym and subcellular_localisation lists
                    syn_lists = df[column].values
                    try:
                        flat_list = [item for sublist in syn_lists for item in sublist]
                        final_df[column] = [list({value: "" for value in flat_list})]

                    except:
                        final_df[column] = [syn_lists]

                else:
                    val_list = df[column].values
                    try:
                        final_df[column] = [list({value: "" for value in val_list})]
                    except:
                        final_df[column] = [val_list]

                # Try to clean up the entries (so they are not a bunch of lists of one item)
                # I will not do this with the UniProt synonyms so I can later find the set between NCBI and UniProt synonyms
                if len(final_df[column]) == 1 and column != "uni_synonyms":
                    try:
                        final_df[column] = final_df[column][0]
                    except:
                        None

            return final_df

        # If a single result was found, return df
        else:
            return df

    else:
        return None


# This function was replaced by the faster and more complete PDB API (see get_pdb_ids below)
# def get_pdb_ids(uniprot_ids):
#     """
#     Function to fetch all PDB IDs linked to a list of UniProt IDs
#     using UniProt's ID mapping (https://www.uniprot.org/help/id_mapping#example_uniparc).
#     """

#     retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
#     session = requests.Session()
#     session.mount("https://", HTTPAdapter(max_retries=retries))

#     # Sleep interval between
#     POLLING_INTERVAL = 3

#     def check_response(response):
#         """
#         Check response of HTTP request.
#         """
#         try:
#             response.raise_for_status()
#         except requests.HTTPError:
#             logger.error(
#                 f"UniProt ID mapping to fetch PDB IDs returned HTTP Error:\n{response.json()}"
#             )
#             return

#     def post_mapping_request(uniprot_ids, from_id_type, to_id_type):
#         """
#         Post ID mapping request to UniProt ID mapping API.
#         """
#         # Post ID mapping request
#         request = requests.post(
#             f"{UNIPROT_IDMAPPING_API}/run",
#             data={"from": from_id_type, "to": to_id_type, "ids": ",".join(uniprot_ids)},
#         )

#         check_response(request)

#         # Return job ID
#         return request.json()["jobId"]

#     def check_id_mapping_results_ready(job_id):
#         """
#         Poll for status of mapping request.
#         """
#         while True:
#             request = session.get(f"{UNIPROT_IDMAPPING_API}/status/{job_id}")
#             check_response(request)
#             j = request.json()
#             if "jobStatus" in j:
#                 # Sleep for POLLING_INTERVAL seconds if job status if "RUNNING"
#                 if j["jobStatus"] == "RUNNING":
#                     logger.info("Checking if PDB IDs are available...")
#                     time.sleep(POLLING_INTERVAL)
#                 else:
#                     # Raise error if job status is other than "RUNNING"
#                     raise Exception(request["jobStatus"])
#             else:
#                 # Return True if job is done running
#                 return True

#     def get_id_mapping_results_search(job_id):
#         """
#         Get mapping results from job ID.
#         """
#         # Get ID mapping results
#         request = requests.get(
#             f"{UNIPROT_IDMAPPING_API}/results/{job_id}",
#         )

#         check_response(request)

#         # Return job ID
#         return request.json()

#     # Post mapping request
#     job_id = post_mapping_request(
#         uniprot_ids=uniprot_ids, from_id_type="UniProtKB_AC-ID", to_id_type="PDB"
#     )

#     # Get results once they are ready
#     if check_id_mapping_results_ready(job_id):
#         results = get_id_mapping_results_search(job_id)

#     # Get list of PDB IDs
#     pdb_ids = []
#     if results["results"]:
#         pdb_ids = results["results"][0].values()

#     return list(pdb_ids)


def get_pdb_ids(ens_id):
    """
    Function to fetch all PDB IDs linked to an Ensembl ID.
    using the PDBe API https://wwwdev.ebi.ac.uk/pdbe/aggregated-api/mappings/ensembl_to_pdb/[ens_id]

    API documentation:
    https://www.ebi.ac.uk/pdbe/aggregated-api/#/SIFTS/get_ensembl_to_pdb_mappings_api_mappings_ensembl_to_pdb__gene_id__get
    """

    res = requests.get(ENS_TO_PDB_API + ens_id)

    if not res.ok:
        # If no PDB IDs were found, return None
        return None

    try:
        pdb_dict = res.json()[ens_id]["mappings"]
    # If no PDB IDs were found, return None
    except KeyError:
        return None

    pdb_ids = []
    for entry in pdb_dict:
        pdb_ids.append(entry["pdb_id"])

    return sorted(list(set(pdb_ids)))


def wrap_cols_func(df, cols):
    """
    Function to wrap columns cols of a
    data frame df for easier reading.
    """
    for col in cols:
        df.loc[:, col] = df[col].str.wrap(30)

    return display(HTML(df.to_html().replace("\\n", "<br>")))


def rest_query(server, query, content_type):
    """
    Function to perform a REST API query.

    Args:
    - server        Server to query.
    - Query         Query that is passed to server.
    - content_type  Content type requested from the server.

    Returns server output.
    """

    r = requests.get(server + query, headers={"Content-Type": content_type})

    if not r.ok:
        raise RuntimeError(
            f"{server} returned error status code {r.status_code}. "
            "Please double-check arguments and try again.\n"
        )

    if content_type == "application/json":
        return r.json()
    else:
        return r.text


def post_query(server, endpoint, query):
    """
    Function to perform a POST API query.

    :param server:  Server to query .
    :param endpoint: Server endpoint
    :param query:   Query that is passed to server.

    :return: server output
    """

    r = requests.post(
        server + endpoint, json=query, headers={"Content-Type": "application/json"}
    )

    if not r.ok:
        raise RuntimeError(
            f"{server} returned error status code {r.status_code}. "
            "Please double-check arguments and try again.\n"
        )

    return r.json()


def graphql_query(server, query, variables):
    """
    Function to perform a GraphQL API query.

    Args:
    - server        Server to query.
    - query         Query that is passed to server.
    - variables     Variables that are passed to server.

    Returns server output.
    """

    r = requests.post(server, json={"query": query, "variables": variables})

    if not r.ok:
        logger.debug(
            f"Server: {server}, Query: {query}, Variables: {variables}, Response: {r.text}"
        )
        raise RuntimeError(
            f"{server} returned error status code {r.status_code}. "
            "Please double-check arguments and try again.\n"
        )

    return r.json()


def find_latest_ens_rel(database=ENSEMBL_FTP_URL):
    """
    Returns the latest Ensembl release number.

    Args:
    - database    Link to Ensembl database.
    """
    # html = requests.get(database)

    # # Raise error if status code not "OK" Response
    # if html.status_code != 200:
    #     raise RuntimeError(
    #         f"The Ensembl FTP server returned error status code {html.status_code}. Please try again."
    #     )

    # soup = BeautifulSoup(html.text, "html.parser")
    # # Find all releases
    # releases = soup.body.findAll(string=re.compile("release-"))
    # # Get release numbers
    # rels = []
    # for rel in releases:
    #     rels.append(rel.split("/")[0].split("-")[-1])

    # # Find highest release number (= latest release)
    # ENS_rel = np.array(rels).astype(int).max()

    html = requests.get(database + "VERSION")
    if html.status_code != 200:
        raise RuntimeError(
            f"The Ensembl FTP server returned error status code {html.status_code}. Please try again."
        )
    ENS_rel = int(html.text)

    return ENS_rel


def search_species_options(database=ENSEMBL_FTP_URL, release=None):
    """
    Function to find all available species core databases for gget search.

    Args:
    - release   Ensembl release for which the databases are fetched.
                (Default: latest release.)

    Returns list of available core databases.
    """
    # Find latest Ensembl release
    ENS_rel = find_latest_ens_rel(database)

    # If release != None, use user-defined Ensembl release
    if release != None:
        if release > ENS_rel:
            logger.warning(
                f"Provided Ensembl release number {release} is greater than the latest release ({ENS_rel})."
            )
        ENS_rel = release

    ## Find all available databases
    # Handle structure of invertebrate database
    if "ensemblgenomes" in database:
        databases = []
        kds = ["plants", "protists", "metazoa", "fungi"]
        for kingdom in kds:
            url = database + f"release-{ENS_rel}/{kingdom}/mysql/"
            html = requests.get(url)

            # Raise error if status code not "OK" Response
            if html.status_code != 200:
                raise RuntimeError(
                    f"The Ensembl server returned error status code {html.status_code}. Please try again."
                )

            soup = BeautifulSoup(html.text, "html.parser")

            # Find all available databases
            for subsoup in soup.body.findAll("a"):
                if "core" in subsoup["href"]:
                    databases.append(subsoup["href"].split("/")[0])

    else:
        url = database + f"release-{ENS_rel}/mysql/"
        html = requests.get(url)

        # Raise error if status code not "OK" Response
        if html.status_code != 200:
            raise RuntimeError(
                f"The Ensembl server returned error status code {html.status_code}. Please try again."
            )

        soup = BeautifulSoup(html.text, "html.parser")

        # Return list of all available databases
        databases = []
        for subsoup in soup.body.findAll("a"):
            if "core" in subsoup["href"]:
                databases.append(subsoup["href"].split("/")[0])

    return databases


def find_nv_kingdom(species, release):
    kds = ["plants", "protists", "metazoa", "fungi"]
    for kingdom in kds:
        url = ENSEMBL_FTP_URL_NV + f"release-{release}/{kingdom}/fasta/"
        html = requests.get(url)

        # Raise error if status code not "OK" Response
        if html.status_code != 200:
            raise RuntimeError(
                f"The Ensembl server returned error status code {html.status_code}. Please try again."
            )

        # Parse the html and generate a clean list of the available genomes
        soup = BeautifulSoup(html.text, "html.parser")

        sps = []
        for subsoup in soup.body.findAll("a"):
            sps.append(subsoup["href"].split("/")[0])

        # Return kingdom if species was found
        if species in sps[5:]:
            return kingdom


def ref_species_options(which, database=ENSEMBL_FTP_URL, release=None):
    """
    Function to find all available species for gget ref.

    Args:
    - which     Which type of file to check for.
                Possible entries: 'dna', 'cdna', 'gtf'.
    - database  Link to Ensembl database.
    - release   Ensembl release for which available species should be fetched.

    Returns list of available species.
    """
    # Find latest Ensembl release
    ENS_rel = find_latest_ens_rel(database)

    # If release != None, use user-defined Ensembl release
    if release != None:
        # Warn user if user-defined release is higher than the latest release
        if release > ENS_rel:
            logger.warning(
                f"Provided Ensembl release number {release} is greater than the latest release ({ENS_rel})."
            )
        ENS_rel = release

    # Handle structure of non-vertebrate database
    if "ensemblgenomes" in database:
        species_list = []
        kds = ["plants", "protists", "metazoa", "fungi"]
        for kingdom in kds:
            # Find all available species for this release and FTP type
            if which == "gtf":
                url = database + f"release-{ENS_rel}/{kingdom}/gtf/"
            elif which in ("dna", "cdna"):
                url = database + f"release-{ENS_rel}/{kingdom}/fasta/"
            html = requests.get(url)

            # Raise error if status code not "OK" Response
            if html.status_code != 200:
                raise RuntimeError(
                    f"The Ensembl server returned error status code {html.status_code}. Please try again."
                )

            # Parse the html and generate a clean list of the available genomes
            soup = BeautifulSoup(html.text, "html.parser")

            sps = []
            for subsoup in soup.body.findAll("a"):
                sps.append(subsoup["href"].split("/")[0])

            species_list.append(sps[5:])

        species_list = flatten(species_list)

    else:
        # Find all available species for this release and FTP type
        if which == "gtf":
            url = database + f"release-{ENS_rel}/gtf/"
        elif which in ("dna", "cdna"):
            url = database + f"release-{ENS_rel}/fasta/"
        html = requests.get(url)

        # Raise error if status code not "OK" Response
        if html.status_code != 200:
            raise RuntimeError(
                f"The Ensembl server returned error status code {html.status_code}. Please try again."
            )

        # Parse the html and generate a clean list of the available genomes
        soup = BeautifulSoup(html.text, "html.parser")

        sps = []
        for subsoup in soup.body.findAll("a"):
            sps.append(subsoup["href"].split("/")[0])

        species_list = sps[5:]

    # Return list of all available species
    return sorted(species_list)


def parse_blast_ref_page(handle):
    """
    Extract RID and RTOE from the NCBI 'please wait' page (handle).
    RTOE = 'Estimated time fo completion.'
    RID = 'Request ID'.

    Returns RID, RTOE

    Code partly adapted from the Biopython BLAST NCBIWWW project written
    by Jeffrey Chang (Copyright 1999), Brad Chapman, and Chris Wroe distributed under the
    Biopython License Agreement and BSD 3-Clause License
    https://github.com/biopython/biopython/blob/171697883aca6894f8367f8f20f1463ce7784d0c/LICENSE.rst
    """

    # Decode handle
    string = handle.read().decode()

    # Find RID
    idx = string.find("RID =")
    if idx == -1:
        rid = None
    else:
        jdx = string.find("\n", idx)
        rid = string[idx + len("RID =") : jdx].strip()

    # Find RTOE
    rtoe_idx = string.find("RTOE =")
    if rtoe_idx == -1:
        rtoe = None
    else:
        rtoe_jdx = string.find("\n", rtoe_idx)
        rtoe = string[rtoe_idx + len("RTOE =") : rtoe_jdx].strip()

    # If neither RID, nor RTOE were found, try to extract error message from HTML page
    if not rid and not rtoe:
        # Search for 'error msInf' div class
        i = string.find('<div class="error msInf">')
        if i != -1:
            msg = string[i + len('<div class="error msInf">') :].strip()
            msg = msg.split("</div>", 1)[0].split("\n", 1)[0].strip()
            if msg:
                raise ValueError(f"Error message from NCBI: {msg}")
        # Search for 'error' class
        i = string.find('<p class="error">')
        if i != -1:
            msg = string[i + len('<p class="error">') :].strip()
            msg = msg.split("</p>", 1)[0].split("\n", 1)[0].strip()
            if msg:
                raise ValueError(f"Error message from NCBI: {msg}")
        # Generic search for error messages
        i = string.find("Message ID#")
        if i != -1:
            # Break the message at the first HTML tag
            msg = string[i:].split("<", 1)[0].split("\n", 1)[0].strip()
            raise ValueError(f"Error message from NCBI: {msg}")
        # Raise general error, if the error layout was not recognized
        raise ValueError(
            "No request ID and no estimated time to completion were found in the NCBI 'please wait' page."
        )
    # Raise error if RTOE was found but RID was not
    elif not rid:
        raise ValueError(
            f"No request ID (RID) was found in the NCBI 'please wait' page. (Although estimated time to completion = {rtoe}.)"
        )
    # Raise error if RTOE was found but RID was not
    elif not rtoe:
        raise ValueError(
            f"No estimated time to completion was found in the NCBI 'please wait' page. (Although request ID = {rid}.)"
        )

    try:
        return rid, int(rtoe)
    except ValueError:
        raise ValueError(
            f"A non-integer estimated time to completion was found in the NCBI 'please wait' page: '{rtoe}'."
        )


def tsv_to_df(tsv_file, headers=None, skiprows=None):
    """
    Convert tsv file to dataframe format.

    Args:
    - tsv_file      File to be converted

    Returns data frame.
    """
    try:
        df = pd.read_csv(tsv_file, sep="\t", names=headers, skiprows=skiprows)
        return df

    except pd.errors.EmptyDataError:
        raise RuntimeError(f"tsv to data frame reformatting failed.")


def create_tmp_fasta(sequences):
    """
    Create temporary FASTA file from str or list of sequences.

    Args:
    - sequences     List of user input amino acid sequences

    Returns: Absolute path to temoprary FASTA file.
    """
    # Generate random ID
    random_id = str(uuid.uuid4())

    if type(sequences) == str:
        sequences = [sequences]

    with open(f"tmp_{random_id}.fa", "w") as f:
        for idx, seq in enumerate(sequences):
            f.write(f">Seq{idx}\n" + seq + "\n")

    return os.path.abspath(f"tmp_{random_id}.fa")


def remove_temp_files(files_to_delete):
    """
    Delete temporary files.

    Args:
    - files_to_delete   List of paths to files to delete.
    """
    for file in files_to_delete:
        if os.path.exists(file):
            os.remove(file)


def json_list_to_df(json_list, columns) -> pd.DataFrame:
    """
    Convert list of JSON objects to data frame.

    Args:

    - json_list     List of JSON objects.
    - columns       List of (column, key) pairs in the data frame.
                    Format of each key is "key.subkey.etc", which corresponds to json_list[i]["key"]["subkey"]["etc"].

    Returns data frame with columns as specified in keys.
    """

    tmp_columns = [[] for _ in range(len(columns))]

    for json_obj in json_list:
        for i, column_key in enumerate(columns):
            _, key = column_key
            keys_split = key.split(".")
            value = json_obj
            for k in keys_split:
                if value is None:
                    break
                value = value[k]
            tmp_columns[i].append(value)

    return pd.DataFrame({key[0]: value for key, value in zip(columns, tmp_columns)})
