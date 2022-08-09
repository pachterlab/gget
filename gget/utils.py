from bs4 import BeautifulSoup
import requests
import re
import pandas as pd
import numpy as np
import urllib
from io import StringIO
from IPython.display import display, HTML
import logging

# Add and format time stamp in logging messages
logging.basicConfig(
    format="%(asctime)s %(levelname)s %(message)s",
    level=logging.INFO,
    datefmt="%c",
)
# Mute numexpr threads info
logging.getLogger("numexpr").setLevel(logging.WARNING)

from .constants import ENSEMBL_FTP_URL


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
            logging.error(
                f"UniProt server request returned with error status code: {r.status_code}. Please double-check arguments or try again later."
            )
        # Convert to json
        json = r.json()

        # If no reviewed results were found, try again for unreviewed results
        if not len(json["results"]) > 0:
            # Submit server request
            r = requests.get(server + id_)
            if not r.ok:
                logging.error(
                    f"UniProt server request returned with error status code: {r.status_code}. Please double-check arguments or try again later."
                )
            # Convert to json
            json = r.json()

            # Warn user if unreviewed results were found
            if len(json["results"]) > 0:
                logging.warning(
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
                    gene_names.append(np.NaN)
            df["gene_name"] = gene_names
            df["query"] = id_

            # Append results for this ID to master data frame
            master_df = pd.concat([master_df, df], axis=0)

        else:
            # If no results were found, warn user and do nothing -> returns empty df
            logging.warning(f"No UniProt sequences were found for ID {id_}.")

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
        logging.error(
            f"UniProt server request returned with error status code: {r.status_code}. Please double-check arguments or try again later."
        )
    # Convert to json
    json = r.json()

    # If no reviewed entries were found, try again for unreviewed entries
    if not len(json["results"]) > 0:
        # Submit server request
        r = requests.get(server + ensembl_id)
        if not r.ok:
            logging.error(
                f"UniProt server request returned with error status code: {r.status_code}. Please double-check arguments or try again later."
            )
        # Convert to json
        json = r.json()

        # Warn user if unreviewed results were found
        if len(json["results"]) > 0:
            if verbose is True:
                logging.warning(
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
                gene_names.append(np.NaN)
        df["primary_gene_name"] = gene_names

        # Get synonyms for each result
        uni_synonyms = []
        for i in np.arange(len(json["results"])):
            uni_syn_temp = []
            try:
                for syn in json["results"][i]["genes"][0]["synonyms"]:
                    uni_syn_temp.append(syn["value"])
            except:
                uni_syn_temp.append(np.NaN)
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
                protein_names.append(np.NaN)
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
                des_temp.append(np.NaN)

            descriptions.append(des_temp)
        df["uniprot_description"] = descriptions

        # Add query colunm
        df["query"] = ensembl_id

        # Return set of all results if more than one UniProt entry was found for this Ensembl ID
        if len(df) > 1:
            final_df = pd.DataFrame()
            for column in df.columns:
                if column == "uni_synonyms":
                    # Flatten synonym lists
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


def find_latest_ens_rel(database=ENSEMBL_FTP_URL):
    """
    Returns the latest Ensembl release number.

    Args:
    - database    Link to Ensembl database.
    """
    html = requests.get(database)

    # Raise error if status code not "OK" Response
    if html.status_code != 200:
        raise RuntimeError(
            f"The Ensembl FTP server returned error status code {html.status_code}. Please try again."
        )

    soup = BeautifulSoup(html.text, "html.parser")
    # Find all releases
    releases = soup.body.findAll(text=re.compile("release-"))
    # Get release numbers
    rels = []
    for rel in releases:
        rels.append(rel.split("/")[0].split("-")[-1])

    # Find highest release number (= latest release)
    ENS_rel = np.array(rels).astype(int).max()

    return ENS_rel


def gget_species_options(release=None):
    """
    Function to find all available species core databases for gget.

    Args:
    - release   Ensembl release for which the databases are fetched.
                (Default: latest release.)

    Returns list of available core databases.
    """
    # Find latest Ensembl release
    ENS_rel = find_latest_ens_rel()

    # If release != None, use user-defined Ensembl release
    if release != None:
        # Do not allow user-defined release if it is higher than the latest release
        if release > ENS_rel:
            raise ValueError(
                "Defined Ensembl release number cannot be greater than latest release."
            )
        else:
            ENS_rel = release

    # Find all available databases
    url = ENSEMBL_FTP_URL + f"release-{ENS_rel}/mysql/"
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


def ref_species_options(which, database=ENSEMBL_FTP_URL, release=None):
    """
    Function to find all available species for gget ref.

    Args:
    - which     Which type of FTP. Possible entries: 'dna', 'cdna', 'gtf'.
    - database  Link to Ensembl database.
    - release   Ensembl release for which available species should be fetched.

    Returns list of available species.
    """
    # Find latest Ensembl release
    ENS_rel = find_latest_ens_rel(database)

    # If release != None, use user-defined Ensembl release
    if release != None:
        # Do not allow user-defined release if it is higher than the latest release
        if release > ENS_rel:
            raise ValueError(
                "Defined Ensembl release number cannot be greater than latest release."
            )
        else:
            ENS_rel = release

    # Find all available species for this release and FTP type
    if which == "gtf":
        url = database + f"release-{ENS_rel}/gtf/"
    if which == "dna" or which == "cdna":
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
    return species_list


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
