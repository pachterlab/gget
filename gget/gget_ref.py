from bs4 import BeautifulSoup
import requests
import re
import numpy as np
import json
import logging

# Add and format time stamp in logging messages
logging.basicConfig(
    format="%(asctime)s %(levelname)s %(message)s",
    level=logging.INFO,
    datefmt="%c",
)
# Custom functions
from .utils import ref_species_options, find_latest_ens_rel

from .constants import ENSEMBL_FTP_URL


def ref(species, which="all", release=None, ftp=False, save=False, list_species=False):
    """
    Fetch FTPs for reference genomes and annotations by species.

    Args:
    - species       Defines the species for which the reference should be fetched in the format "<genus>_<species>",
                    e.g. species = "homo_sapiens".
    - which         Defines which results to return.
                    Default: 'all' -> Returns all available results.
                    Possible entries are one or a combination (as a list of strings) of the following:
                    'gtf' - Returns the annotation (GTF).
                    'cdna' - Returns the trancriptome (cDNA).
                    'dna' - Returns the genome (DNA).
                    'cds - Returns the coding sequences corresponding to Ensembl genes. (Does not contain UTR or intronic sequence.)
                    'cdrna' - Returns transcript sequences corresponding to non-coding RNA genes (ncRNA).
                    'pep' - Returns the protein translations of Ensembl genes.
    - release       Defines the Ensembl release number from which the files are fetched, e.g. release = 104.
                    (Ensembl releases earlier than release 48 are not supported.)
                    By default, the latest Ensembl release is used.
    - ftp           Return only the requested FTP links in a list (default: False).
    - save          Save the results in the local directory (default: False).

    - list_species  If True and `species=None`, returns a list of all available species (default: False).
                    (Can be combined with `release` to get the available species from a specific Ensembl release.)

    Returns a dictionary containing the requested URLs with their respective Ensembl version and release date and time.
    (If FTP=True, returns a list containing only the URLs.)
    """

    # Return list of all available species
    if list_species is True:
        # Find all available species for GTFs for this Ensembl release
        species_list_gtf = ref_species_options("gtf", release=release)
        # Find all available species for FASTAs for this Ensembl release
        species_list_dna = ref_species_options("dna", release=release)

        # Find intersection of the two lists
        # (Only species which have GTF and FASTAs available can continue)
        species_list = list(set(species_list_gtf) & set(species_list_dna))

        if release is None:
            logging.info(
                f"Fetching available genomes (GTF and FASTAs present) from Ensembl release {find_latest_ens_rel()} (latest)."
            )

        else:
            logging.info(
                f"Fetching available genomes (GTF and FASTAs present) from Ensembl release {release}."
            )

        return sorted(species_list)

    # Species shortcuts
    if species == "human":
        species = "homo_sapiens"
    if species == "mouse":
        species = "mus_musculus"

    # In case species was passed with upper case letters
    species = species.lower()

    ## Find latest Ensembl release
    url = ENSEMBL_FTP_URL
    html = requests.get(url)

    # Raise error if status code not "OK" Response
    if html.status_code != 200:
        raise RuntimeError(
            f"Ensembl FTP site returned error status code {html.status_code}. Please try again."
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

    # If release != None, use user-defined Ensembl release
    if release != None:
        # Do not allow user-defined release if it is higher than the latest release
        if release > ENS_rel:
            raise ValueError(
                f"Defined Ensembl release number cannot be greater than latest release ({ENS_rel}).\n"
            )
        else:
            ENS_rel = release

    ## Raise error if species not found
    # Find all available species for GTFs for this Ensembl release
    species_list_gtf = ref_species_options("gtf", release=ENS_rel)
    # Find all available species for FASTAs for this Ensembl release
    species_list_dna = ref_species_options("dna", release=ENS_rel)

    # Find intersection of the two lists
    # (Only species which have GTF and FASTAs available can continue)
    species_list = list(set(species_list_gtf) & set(species_list_dna))

    if species not in species_list:
        raise ValueError(
            f"Species does not match any available species for Ensembl release {ENS_rel}.\n"
            "Please double-check spelling.\n"
            "'gget ref --list' -> lists out all available species (Python: 'gget.ref(None, list_species=True)').\n"
            "Combine with `release` argument to define specific Ensembl release (default: latest).\n"
        )

    ## Get GTF link for this species and release
    url = ENSEMBL_FTP_URL + f"release-{ENS_rel}/gtf/{species}/"
    html = requests.get(url)

    # Raise error if status code not "OK" Response
    if html.status_code != 200:
        raise RuntimeError(
            f"HTTP response status code {html.status_code}. Please try again.\n"
        )

    soup = BeautifulSoup(html.text, "html.parser")

    # The url can be found under an <a> object tag in the html,
    # but the date and size do not have an object tag (element=None)
    nones = []
    a_elements = []
    pre = soup.find("pre")
    for element in pre.descendants:
        if element.name == "a":
            a_elements.append(element)
        elif element.name != "a":
            nones.append(element)

    # Find the <a> element containing the url
    for i, string in enumerate(a_elements):
        if f"{ENS_rel}.gtf.gz" in string["href"]:
            gtf_str = string
            # Get release date and time from <None> elements (since there are twice as many, 2x and +1 to move from string to date)
            gtf_date_size = nones[i * 2 + 1]

    gtf_url = ENSEMBL_FTP_URL + f"release-{ENS_rel}/gtf/{species}/{gtf_str['href']}"

    gtf_date = gtf_date_size.strip().split("  ")[0]
    gtf_size = gtf_date_size.strip().split("  ")[-1]

    ## Get cDNA FASTA link for this species and release
    url = ENSEMBL_FTP_URL + f"release-{ENS_rel}/fasta/{species}/cdna"
    html = requests.get(url)

    # Raise error if status code not "OK" Response
    if html.status_code != 200:
        raise RuntimeError(
            f"HTTP response status code {html.status_code}. Please try again.\n"
        )

    soup = BeautifulSoup(html.text, "html.parser")

    # The url can be found under an <a> object tag in the html,
    # but the date and size do not have an object tag (element=None)
    nones = []
    a_elements = []
    pre = soup.find("pre")
    for element in pre.descendants:
        if element.name == "a":
            a_elements.append(element)
        elif element.name != "a":
            nones.append(element)

    # Find the <a> element containing the url
    for i, string in enumerate(a_elements):
        if "cdna.all.fa" in string["href"]:
            cdna_str = string
            # Get release date and time from <None> elements (since there are twice as many, 2x and +1 to move from string to date)
            cdna_date_size = nones[i * 2 + 1]

    cdna_url = (
        ENSEMBL_FTP_URL + f"release-{ENS_rel}/fasta/{species}/cdna/{cdna_str['href']}"
    )

    cdna_date = cdna_date_size.strip().split("  ")[0]
    cdna_size = cdna_date_size.strip().split("  ")[-1]

    ## Get DNA FASTA link for this species and release
    url = ENSEMBL_FTP_URL + f"release-{ENS_rel}/fasta/{species}/dna"
    html = requests.get(url)

    # Raise error if status code not "OK" Response
    if html.status_code != 200:
        raise RuntimeError(
            f"HTTP response status code {html.status_code}. Please try again.\n"
        )

    soup = BeautifulSoup(html.text, "html.parser")

    # The url can be found under an <a> object tag in the html,
    # but the date and size do not have an object tag (element=None)
    nones = []
    a_elements = []
    pre = soup.find("pre")

    for element in pre.descendants:
        if element.name == "a":
            a_elements.append(element)
        elif element.name != "a":
            nones.append(element)

    # Get primary assembly if available, otherwise toplevel assembly
    dna_str = None
    for i, string in enumerate(a_elements):
        if ".dna.primary_assembly.fa" in string["href"]:
            dna_str = string
            # Get date from non-assigned values (since there are twice as many, 2x and +1 to move from string to date)
            dna_date_size = nones[i * 2 + 1]

    # Find the <a> element containing the url
    if dna_str is None:
        for i, string in enumerate(a_elements):
            if ".dna.toplevel.fa" in string["href"]:
                dna_str = string
                # Get date from non-assigned values (since there are twice as many, 2x and +1 to move from string to date)
                dna_date_size = nones[i * 2 + 1]

    dna_url = (
        ENSEMBL_FTP_URL + f"release-{ENS_rel}/fasta/{species}/dna/{dna_str['href']}"
    )

    dna_date = dna_date_size.strip().split("  ")[0]
    dna_size = dna_date_size.strip().split("  ")[-1]
    # Strip again to remove any extra spaces
    dna_date = dna_date.strip()
    dna_size = dna_size.strip()

    ## Get CDS FASTA link for this species and release
    url = ENSEMBL_FTP_URL + f"release-{ENS_rel}/fasta/{species}/cds"
    html = requests.get(url)

    # Raise error if status code not "OK" Response
    if html.status_code != 200:
        raise RuntimeError(
            f"HTTP response status code {html.status_code}. Please try again.\n"
        )

    soup = BeautifulSoup(html.text, "html.parser")

    # The url can be found under an <a> object tag in the html,
    # but the date and size do not have an object tag (element=None)
    nones = []
    a_elements = []
    pre = soup.find("pre")
    for element in pre.descendants:
        if element.name == "a":
            a_elements.append(element)
        elif element.name != "a":
            nones.append(element)

    # Find the <a> element containing the url
    for i, string in enumerate(a_elements):
        if "cds.all.fa" in string["href"]:
            cds_str = string
            # Get release date and time from <None> elements (since there are twice as many, 2x and +1 to move from string to date)
            cds_date_size = nones[i * 2 + 1]

    cds_url = (
        ENSEMBL_FTP_URL + f"release-{ENS_rel}/fasta/{species}/cds/{cds_str['href']}"
    )

    cds_date = cds_date_size.strip().split("  ")[0]
    cds_size = cds_date_size.strip().split("  ")[-1]

    ## Get ncRNA FASTA link for this species and release
    url = ENSEMBL_FTP_URL + f"release-{ENS_rel}/fasta/{species}/ncrna"
    html = requests.get(url)

    # Raise error if status code not "OK" Response
    if html.status_code != 200:
        raise RuntimeError(
            f"HTTP response status code {html.status_code}. Please try again.\n"
        )

    soup = BeautifulSoup(html.text, "html.parser")

    # The url can be found under an <a> object tag in the html,
    # but the date and size do not have an object tag (element=None)
    nones = []
    a_elements = []
    pre = soup.find("pre")
    for element in pre.descendants:
        if element.name == "a":
            a_elements.append(element)
        elif element.name != "a":
            nones.append(element)

    # Find the <a> element containing the url
    for i, string in enumerate(a_elements):
        if ".ncrna.fa" in string["href"]:
            ncrna_str = string
            # Get release date and time from <None> elements (since there are twice as many, 2x and +1 to move from string to date)
            ncrna_date_size = nones[i * 2 + 1]

    ncrna_url = (
        ENSEMBL_FTP_URL + f"release-{ENS_rel}/fasta/{species}/ncrna/{ncrna_str['href']}"
    )

    ncrna_date = ncrna_date_size.strip().split("  ")[0]
    ncrna_size = ncrna_date_size.strip().split("  ")[-1]

    ## Get pep FASTA link for this species and release
    url = ENSEMBL_FTP_URL + f"release-{ENS_rel}/fasta/{species}/pep"
    html = requests.get(url)

    # Raise error if status code not "OK" Response
    if html.status_code != 200:
        raise RuntimeError(
            f"HTTP response status code {html.status_code}. Please try again.\n"
        )

    soup = BeautifulSoup(html.text, "html.parser")

    # The url can be found under an <a> object tag in the html,
    # but the date and size do not have an object tag (element=None)
    nones = []
    a_elements = []
    pre = soup.find("pre")
    for element in pre.descendants:
        if element.name == "a":
            a_elements.append(element)
        elif element.name != "a":
            nones.append(element)

    # Find the <a> element containing the url
    for i, string in enumerate(a_elements):
        if ".pep.all.fa" in string["href"]:
            pep_str = string
            # Get release date and time from <None> elements (since there are twice as many, 2x and +1 to move from string to date)
            pep_date_size = nones[i * 2 + 1]

    pep_url = (
        ENSEMBL_FTP_URL + f"release-{ENS_rel}/fasta/{species}/pep/{pep_str['href']}"
    )

    pep_date = pep_date_size.strip().split("  ")[0]
    pep_size = pep_date_size.strip().split("  ")[-1]

    ## Return results
    # If single which passed as string, convert to list
    if type(which) == str:
        which = [which]

    # Raise error if several values are passed and 'all' is included
    if len(which) > 1 and "all" in which:
        raise ValueError(
            "Parameter 'which' must be 'all', or any one or a combination of the following: 'gtf', 'cdna', 'dna', 'cds', 'ncrna', 'pep'.\n"
        )

    # If FTP=False, return dictionary/json of specified results
    if ftp is False:
        ref_dict = {species: {}}
        for return_val in which:
            if return_val == "all":
                ref_dict = {
                    species: {
                        "transcriptome_cdna": {
                            "ftp": cdna_url,
                            "ensembl_release": int(ENS_rel),
                            "release_date": cdna_date.split(" ")[0],
                            "release_time": cdna_date.split(" ")[1],
                            "bytes": cdna_size,
                        },
                        "genome_dna": {
                            "ftp": dna_url,
                            "ensembl_release": int(ENS_rel),
                            "release_date": dna_date.split(" ")[0],
                            "release_time": dna_date.split(" ")[1],
                            "bytes": dna_size,
                        },
                        "annotation_gtf": {
                            "ftp": gtf_url,
                            "ensembl_release": int(ENS_rel),
                            "release_date": gtf_date.split(" ")[0],
                            "release_time": gtf_date.split(" ")[1],
                            "bytes": gtf_size,
                        },
                        "coding_seq_cds": {
                            "ftp": cds_url,
                            "ensembl_release": int(ENS_rel),
                            "release_date": cds_date.split(" ")[0],
                            "release_time": cds_date.split(" ")[1],
                            "bytes": cds_size,
                        },
                        "non-coding_seq_ncRNA": {
                            "ftp": ncrna_url,
                            "ensembl_release": int(ENS_rel),
                            "release_date": ncrna_date.split(" ")[0],
                            "release_time": ncrna_date.split(" ")[1],
                            "bytes": ncrna_size,
                        },
                        "protein_translation_pep": {
                            "ftp": pep_url,
                            "ensembl_release": int(ENS_rel),
                            "release_date": pep_date.split(" ")[0],
                            "release_time": pep_date.split(" ")[1],
                            "bytes": pep_size,
                        },
                    }
                }
            elif return_val == "gtf":
                dict_temp = {
                    "annotation_gtf": {
                        "ftp": gtf_url,
                        "ensembl_release": int(ENS_rel),
                        "release_date": gtf_date.split(" ")[0],
                        "release_time": gtf_date.split(" ")[1],
                        "bytes": gtf_size,
                    },
                }
                ref_dict[species].update(dict_temp)
            elif return_val == "cdna":
                dict_temp = {
                    "transcriptome_cdna": {
                        "ftp": cdna_url,
                        "ensembl_release": int(ENS_rel),
                        "release_date": cdna_date.split(" ")[0],
                        "release_time": cdna_date.split(" ")[1],
                        "bytes": cdna_size,
                    },
                }
                ref_dict[species].update(dict_temp)
            elif return_val == "dna":
                dict_temp = {
                    "genome_dna": {
                        "ftp": dna_url,
                        "ensembl_release": int(ENS_rel),
                        "release_date": dna_date.split(" ")[0],
                        "release_time": dna_date.split(" ")[1],
                        "bytes": dna_size,
                    },
                }
                ref_dict[species].update(dict_temp)
            elif return_val == "cds":
                dict_temp = {
                    "coding_seq_cds": {
                        "ftp": cds_url,
                        "ensembl_release": int(ENS_rel),
                        "release_date": cds_date.split(" ")[0],
                        "release_time": cds_date.split(" ")[1],
                        "bytes": cds_size,
                    },
                }
                ref_dict[species].update(dict_temp)
            elif return_val == "ncrna":
                dict_temp = {
                    "non-coding_seq_ncRNA": {
                        "ftp": ncrna_url,
                        "ensembl_release": int(ENS_rel),
                        "release_date": ncrna_date.split(" ")[0],
                        "release_time": ncrna_date.split(" ")[1],
                        "bytes": ncrna_size,
                    },
                }
                ref_dict[species].update(dict_temp)
            elif return_val == "pep":
                dict_temp = {
                    "protein_translation_pep": {
                        "ftp": pep_url,
                        "ensembl_release": int(ENS_rel),
                        "release_date": pep_date.split(" ")[0],
                        "release_time": pep_date.split(" ")[1],
                        "bytes": pep_size,
                    },
                }
                ref_dict[species].update(dict_temp)
            else:
                raise ValueError(
                    "Parameter 'which' must be 'all', or any one or a combination of the following: 'gtf', 'cdna', 'dna', 'cds', 'ncrna', 'pep'.\n"
                )

        if save:
            with open("ref_results.json", "w", encoding="utf-8") as file:
                json.dump(ref_dict, file, ensure_ascii=False, indent=4)

        logging.info(
            f"Fetching reference information for {species} from Ensembl release: {ENS_rel}."
        )
        return ref_dict

    # If FTP==True, return only the specified URLs as a list
    if ftp:
        logging.info(
            f"Fetching reference information for {species} from Ensembl release: {ENS_rel}."
        )
        results = []
        for return_val in which:
            if return_val == "all":
                results.append(gtf_url)
                results.append(cdna_url)
                results.append(dna_url)
                results.append(cds_url)
                results.append(ncrna_url)
                results.append(pep_url)
            elif return_val == "gtf":
                results.append(gtf_url)
            elif return_val == "cdna":
                results.append(cdna_url)
            elif return_val == "dna":
                results.append(dna_url)
            elif return_val == "cds":
                results.append(cds_url)
            elif return_val == "ncrna":
                results.append(ncrna_url)
            elif return_val == "pep":
                results.append(pep_url)
            else:
                raise ValueError(
                    "Parameter 'which' must be 'all', or any one or a combination of the following: 'gtf', 'cdna', 'dna', 'cds', 'ncrna', 'pep'.\n"
                )

        if save:
            with open("ref_results.json", "w", encoding="utf-8") as f:
                json.dump(results_dict, f, ensure_ascii=False, indent=4)

        return results
