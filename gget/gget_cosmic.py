import requests
import pandas as pd
import subprocess
import os
import re
import json as json_package
import base64
import shutil
import tarfile
import gzip
import getpass

# Constants
from .constants import COSMIC_GET_URL
from .utils import set_up_logger, get_latest_cosmic

logger = set_up_logger()


def is_valid_email(email):
    """
    Check if an e-mail address is valid.
    """
    email_pattern = re.compile(r"(^[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+$)")

    return re.match(email_pattern, email) is not None


def download_reference(download_link, tar_folder_path, file_path, verbose):
    email = input("Please enter your COSMIC email: ")
    if not is_valid_email(email):
        raise ValueError("The email address is not valid.")
    password = getpass.getpass("Please enter your COSMIC password: ")

    # Concatenate the email and password with a colon
    input_string = f"{email}:{password}\n"

    encoded_bytes = base64.b64encode(input_string.encode("utf-8"))
    encoded_string = encoded_bytes.decode("utf-8")
    curl_command = [
        "curl",
        "-H",
        f"Authorization: Basic {encoded_string}",
        download_link,
    ]

    if verbose:
        logger.info("Downloading data...")

    result = subprocess.run(curl_command, capture_output=True, text=True)

    try:
        response_data = json_package.loads(result.stdout)
    except ValueError:
        raise RuntimeError(
            "Failed to download file. Please double-check arguments (especially cosmic_version) and try again."
        )
    try:
        true_download_url = response_data.get("url")
    except AttributeError:
        raise AttributeError("Invalid username or password.")

    curl_command2 = ["curl", true_download_url, "--output", f"{tar_folder_path}.tar"]
    result2 = subprocess.run(curl_command2, capture_output=True, text=True)

    if result2.returncode != 0:
        raise RuntimeError(
            f"Failed to download file. Return code: {result.returncode}\n{result.stderr}"
        )

    with tarfile.open(f"{tar_folder_path}.tar", "r") as tar:
        tar.extractall(path=tar_folder_path)
        if verbose:
            logger.info(f"Extracted tar file to {tar_folder_path}")

    with gzip.open(f"{file_path}.gz", "rb") as f_in:
        with open(file_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        if verbose:
            logger.info(f"Unzipped file to {file_path}")


def select_reference(
    mutation_class, reference_dir, grch_version, cosmic_version, verbose
):
    # if mutation_class == "transcriptome":
    #     download_link = f"https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch{grch_version}/cosmic/v{cosmic_version}/Cosmic_Genes_Fasta_v{cosmic_version}_GRCh{grch_version}.tar&bucket=downloads"
    #     tarred_folder = f"Cosmic_Genes_Fasta_v{cosmic_version}_GRCh{grch_version}"
    #     contained_file = f"Cosmic_Genes_v{cosmic_version}_GRCh{grch_version}.fasta"

    if mutation_class == "cancer":
        if grch_version == 38:
            logger.error(
                "CancerMutationCensus data is only available for GRCh37. Define grch_version=37."
            )
        download_link = f"https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=GRCh{grch_version}/cmc/v{cosmic_version}/CancerMutationCensus_AllData_Tsv_v{cosmic_version}_GRCh{grch_version}.tar&bucket=downloads"
        tarred_folder = (
            f"CancerMutationCensus_AllData_Tsv_v{cosmic_version}_GRCh{grch_version}"
        )
        contained_file = (
            f"CancerMutationCensus_AllData_v{cosmic_version}_GRCh{grch_version}.tsv"
        )

    elif mutation_class == "cell_line":
        download_link = f"https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch{grch_version}/cell_lines/v{cosmic_version}/CellLinesProject_GenomeScreensMutant_Tsv_v{cosmic_version}_GRCh{grch_version}.tar&bucket=downloads"
        tarred_folder = f"CellLinesProject_GenomeScreensMutant_Tsv_v{cosmic_version}_GRCh{grch_version}"
        contained_file = f"CellLinesProject_GenomeScreensMutant_v{cosmic_version}_GRCh{grch_version}.tsv"

    elif mutation_class == "census":
        download_link = f"https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch{grch_version}/cosmic/v{cosmic_version}/Cosmic_MutantCensus_Tsv_v{cosmic_version}_GRCh{grch_version}.tar&bucket=downloads"
        tarred_folder = f"Cosmic_MutantCensus_Tsv_v{cosmic_version}_GRCh{grch_version}"
        contained_file = f"Cosmic_MutantCensus_v{cosmic_version}_GRCh{grch_version}.tsv"

    elif mutation_class == "resistance":
        download_link = f"https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch{grch_version}/cosmic/v{cosmic_version}/Cosmic_ResistanceMutations_Tsv_v{cosmic_version}_GRCh{grch_version}.tar&bucket=downloads"
        tarred_folder = (
            f"Cosmic_ResistanceMutations_Tsv_v{cosmic_version}_GRCh{grch_version}"
        )
        contained_file = (
            f"Cosmic_ResistanceMutations_v{cosmic_version}_GRCh{grch_version}.tsv"
        )

    elif mutation_class == "genome_screen":
        download_link = f"https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch{grch_version}/cosmic/v{cosmic_version}/Cosmic_GenomeScreensMutant_Tsv_v{cosmic_version}_GRCh{grch_version}.tar&bucket=downloads"
        tarred_folder = (
            f"Cosmic_GenomeScreensMutant_Tsv_v{cosmic_version}_GRCh{grch_version}"
        )
        contained_file = (
            f"Cosmic_GenomeScreensMutant_v{cosmic_version}_GRCh{grch_version}.tsv"
        )

    elif mutation_class == "targeted_screen":
        download_link = f"https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch{grch_version}/cosmic/v{cosmic_version}/Cosmic_CompleteTargetedScreensMutant_Tsv_v{cosmic_version}_GRCh{grch_version}.tar&bucket=downloads"
        tarred_folder = f"Cosmic_CompleteTargetedScreensMutant_Tsv_v{cosmic_version}_GRCh{grch_version}"
        contained_file = f"Cosmic_CompleteTargetedScreensMutant_v{cosmic_version}_GRCh{grch_version}.tsv"

    # Only available for the latest COSMIC version
    elif mutation_class == "cancer_example":
        download_link = f"https://cog.sanger.ac.uk/cosmic-downloads-production/taster/example_grch{grch_version}.tar"
        tarred_folder = f"example_GRCh{grch_version}"
        contained_file = f"CancerMutationCensus_AllData_v{get_latest_cosmic()}_GRCh{grch_version}.tsv"

    tar_folder_path = os.path.join(reference_dir, tarred_folder)
    file_path = os.path.join(tar_folder_path, contained_file)

    overwrite = True
    if os.path.exists(file_path):
        proceed = (
            input(
                "The requested COSMIC database already exists at the destination. Would you like to overwrite the existing files (y/n)? "
            )
            .strip()
            .lower()
        )
        if proceed in ["yes", "y"]:
            overwrite = True
        else:
            overwrite = False

    if overwrite:
        # Only the example database can be downloaded directly (without an account)
        if mutation_class == "cancer_example":
            curl_command = [
                "curl",
                "-L",
                "--output",
                f"{tar_folder_path}.tar",
                download_link,
            ]
            result = subprocess.run(curl_command, capture_output=True, text=True)

            with tarfile.open(f"{tar_folder_path}.tar", "r") as tar:
                tar.extractall(path=tar_folder_path)
                if verbose:
                    logger.info(f"Extracted tar file to {tar_folder_path}")

        # Download full databases
        else:
            proceed = (
                input(
                    "Downloading complete databases from COSMIC requires an account (https://cancer.sanger.ac.uk/cosmic/register; free for academic use, license for commercial use).\nWould you like to proceed (y/n)? "
                )
                .strip()
                .lower()
            )
            if proceed in ["yes", "y"]:
                download_reference(download_link, tar_folder_path, file_path, verbose)
            else:
                raise KeyboardInterrupt(
                    f"Database download canceled. Learn more about COSMIC at https://cancer.sanger.ac.uk/cosmic/download/cosmic."
                )

    return file_path, overwrite


def convert_chromosome_value_to_int_when_possible(val):
    try:
        # Try to convert the value to a float, then to an int, and finally to a string
        return str(int(float(val)))
    except ValueError:
        # If conversion fails, keep the value as it is
        return val


def cosmic(
    searchterm,
    entity="mutations",
    limit=100,
    json=False,
    download_cosmic=False,
    mutation_class="cancer",
    cosmic_version=None,
    grch_version=37,
    gget_mutate=True,
    keep_genome_info=False,
    remove_duplicates=False,
    out=None,
    verbose=True,
):
    """
    Search for genes, mutations, etc associated with cancers using the COSMIC
    (Catalogue Of Somatic Mutations In Cancer) database
    (https://cancer.sanger.ac.uk/cosmic).
    NOTE: Licence fees apply for the commercial use of COSMIC.

    Args for querying information about specific cancers/genes/etc:
    - searchterm      (str) Search term, which can be a mutation, gene name (or Ensembl ID), sample, etc.
                      Examples for the searchterm and entitity arguments:

                      | searchterm   | entitity    |
                      |--------------|-------------|
                      | EGFR         | mutations   | -> Find mutations in the EGFR gene that are associated with cancer
                      | v600e        | mutations   | -> Find genes for which a v600e mutation is associated with cancer
                      | COSV57014428 | mutations   | -> Find mutations associated with this COSMIC mutations ID
                      | EGFR         | genes       | -> Get the number of samples, coding/simple mutations, and fusions observed in COSMIC for EGFR
                      | prostate     | cancer      | -> Get number of tested samples and mutations for prostate cancer
                      | prostate     | tumour_site | -> Get number of tested samples, genes, mutations, fusions, etc. with 'prostate' as primary tissue site
                      | ICGC         | studies     | -> Get project code and descriptions for all studies from the ICGC (International Cancer Genome Consortium)
                      | EGFR         | pubmed      | -> Find PubMed publications on EGFR and cancer
                      | ICGC         | samples     | -> Get metadata on all samples from the ICGC (International Cancer Genome Consortium)
                      | COSS2907494  | samples     | -> Get metadata on this COSMIC sample ID (cancer type, tissue, # analyzed genes, # mutations, etc.)

                      NOTE: Set to None when downloading COSMIC databases with download_cosmic=True.
    - entity          (str) Defines the type of the results to return. One of the following:
                      'mutations' (default), 'genes', 'cancer', 'tumour_site', 'studies', 'pubmed', or 'samples'.
    - limit           (int) Number of hits to return. Default: 100
    - json            (True/False) If True, returns results in json format instead of data frame. Default: False

    Returns a data frame with the requested results.

    Args for downloading COSMIC databases:
    NOTE: Downloading complete databases from COSMIC requires an account (https://cancer.sanger.ac.uk/cosmic/register; free for academic use, license for commercial use)
    - download_cosmic   (True/False) whether to switch into database download mode. Default: False
    - mutation_class    (str) Type of COSMIC database to download. One of the following:
                        'cancer' (default), 'cell_line', 'census', 'resistance', 'genome_screen', 'targeted_screen', 'cancer_example'
    - cosmic_version    (int) Version of the COSMIC database. Default: None -> Defaults to latest version.
    - grch_version      (int) Version of the human GRCh reference genome the COSMIC database was based on (37 or 38). Default: 37
    - gget_mutate       (True/False) Whether to create a modified version of the database for use with gget mutate. Default: True
    - keep_genome_info  (True/False) Whether to keep genome information (e.g. location of mutation in the genome) in the modified database for use with gget mutate. Default: False
    - remove_duplicates (True/False) Whether to remove duplicate rows from the modified database for use with gget mutate. Default: False

    General args:
    - out             (str) Path to the file (or folder when downloading databases with the download_cosmic flag) the results will be saved in, e.g. 'path/to/results.json'.
                      Default: None
                      -> When download_cosmic=False: Results will be returned to standard out
                      -> When download_cosmic=True: Database will be downloaded into current working directory
    - verbose         (True/False) whether to print progress information. Default: True

    Saves the requested database into the specified folder (or current working directory if out=None).
    """

    if verbose:
        logger.info("NOTE: Licence fees apply for the commercial use of COSMIC.")

    ## Database download
    if download_cosmic:
        mut_class_allowed = [
            "cancer",
            "cell_line",
            "census",
            "resistance",
            "genome_screen" "targeted_screen",
            "cancer_example",
        ]
        if mutation_class not in mut_class_allowed:
            raise ValueError(
                f"Parameter 'mutation_class' must be one of the following: {', '.join(mut_class_allowed)}.\n"
            )

        grch_allowed = [37, 38]
        if grch_version not in grch_allowed:
            raise ValueError(
                f"Parameter 'grch_version' must be one of the following: {', '.join(grch_allowed)}.\n"
            )

        if not out:
            out = os.getcwd()

        if not os.path.exists(out):
            os.makedirs(out)

        if not cosmic_version:
            cosmic_version = get_latest_cosmic()
            if verbose:
                logger.info(
                    f"Downloading data from latest COSMIC version (v{cosmic_version})."
                )

        ## Download requested database
        mutation_tsv_file, overwrite = select_reference(
            mutation_class, out, grch_version, cosmic_version, verbose
        )

        if gget_mutate and overwrite:
            ## Create copy of results formatted for further use by gget mutate
            if verbose:
                logger.info(
                    "Creating modified mutations file for use with gget mutate..."
                )

            if mutation_class == "cancer" or mutation_class == "cancer_example":
                relevant_cols = [
                    "GENE_NAME",
                    "ACCESSION_NUMBER",
                    "MUTATION_URL",
                    "Mutation CDS",
                    "Mutation AA",
                ]
                if keep_genome_info:
                    relevant_cols.extend(
                        [
                            "Mutation genome position GRCh37",
                            "GENOMIC_WT_ALLELE_SEQ",
                            "GENOMIC_MUT_ALLELE_SEQ",
                            "GENOMIC_MUTATION_ID",
                        ]
                    )  # * uncomment to include strand information (tested not to be accurate for CMC)
                    # relevant_cols.extend(['Mutation genome position GRCh37', 'GENOMIC_MUTATION_ID'])    #* erase to include strand information (tested not to be accurate for CMC)
            else:
                relevant_cols = [
                    "GENE_SYMBOL",
                    "TRANSCRIPT_ACCESSION",
                    "MUTATION_ID",
                    "MUTATION_CDS",
                    "MUTATION_AA",
                ]
                if keep_genome_info:
                    relevant_cols.extend(["HGVSG", "STRAND", "GENOMIC_MUTATION_ID"])

            df = pd.read_csv(mutation_tsv_file, usecols=relevant_cols, sep="\t")

            # Get seq_ID and mutation columns
            if mutation_class == "cancer" or mutation_class == "cancer_example":
                df["MUTATION_URL"] = df["MUTATION_URL"].str.extract(r"id=(\d+)")
                df = df.rename(
                    columns={
                        "ACCESSION_NUMBER": "seq_ID",
                        "Mutation CDS": "mutation",
                        "MUTATION_URL": "MUTATION_ID",
                        "Mutation AA": "mutation_aa",
                    }
                )

                if keep_genome_info:
                    # * erase to include strand information (tested not to be accurate for CMC)
                    # df = df.rename(
                    #     columns={
                    #         "Mutation genome position GRCh37": "position_genome",
                    #     }
                    # )

                    from gget.gget_mutate import (
                        mutation_pattern,
                        convert_chromosome_value_to_int_when_possible,
                    )
                    import numpy as np

                    # * uncomment to include strand information (tested not to be accurate for CMC)
                    df[["chromosome", "GENOME_POS"]] = df[
                        "Mutation genome position GRCh37"
                    ].str.split(":", expand=True)
                    df["chromosome"] = df["chromosome"].apply(
                        convert_chromosome_value_to_int_when_possible
                    )
                    df[["GENOME_START", "GENOME_STOP"]] = df["GENOME_POS"].str.split(
                        "-", expand=True
                    )

                    df[["nucleotide_positions", "actual_mutation"]] = df[
                        "mutation"
                    ].str.extract(mutation_pattern)

                    sub_mask = df["actual_mutation"].str.contains(">")
                    ins_mask = (df["actual_mutation"].str.contains("ins")) & (
                        ~df["actual_mutation"].str.contains("delins")
                    )
                    delins_mask = df["actual_mutation"].str.contains("delins")
                    ins_delins_mask = ins_mask | delins_mask
                    sub_ins_delins_mask = sub_mask | ins_delins_mask

                    df.loc[sub_mask, "wt_allele_cds"] = (
                        df.loc[sub_mask, "actual_mutation"].str.split(">").str[0]
                    )
                    df.loc[sub_mask, "mut_allele_cds"] = (
                        df.loc[sub_mask, "actual_mutation"].str.split(">").str[1]
                    )

                    df.loc[ins_delins_mask, "mut_allele_cds"] = df.loc[
                        ins_delins_mask, "actual_mutation"
                    ].str.extract(r"ins(.+)")[0]

                    df["strand"] = np.nan

                    df.loc[sub_ins_delins_mask, "strand"] = np.where(
                        pd.isna(df.loc[sub_ins_delins_mask, "GENOMIC_MUT_ALLELE_SEQ"]),
                        np.nan,
                        np.where(
                            df.loc[sub_ins_delins_mask, "mut_allele_cds"]
                            != df.loc[sub_ins_delins_mask, "GENOMIC_MUT_ALLELE_SEQ"],
                            "-",
                            "+",
                        ),
                    )

                    df.loc[sub_mask, "actual_mutation_updated"] = (
                        df.loc[sub_mask, "GENOMIC_WT_ALLELE_SEQ"]
                        + ">"
                        + df.loc[sub_mask, "GENOMIC_MUT_ALLELE_SEQ"]
                    )
                    df.loc[ins_mask, "actual_mutation_updated"] = (
                        "ins" + df.loc[ins_mask, "GENOMIC_MUT_ALLELE_SEQ"]
                    )
                    df.loc[delins_mask, "actual_mutation_updated"] = (
                        "delins" + df.loc[delins_mask, "GENOMIC_MUT_ALLELE_SEQ"]
                    )

                    df.loc[~sub_ins_delins_mask, "actual_mutation_final"] = df.loc[
                        ~sub_ins_delins_mask, "actual_mutation"
                    ]

                    df.loc[sub_ins_delins_mask, "actual_mutation_final"] = np.where(
                        pd.isna(df.loc[sub_ins_delins_mask, "strand"]),
                        np.nan,
                        np.where(
                            df.loc[sub_ins_delins_mask, "strand"] == "+",
                            df.loc[sub_ins_delins_mask, "actual_mutation"],
                            df.loc[sub_ins_delins_mask, "actual_mutation_updated"],
                        ),
                    )

                    df["mutation_genome"] = np.where(
                        df["GENOME_START"] != df["GENOME_STOP"],
                        "g."
                        + df["GENOME_START"].astype(str)
                        + "_"
                        + df["GENOME_STOP"].astype(str)
                        + df["actual_mutation_final"],
                        "g."
                        + df["GENOME_START"].astype(str)
                        + df["actual_mutation_final"],
                    )

                    df.loc[
                        df["Mutation genome position GRCh37"].isna(), "mutation_genome"
                    ] = np.nan

                    df.drop(
                        columns=[
                            "GENOME_POS",
                            "GENOME_START",
                            "GENOME_STOP",
                            "nucleotide_positions",
                            "actual_mutation",
                            "actual_mutation_updated",
                            "actual_mutation_final",
                            "Mutation genome position GRCh37",
                            "wt_allele_cds",
                            "mut_allele_cds",
                            "GENOMIC_WT_ALLELE_SEQ",
                            "GENOMIC_MUT_ALLELE_SEQ",
                        ],
                        inplace=True,
                    )

            else:
                df = df.rename(
                    columns={
                        "GENE_SYMBOL": "GENE_NAME",
                        "TRANSCRIPT_ACCESSION": "seq_ID",
                        "MUTATION_CDS": "mutation",
                        "MUTATION_AA": "mutation_aa",
                    }
                )

                if keep_genome_info:
                    df["mutation_genome"] = df["HGVSG"].str.split(":").str[1]

                    df.drop(columns=["HGVSG"], inplace=True)

                    df = df.rename(
                        columns={
                            "CHROMOSOME": "chromosome",
                            "STRAND": "strand",
                        }
                    )

            # Remove version numbers from Ensembl IDs
            df["seq_ID"] = df["seq_ID"].str.split(".").str[0]

            df["gene_name"] = df["GENE_NAME"].astype(str)
            df["mutation_id"] = df["MUTATION_ID"].astype(str)

            df = df.drop(columns=["GENE_NAME", "MUTATION_ID"])

            if remove_duplicates:
                duplicate_count = (
                    df.duplicated(subset=["seq_ID", "mutation"], keep=False).sum() // 2
                )
                print(
                    f"Removing {duplicate_count} duplicate entries from the COSMIC csv for gget mutate: {duplicate_count}"
                )
                df["non_na_count"] = df.notna().sum(axis=1)
                df = df.sort_values(by="non_na_count", ascending=False)
                df = df.drop_duplicates(subset=["seq_ID", "mutation"], keep="first")
                df = df.drop(columns=["non_na_count"])

            mutate_csv_out = mutation_tsv_file.replace(".tsv", "_gget_mutate.csv")
            df.to_csv(mutate_csv_out, index=False)

            if verbose:
                logger.info(
                    f"Modified mutations file for use with gget mutate created at {mutate_csv_out}"
                )

    else:
        # Check if 'entity' argument is valid
        sps = [
            "mutations",
            "pubmed",
            "genes",
            "studies",
            "samples",
            "cancer",
            "tumour_site",
        ]
        if entity not in sps:
            raise ValueError(
                f"'entity' argument specified as {entity}. Expected one of: {', '.join(sps)}"
            )

        # Translate categories to match COSMIC data table IDs
        if entity == "cancer":
            entity = "disease"

        if entity == "tumour_site":
            entity = "tumour"

        r = requests.get(
            url=COSMIC_GET_URL + entity + "?q=" + searchterm + "&export=json"
        )

        # Check if the request returned an error (e.g. gene not found)
        if not r.ok:
            raise RuntimeError(
                f"COSMIC API request returned error {r.status_code}. "
                "Please double-check the arguments and try again.\n"
            )

        if r.text == "\n":
            logger.warning(
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
                        parsing_mutations[1]
                        .replace('" ', "")
                        .replace('"', "")
                        .capitalize()
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
            if out:
                # Create saving directory
                directory = "/".join(out.split("/")[:-1])
                if directory != "":
                    os.makedirs(directory, exist_ok=True)

                json_out = os.path.join(out, f"gget_cosmic_{entity}_{searchterm}.json")
                with open(json_out, "w", encoding="utf-8") as f:
                    json_package.dump(results_dict, f, ensure_ascii=False, indent=4)

            else:
                return results_dict

        else:
            if out:
                # Create saving directory
                directory = "/".join(out.split("/")[:-1])
                if directory != "":
                    os.makedirs(directory, exist_ok=True)

                df_out = os.path.join(out, f"gget_cosmic_{entity}_{searchterm}.csv")
                corr_df.to_csv(df_out, index=False)

            else:
                return corr_df
