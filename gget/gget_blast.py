import pandas as pd
import json as json_package
import time
from bs4 import BeautifulSoup
import logging

# Add and format time stamp in logging messages
logging.basicConfig(
    format="%(asctime)s %(levelname)s %(message)s",
    level=logging.INFO,
    datefmt="%c",
)
# Using urllib instead of requests here because requests does not
# support long queries (queries very long here due to input sequence)
from urllib.request import urlopen, Request
from urllib.parse import urlencode

# Custom functions
from .utils import parse_blast_ref_page, wrap_cols_func

# Constants
from .constants import (
    BLAST_URL,
    BLAST_CLIENT,
)


def blast(
    sequence,
    program="default",
    database="default",
    limit=50,
    expect=10.0,
    low_comp_filt=False,
    megablast=True,
    verbose=True,
    wrap_text=False,
    json=False,
    save=False,
):
    """
    BLAST a nucleotide or amino acid sequence against any BLAST DB.
    Args:
     - sequence       Sequence (str) or path to FASTA file.
                      (If more than one sequence in FASTA file, only the first will be submitted to BLAST.)
     - program        'blastn', 'blastp', 'blastx', 'tblastn', or 'tblastx'.
                      Default: 'blastn' for nucleotide sequences; 'blastp' for amino acid sequences.
     - database       'nt', 'nr', 'refseq_rna', 'refseq_protein', 'swissprot', 'pdbaa', or 'pdbnt'.
                      Default: 'nt' for nucleotide sequences; 'nr' for amino acid sequences.
                      More info on BLAST databases: https://ncbi.github.io/blast-cloud/blastdb/available-blastdbs.html
     - limit          Limits number of hits to return. Default 50.
     - expect         float or None. An expect value cutoff. Default 10.0.
     - low_comp_filt  True/False whether to apply low complexity filter. Default False.
     - megablast      True/False whether to use the MegaBLAST algorithm (blastn only). Default True.
     - verbose        True/False whether to print progress information. Default True.
     - wrap_text      If True, displays data frame with wrapped text for easy reading. Default: False.
     - json           If True, returns results in json format instead of data frame. Default: False.
     - save           If True, the data frame is saved as a csv in the current directory (default: False).

    Returns a data frame with the BLAST results.

    NCBI server rule:
    Run scripts weekends or between 9 pm and 5 am Eastern time
    on weekdays if more than 50 searches will be submitted.

    Note: This function does not check the validity of the arguments
    and passes the values to the server as is.
    """
    # Server rules:
    # 1. Do not contact the server more often than once every 10 seconds.
    # 2. Do not poll for any single RID more often than once a minute.
    # 3. Use the URL parameter email and tool, so that the NCBI
    #    can contact you if there is a problem.
    # 4. Run scripts weekends or between 9 pm and 5 am Eastern time
    #    on weekdays if more than 50 searches will be submitted.
    # Reference: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=DeveloperInfo

    # Please note that NCBI uses the new Common URL API for BLAST searches
    # on the internet (http://ncbi.github.io/blast-cloud/dev/api.html). Thus,
    # some of the arguments used by this function are not (or are no longer)
    # officially supported by NCBI. Although they are still functioning, this
    # may change in the future.

    # Define server URL and client
    url = BLAST_URL
    client = BLAST_CLIENT

    ## Clean up arguments
    # If the path to a fasta file was provided instead of a nucleotide sequence,
    # read the file and extract the first sequence
    if "." in sequence:
        if ".txt" in sequence:
            # Read the text file
            titles = []
            seqs = []
            with open(sequence) as text_file:
                for i, line in enumerate(text_file):
                    # Recognize a title line by the '>' character
                    if line[0] == ">":
                        # Append title line to titles list
                        titles.append(line.strip())
                    else:
                        seqs.append(line.strip())

        elif ".fa" in sequence:
            # Read the FASTA
            titles = []
            seqs = []
            with open(sequence) as fasta_file:
                for i, line in enumerate(fasta_file):
                    # Each second line will be a title line
                    if i % 2 == 0:
                        if line[0] != ">":
                            raise ValueError(
                                "Expected FASTA to start with a '>' character. "
                            )
                        else:
                            # Append title line to titles list
                            titles.append(line.strip())
                    else:
                        if line[0] == ">":
                            raise ValueError(
                                "FASTA contains two lines starting with '>' in a row -> missing sequence line. "
                            )
                        # Append sequences line to seqs list
                        else:
                            seqs.append(line.strip())
        else:
            raise ValueError(
                "File format not recognized. gget BLAST currently only supports '.txt' or '.fa' files. "
            )

        # Set the first sequence from the fasta file as 'sequence'
        sequence = seqs[0]
        if len(seqs) > 1:
            logging.warning(
                "File contains more than one sequence. Only the first sequence will be submitted to BLAST."
            )

    ## Set program and database

    # Convert program and database to lower case
    program = program.lower()
    database = database.lower()
    # Valid program and database options
    programs = ["blastn", "blastp", "blastx", "tblastn", "tblastx"]
    dbs = ["nt", "nr", "refseq_rna", "refseq_protein", "swissprot", "pdbaa", "pdbnt"]

    # If user does not specify the program,
    # check if a nulceotide or amino acid sequence was passed
    if program == "default":
        # Set of all possible nucleotides and amino acids
        nucleotides = set("ATGC")
        amino_acids = set("ARNDCQEGHILKMFPSTWYVBZ")

        # If sequence is a nucleotide sequence, set program to blastn
        if set(sequence) <= nucleotides:
            program = "blastn"

            # Set database to nt (unless user specified another database)
            if database == "default":
                database = "nt"
                if verbose:
                    logging.info("Sequence recognized as nucleotide sequence.")
                    logging.info("BLAST will use program 'blastn' with database 'nt'.")
            else:
                # Check if the user specified database is valid
                if database not in dbs:
                    raise ValueError(
                        f"Database specified is {database}. Expected one of: {', '.join(dbs)}"
                    )

                else:
                    if verbose:
                        logging.info("Sequence recognized as nucleotide sequence.")
                        logging.info(
                            "BLAST will use program 'blastn' with user-specified database."
                        )
        # If sequence is an amino acid sequence, set program to blastp
        elif set(sequence) <= amino_acids:
            program = "blastp"

            # Set database to nr (unless user specified another database)
            if database == "default":
                database = "nr"
                if verbose:
                    logging.info("Sequence recognized as amino acid sequence.")
                    logging.info("BLAST will use program 'blastp' with database 'nr'.")
            else:
                # Check if the user specified database is valid
                if database not in dbs:
                    raise ValueError(
                        f"Database specified is {database}. Expected one of: {', '.join(dbs)}"
                    )

                else:
                    if verbose:
                        logging.info("Sequence recognized as amino acid sequence.")
                        logging.info(
                            "BLAST will use program 'blastp' with user-specified database."
                        )
        else:
            raise ValueError(
                f"""
                Sequence not automatically recognized as a nucleotide or amino acid sequence.
                Please specify 'program' and 'database'.
                Program options: {', '.join(programs)} 
                Database options:  {', '.join(dbs)} 
                """
            )

    else:
        # Check if the user specified program is valid
        if program not in programs:
            raise ValueError(
                f"Program specified is {program}. Expected one of: {', '.join(programs)}"
            )

        # Ask user to also specify database
        if database == "default":
            raise ValueError(
                f"""
                User-specified program requires user-specified database. Please also specify argument 'database'. 
                Database options:  {', '.join(dbs)}
                """
            )
        else:
            # Check if the user specified database is valid
            if database not in dbs:
                raise ValueError(
                    f"Database specified is {database}. Expected one of: {', '.join(dbs)}"
                )

    ## Translate filter arguments
    if low_comp_filt is False:
        low_comp_filt = None
    else:
        low_comp_filt = "T"

    if megablast is False:
        megablast = None
    else:
        megablast = "on"

    ## Submit search
    #  The following code was partly adapted from the Biopython BLAST NCBIWWW project written
    #  by Jeffrey Chang (Copyright 1999), Brad Chapman, and Chris Wroe distributed under the
    #  Biopython License Agreement and BSD 3-Clause License
    #  https://github.com/biopython/biopython/blob/171697883aca6894f8367f8f20f1463ce7784d0c/LICENSE.rst

    # Args for the PUT command
    put_args = [
        ("PROGRAM", program),
        ("DATABASE", database),
        ("QUERY", sequence),
        ("DESCRIPTIONS", limit),
        ("HITLIST_SIZE", limit),
        ("ALIGNMENTS", 0),
        ("EXPECT", expect),
        ("FILTER", low_comp_filt),
        ("MEGABLAST", megablast),
        ("CMD", "Put"),
    ]

    # Define query
    put_query = [x for x in put_args if x[1] is not None]
    put_message = urlencode(put_query).encode()

    # Submit search to server
    request = Request(url, put_message, {"User-Agent": client})
    handle = urlopen(request)

    ## Fetch Request ID (RID) and estimated time to completion (RTOE)
    RID, RTOE = parse_blast_ref_page(handle)

    # Wait for search to complete
    # (At least 11 seconds to comply with server rule 1)
    if RTOE < 11:
        # Communicate RTOE
        if verbose:
            logging.info(f"BLAST initiated. Estimated time to completion: 11 seconds.")
        time.sleep(11)
    else:
        # Communicate RTOE
        if verbose:
            logging.info(
                f"BLAST initiated with search ID {RID}. Estimated time to completion: {RTOE} seconds."
            )
        time.sleep(int(RTOE))

    ## Poll server for status and fetch search results
    # Args for the GET command
    get_args = [
        ("RID", RID),
        ("DESCRIPTIONS", limit),
        ("HITLIST_SIZE", limit),
        ("ALIGNMENTS", 0),
        ("FORMAT_TYPE", "HTML"),
        ("CMD", "Get"),
    ]
    get_query = [x for x in get_args if x[1] is not None]
    get_message = urlencode(get_query).encode()

    ## Poll NCBI until the results are ready
    searching = True
    i = 0
    while searching:
        if i > 0:
            # Sleep for 61 seconds if first fetch was not succesful
            # to comply with server rules
            time.sleep(61)

        # Query for search status
        request = Request(url, get_message, {"User-Agent": client})
        handle = urlopen(request)
        results = handle.read().decode()

        # Fetch search status
        i = results.index("Status=")
        j = results.index("\n", i)
        status = results[i + len("Status=") : j].strip()

        if status == "WAITING":
            if verbose:
                logging.info("BLASTING...")
            i += 1
            continue

        elif status == "FAILED":
            logging.error(
                f"Search {RID} failed; please try again and/or report to blast-help@ncbi.nlm.nih.gov."
            )
            return

        elif status == "UNKNOWN":
            logging.error(f"NCBI status {status}. Search {RID} expired.")
            return

        elif status == "READY":
            if verbose:
                logging.info("Retrieving results...")
            # Stop search
            searching = False

            ## Return results
            # Parse HTML results
            soup = BeautifulSoup(results, "html.parser")
            # Get the descriptions table
            dsc_table = soup.find(
                lambda tag: tag.name == "table"
                and tag.has_attr("id")
                and tag["id"] == "dscTable"
            )

            if dsc_table is None:
                logging.error(
                    f"No significant similarity found for search {RID}. If your sequence is very short, try increasing the 'expect' argument."
                )
                return

            results_df = pd.read_html(str(dsc_table))[0]
            # Drop the first column
            results_df = results_df.iloc[:, 1:]

            if wrap_text:
                df_wrapped = results_df.copy()
                wrap_cols_func(df_wrapped, ["Description"])

            if json:
                results_dict = json_package.loads(results_df.to_json(orient="records"))
                if save:
                    with open("gget_blast_results.json", "w", encoding="utf-8") as f:
                        json_package.dump(results_dict, f, ensure_ascii=False, indent=4)

                return results_dict

            else:
                # Save
                if save:
                    results_df.to_csv("gget_blast_results.csv", index=False)

                return results_df

        else:
            logging.error(
                f"Something unexpected happened. Search {RID} possibly failed; please try again and/or report to blast-help@ncbi.nlm.nih.gov"
            )
            return
