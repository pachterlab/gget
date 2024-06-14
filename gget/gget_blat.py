import json as json_package
from json.decoder import JSONDecodeError
import pandas as pd
from urllib.request import urlopen

from .utils import set_up_logger, read_fasta

logger = set_up_logger()


def blat(
    sequence,
    seqtype="default",
    assembly="human",
    json=False,
    save=False,
    verbose=True,
):
    """
    BLAT a nucleotide or amino acid sequence against any BLAT UCSC assembly.

    Args:
     - sequence       Sequence (str) or path to fasta file containing one sequence.
     - seqtype        'DNA', 'protein', 'translated%20RNA', or 'translated%20DNA'.
                      Default: 'DNA' for nucleotide sequences; 'protein' for amino acid sequences.
     - assembly       'human' (hg38) (default), 'mouse' (mm39), 'zebrafinch' (taeGut2),
                      or any of the species assemblies available at https://genome.ucsc.edu/cgi-bin/hgBlat
                      (use short assembly name as listed after the "/").
     - json           If True, returns results in json format instead of data frame. Default: False.
     - save           If True, the data frame is saved as a csv in the current directory (default: False).
     - verbose        True/False whether to print progress information. Default True.

    Returns a data frame with the BLAT results.
    """

    ## Clean up sequence
    # If the path to a fasta file was provided instead of a nucleotide sequence,
    # read the file and extract the first sequence
    if "." in sequence:
        if ".txt" in sequence or ".fa" in sequence:
            _, seqs = read_fasta(sequence)

        else:
            raise ValueError(
                "File format not recognized. gget BLAT currently only supports '.txt' or '.fa' files. "
            )

        # Set the first sequence from the fasta file as 'sequence'
        sequence = seqs[0]
        if len(seqs) > 1:
            if verbose:
                logger.info(
                    "File contains more than one sequence. Only the first sequence will be submitted to BLAT."
                )

    # Shorten sequence to length limit if necessary
    if len(sequence) > 8000:
        if verbose:
            logger.info(
                "Length of sequence is > 8000. Only the fist 8000 characters will be submitted to BLAT."
            )
        sequence = sequence[:8000]

    # Convert sequence to upper case
    sequence = sequence.upper()

    ## Set seqtype
    # Valid seqtype options
    seqtypes = ["DNA", "protein", "translated%20RNA", "translated%20DNA"]

    # If user does not specify the seqtype,
    # check if a nucleotide or amino acid sequence was passed
    if seqtype == "default":
        # Set of all possible nucleotides and amino acids
        nucleotides = set("ATGCN")
        amino_acids = set("ARNDCQEGHILKMFPSTWYVBZXBJZ")

        # If sequence is a nucleotide sequence, set seqtype to DNA
        if set(sequence) <= nucleotides:
            seqtype = "DNA"
            if verbose:
                logger.info(
                    f"Sequence recognized as nucleotide sequence. 'seqtype' will be set as {seqtype}."
                )

        # If sequence is an amino acid sequence, set seqtype to protein
        elif set(sequence) <= amino_acids:
            seqtype = "protein"
            if verbose:
                logger.info(
                    f"Sequence recognized as amino acid sequence. 'seqtype' will be set as {seqtype}."
                )

        else:
            raise ValueError(
                f"""
                Sequence not automatically recognized as a nucleotide or amino acid sequence.
                Please specify 'seqtype'.
                Seqtype options: {', '.join(seqtypes)} 
                """
            )

    else:
        # Check if the user specified seqtype is valid
        if seqtype not in seqtypes:
            raise ValueError(
                f"Seqtype specified is {seqtype}. Expected one of {', '.join(seqtypes)}"
            )

    ## Set assembly
    # Note: If assembly not found, defaults to hg38
    if assembly == "human" or assembly == "homo_sapiens":
        database = "hg38"
    elif assembly == "mouse" or assembly == "mus_musculus":
        database = "mm39"
    elif assembly == "zebrafinch" or assembly == "taeniopygia_guttata":
        database = "taeGut2"
    else:
        database = assembly

    # Define server URL
    url = f"https://genome.ucsc.edu/cgi-bin/hgBlat?userSeq={sequence}&type={seqtype}&db={database}&output=json"

    # Submit URL request
    r = urlopen(url)
    if r.status != 200:
        raise RuntimeError(
            f"HTTP response status code {r.status}. "
            "Please double-check arguments and try again.\n"
        )

    try:
        # Read json results into a dictionary
        results = json_package.load(r)
    except JSONDecodeError:
        logger.error(
            f"""
            BLAT of seqtype '{seqtype}' using assembly '{database}' was unsuccesful. 
            Possible causes: 
            - Sequence possibly too short (required minimum: 20 characters). 
            - Assembly possibly invalid. All available species with their respective assemblies are listed at https://genome.ucsc.edu/cgi-bin/hgBlat.
            """
        )
        return

    if len(results["blat"]) == 0:
        if verbose:
            logger.info(
                f"No {seqtype} BLAT matches were found for this sequence in genome {results['genome']}."
            )
        return

    # Let user know if assembly was not found
    # If this is the case, BLAT automatically defaults to human (hg38)
    if results["genome"] != database:
        logger.warning(
            f"Assembly {database} not recognized. Defaulted to {results['genome']} instead."
        )

    ## Build data frame to resemble BLAT web search results
    # Define dataframe from dictionary
    df_dict = {}

    for field in results["fields"]:
        df_dict.update({field: []})

    for blat_result_list in results["blat"]:
        for field, (i, result) in zip(results["fields"], enumerate(blat_result_list)):
            df_dict[field].append(result)

    df = pd.DataFrame(df_dict)

    # Calculate % aligned sequence of submitted sequence
    aligned_size = df["qEnd"] - df["qStart"]
    df["%_aligned"] = round((100 / df["qSize"]) * aligned_size, 2)
    # Calculate % matched sequence of aligned sequence
    df["%_matched"] = round((100 / aligned_size) * df["matches"], 2)
    # Add genome column
    df["genome"] = results["genome"]

    # Adjust sequence start to match website
    df["qStart"] = df["qStart"] + 1
    df["tStart"] = df["tStart"] + 1

    # Rename columns
    columns_dict = {
        "misMatches": "mismatches",
        "qName": "query",
        "qSize": "query_size",
        "qStart": "aligned_start",
        "qEnd": "aligned_end",
        "tName": "chromosome",
        "tStart": "start",
        "tEnd": "end",
    }

    df = df.rename(columns=columns_dict)

    # Change columns order (this also drops all unmentioned columns)
    df = df.reindex(
        columns=[
            "genome",
            "query_size",
            "aligned_start",
            "aligned_end",
            "matches",
            "mismatches",
            "%_aligned",
            "%_matched",
            "chromosome",
            "strand",
            "start",
            "end",
        ]
    )

    if json:
        results_dict = json_package.loads(df.to_json(orient="records"))
        if save:
            with open("gget_blat_results.json", "w", encoding="utf-8") as f:
                json_package.dump(results_dict, f, ensure_ascii=False, indent=4)

        return results_dict

    else:
        if save:
            df.to_csv("gget_blat_results.csv", index=False)

        return df
