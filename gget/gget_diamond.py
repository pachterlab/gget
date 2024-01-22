import logging
import subprocess
import sys
import platform
import os
import pandas as pd
import uuid
import json as json_package

# Add and format time stamp in logging messages
logging.basicConfig(
    format="%(asctime)s %(levelname)s %(message)s",
    level=logging.INFO,
    datefmt="%c",
)
# Mute numexpr threads info
logging.getLogger("numexpr").setLevel(logging.WARNING)

from .compile import PACKAGE_PATH
from .utils import tsv_to_df, create_tmp_fasta, remove_temp_files

# Path to precompiled diamond binary
if platform.system() == "Windows":
    PRECOMPILED_DIAMOND_PATH = os.path.join(
        PACKAGE_PATH, f"bins/{platform.system()}/diamond.exe"
    )
else:
    PRECOMPILED_DIAMOND_PATH = os.path.join(
        PACKAGE_PATH, f"bins/{platform.system()}/diamond"
    )


def diamond(
    query,
    reference,
    diamond_db=None,
    sensitivity="very-sensitive",
    threads=1,
    diamond_binary=None,
    verbose=True,
    json=False,
    out=None,
):
    """
    Align multiple protein or translated DNA sequences using DIAMOND (https://www.nature.com/articles/nmeth.3176).

    Args:
    - query          Sequences (str or list) or path to FASTA file containing sequences to be aligned against the reference.
    - reference      Reference sequences (str or list) or path to FASTA file containing reference sequences.
    - diamond_db     Path to save DIAMOND database created from reference.
                     Default: None -> Temporary db file will be deleted after alignment or saved in 'out' if 'out' is provided.
    - sensitivity    Sensitivity of DIAMOND alignment.
                     One of the following: fast, mid-sensitive, sensitive, more-sensitive, very-sensitive or ultra-sensitive.
                     Default: "very-sensitive"
    - threads        Number of threads to use for alignment. Default: 1.
    - diamond_binary Path to DIAMOND binary, e.g. path/bins/Linux/diamond. Default: None -> Uses DIAMOND binary installed with gget.
    - verbose        True/False whether to print progress information. Default True.
    - json           If True, returns results in json format instead of data frame. Default: False.
    - out            Path to folder to save DIAMOND results in. Default: Standard out, temporary files are deleted.

    Returns a data frame with the DIAMOND alignment results. (Or JSON formatted dictionary if json=True.)
    """
    # Check argument validity
    supported_sens = [
        "fast",
        "mid-sensitive",
        "sensitive",
        "more-sensitive",
        "very-sensitive",
        "ultra-sensitive",
    ]
    if sensitivity not in supported_sens:
        raise ValueError(
            f"'sensitivity' argument specified as {sensitivity}. Expected one of: {', '.join(supported_sens)}"
        )

    # Handle command line passing path to FASTA as a list
    if isinstance(query, list) and len(query) == 1:
        query = query[0]
    if isinstance(reference, list) and len(reference) == 1:
        reference = reference[0]

    # Define paths to query/reference/db/output files
    files_to_delete = []
    if "." in query:
        input_file = os.path.abspath(query)
    else:
        input_file = create_tmp_fasta(query)
        files_to_delete.append(input_file)

    if "." in reference:
        reference_file = os.path.abspath(reference)
    else:
        reference_file = create_tmp_fasta(reference)
        files_to_delete.append(reference_file)

    # Create out folder if it does not exist
    if out:
        os.makedirs(out, exist_ok=True)

    if out:
        out = os.path.abspath(out)
        output = f"{out}/DIAMOND_results.tsv"
    else:
        output = os.path.abspath(f"tmp_{str(uuid.uuid4())}_out.tsv")
        files_to_delete.append(output)

    if not diamond_db and out:
        diamond_db = f"{out}/DIAMOND_db"
    elif not diamond_db:
        diamond_db = os.path.abspath(f"tmp_db_{str(uuid.uuid4())}")
        files_to_delete.append(diamond_db + ".dmnd")

    if diamond_binary:
        DIAMOND = diamond_binary
    else:
        DIAMOND = PRECOMPILED_DIAMOND_PATH

    # Replace slashes in paths for Windows compatibility
    if platform.system() == "Windows":
        DIAMOND_w = DIAMOND.replace("/", "\\")
        reference_file_w = reference_file.replace("/", "\\")
        diamond_db_w = diamond_db.replace("/", "\\")
        input_file_w = input_file.replace("/", "\\")
        reference_file_w = reference_file.replace("/", "\\")
        output_w = output.replace("/", "\\")

    if platform.system() == "Windows":
        command = f"{DIAMOND} version \
        && {DIAMOND_w} makedb --quiet --in {reference_file_w} --db {diamond_db_w} --threads {threads} \
        && {DIAMOND_w} blastp --outfmt 6 qseqid sseqid pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore \
            --quiet --query {input_file_w} --db {reference_file_w} --out {output_w} --{sensitivity} --threads {threads} --ignore-warnings"
    else:
        command = f"'{DIAMOND}' version \
        && '{DIAMOND}' makedb --quiet --in '{reference_file}' --db '{diamond_db}' --threads {threads} \
        && '{DIAMOND}' blastp --outfmt 6 qseqid sseqid pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore \
            --quiet --query '{input_file}' --db '{reference_file}' --out '{output}' --{sensitivity} --threads {threads} --ignore-warnings"

    # Run DIAMOND
    if verbose:
        logging.info(f"Creating DIAMOND database and initiating alignment...")

    with subprocess.Popen(command, shell=True, stderr=subprocess.PIPE) as process:
        stderr = process.stderr.read().decode("utf-8")
        # Log the standard error if it is not empty
        if stderr:
            sys.stderr.write(stderr)

    # Exit system if the subprocess returned wstdout = sys.stdout
    if process.wait() != 0:
        raise RuntimeError("DIAMOND alignment failed.")
    else:
        if verbose:
            logging.info(f"DIAMOND alignment complete.")

    df_diamond = tsv_to_df(
        output,
        headers=[
            "query_accession",
            "subject_accession",
            "identity_percentage",
            "query_seq_length",
            "subject_seq_length",
            "length",
            "mismatches",
            "gap_openings",
            "query_start",
            "query_end",
            "subject_start",
            "subject_end",
            "e-value",
            "bit_score",
        ],
    )

    # Delete temporary files
    if files_to_delete:
        remove_temp_files(files_to_delete)

    if json:
        results_dict = json_package.loads(df_diamond.to_json(orient="records"))
        if out:
            with open(f"{out}/gget_diamond_results.json", "w", encoding="utf-8") as f:
                json_package.dump(results_dict, f, ensure_ascii=False, indent=4)

        return results_dict

    else:
        if out:
            df_diamond.to_csv(f"{out}/gget_diamond_results.csv", index=False)
        return df_diamond
