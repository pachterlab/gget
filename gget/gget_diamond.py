import logging
import subprocess
import sys
import platform
import os
import pandas as pd

from .compile import PACKAGE_PATH

# Path to precompiled diamond binary
if platform.system() == "Windows":
    PRECOMPILED_DIAMOND_PATH = os.path.join(
        PACKAGE_PATH, f"bins/{platform.system()}/diamond.exe"
    )
else:
    PRECOMPILED_DIAMOND_PATH = os.path.join(
        PACKAGE_PATH, f"bins/{platform.system()}/diamond"
    )

from .gget_setup import (
    ELM_FILES, 
    ELM_INSTANCES_FASTA,
    ELM_CLASSES_TSV,
    ELM_INSTANCES_TSV
)

from .constants import RANDOM_ID



def diamond(input, reference, json=False, verbose=True, out=None, sensitivity="very-sensitive"):
    """
    Perform protein sequence alignment using DIAMOND for multiple sequences

    Args:
     - input          Input sequences path and file name (include ,fa) in FASTA file format
     - reference      Reference file path and file name (include .fa) in FASTA file format
     - json           If True, returns results in json format instead of data frame. Default: False.
     - out            folder name to save two resulting csv files. Default: results (default: None).
     - verbose        True/False whether to print progress information. Default True.
     - sensitivity    The sensitivity can be adjusted using the options --fast, --mid-sensitive, --sensitive, --more-sensitive, --very-sensitive and --ultra-sensitive.

    Returns DIAMOND output in tsv format 
    """
    #TODO: --very_sensitive and makedb --in as args
    # if out is None, create temp file and delete once get dataframe
    # if make
    if out is None:
        command = f"diamond makedb --in {reference} -d reference && diamond blastp -q {input} -d reference -o tmp_out.tsv --{sensitivity}"
    else:
        # The double-quotation marks allow white spaces in the path, but this does not work for Windows
        command = f"diamond makedb --in {reference} -d reference && diamond blastp -q {input} -d reference -o {out}.tsv --{sensitivity}"
    # Run diamond command and write command output
    with subprocess.Popen(command, shell=True, stderr=subprocess.PIPE) as process_2:
        stderr_2 = process_2.stderr.read().decode("utf-8")
        # Log the standard error if it is not empty
        if stderr_2:
            sys.stderr.write(stderr_2)
    # Exit system if the subprocess returned wstdout = sys.stdout

    if process_2.wait() != 0:
        logging.error(
            """
            DIAMOND failed. Please check that you have a diamond executable file for Windows or Linux in the bin folder.
            """
        )
        return
    else:
        logging.info(
            f"DIAMOND run complete."
        )
       
              
    # df_diamond = tsv_to_df("diamond_out.tsv", ["query_accession", "target_accession", "Per. Ident" , "length", "mismatches", "gap_openings", "query_start", "query_end", "target_start", "target_end", "e-value", "bit_score"])
    # return df_diamond
   

              