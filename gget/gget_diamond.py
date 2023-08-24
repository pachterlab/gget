import logging
import subprocess
import sys
import platform
import os
import pandas as pd
import uuid

# DIAMOND and ELM id for temporary files
RANDOM_ID = str(uuid.uuid4())


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


def tsv_to_df(tsv_file, headers = None):
    """
    Convert tsv file to dataframe format

    Args:
    tsv_file - file to be converted 

    Returns:
    df -  dataframe
    
    """
    
    try:
        df = pd.DataFrame()
        if headers:
            df = pd.read_csv(tsv_file, sep="\t", names=headers)
        else:
            # ELM Instances.tsv file had 5 lines before headers and data
            df = pd.read_csv(tsv_file, sep="\t", skiprows=5)
        return df
    

    except pd.errors.EmptyDataError:
        logging.warning(f"Query did not result in any matches.")
        return None

def create_input_file(sequences):
    """
    Copy sequences to a temporary fasta file for DIAMOND alignment

    Args:
    sequences - list of user input amino acid sequences 

    Returns: input file absolute path
    """
    print(sequences)
    with open(f"tmp_{RANDOM_ID}.fa", 'w') as f:
        for idx, seq in enumerate(sequences):
            f.write(f'>Seq {idx}\n{seq}')

    return f"tmp_{RANDOM_ID}.fa"
    # check if correct sequences are written to file
    # try:
    #     with open(f"{os.getcwd()}tmp_{RANDOM_ID}.fa", 'r') as f:
    #         print(f.read())
    # except:
    #     continue

def remove_temp_files():
    """
    Delete temporary files

    Args:
    input       - Input fasta file containing amino acid sequences
    out         - Output tsv file containing the output returned by DIAMOND
    reference   - Reference database binary file produced by DIAMOND

    Returns: 
    None 
    """
    if os.path.exists("tmp_out.tsv"):
        os.remove("tmp_out.tsv")
    if os.path.exists("reference.dmnd"):
        os.remove("reference.dmnd")
    if os.path.exists("tmp_{RANDOM_ID}.fa"):
        os.remove("tmp_{RANDOM_ID}.fa")

def diamond(sequences, reference, json=False, verbose=True, out=None, sensitivity="very-sensitive"):
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
    
    input_file = create_input_file(sequences)
    print(f"input file name {input_file}")

    if out is None:
        command = f"diamond makedb --in {reference} -d reference && diamond blastp -q {input_file} -d reference -o tmp_out.tsv --{sensitivity}"
    else:
        # The double-quotation marks allow white spaces in the path, but this does not work for Windows
        command = f"diamond makedb --in {reference} -d reference && diamond blastp -q {input_file} -d reference -o {out}.tsv --{sensitivity}"
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
       
              
    df_diamond = tsv_to_df("diamond_out.tsv", ["query_accession", "target_accession", "Per. Ident" , "length", "mismatches", "gap_openings", "query_start", "query_end", "target_start", "target_end", "e-value", "bit_score"])
    remove_temp_files()
    return df_diamond
   

              