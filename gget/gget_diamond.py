import logging
import subprocess
import sys

def diamond(input, reference, json=False, verbose=True, out=None):
    """
    Perform protein sequence alignment using DIAMOND for multiple sequences

    Args:
     - input          Input sequences path and file name (include ,fa) in FASTA file format
     - reference      Reference file path and file name (include .fa) in FASTA file format
     - json           If True, returns results in json format instead of data frame. Default: False.
     - out            folder name to save two resulting csv files. Default: results (default: None).
     - verbose        True/False whether to print progress information. Default True.

    Returns DIAMOND output in tsv format 
    """
    
    if out is None:
        command = f"diamond makedb --in {reference} -d reference && diamond blastp -q {input} -d reference -o diamond_out.tsv --very-sensitive"
    else:
        # The double-quotation marks allow white spaces in the path, but this does not work for Windows
        command = f"diamond makedb --in {reference} -d reference && diamond blastp -q {input} -d reference -o {out}.tsv --very-sensitive"
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