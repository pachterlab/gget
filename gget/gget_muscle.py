import sys
import time
import logging
# Add and format time stamp in logging messages
logging.basicConfig(format="%(asctime)s %(message)s", datefmt="%d %b %Y %H:%M:%S")
import os
import platform
import subprocess
import itertools

# Custom functions
from .compile import (
    compile_muscle,
    MUSCLE_PATH,
    PACKAGE_PATH
)
from .utils import (
    aa_colors,
    n_colors
)

# Path to precompiled muscle binary for this system
PRECOMPILED_MUSCLE_PATH = os.path.join(PACKAGE_PATH, f"bins/{platform.system()}/muscle")

def muscle(fasta, 
           super5=False, 
           out="default"
          ):
    """
    Align multiple nucleotide or amino acid sequences against each other (using the Muscle v5 algorithm).
    
    Args:
    - fasta     Path to fasta file containing the sequences to be aligned.
    - super5    True/False (default: False). 
                If True, align input using Super5 algorithm instead of PPP algorithm to decrease time and memory.
                Use for large inputs (a few hundred sequences).
    - out       Path to the 'aligned FASTA' (.afa) file the results will be saved in (default: "muscle_results.afa").
        
    Returns alignment results in an "aligned FASTA" (.afa) file.
    """
    # Muscle v5 documentation: https://drive5.com/muscle5

    # Get absolute path to input fasta file
    abs_fasta_path = os.path.abspath(fasta)
    
    # Get absolute path to output .afa file
    if out == "default":
        abs_out_path = os.path.join(os.getcwd(), "muscle_results.afa")
    else:
        abs_out_path = os.path.abspath(out)
    
    # Compile muscle if it is not already compiled
    if os.path.isfile(PRECOMPILED_MUSCLE_PATH) == False:
        # Compile muscle
        compile_muscle()
        muscle_path = MUSCLE_PATH
        
    else:
        logging.warning(
            "MUSCLE compiled. "
        )
        muscle_path = PRECOMPILED_MUSCLE_PATH
        
    # Assign read, write, and execute permission to muscle binary
    with subprocess.Popen(f"chmod 755 {muscle_path}", shell=True, stderr=subprocess.PIPE) as process_1:
        stderr_1 = process_1.stderr.read().decode("utf-8")
        # Log the standard error if it is not empty
        if stderr_1:
            logging.warning(stderr_1)
    # Exit system if the subprocess returned with an error
    if process_1.wait() != 0:
        return
    
    # Define muscle command
    if super5:
        command = f"{muscle_path} -super5 {abs_fasta_path} -output {abs_out_path}"
    else:
        command = f"{muscle_path} -align {abs_fasta_path} -output {abs_out_path}"
     
    # Record MUSCLE align start
    start_time = time.time()
    
    # Run muscle command and write command output
    with subprocess.Popen(command, shell=True, stderr=subprocess.PIPE) as process_2:
        stderr_2 = process_2.stderr.read().decode("utf-8")
        # Log the standard error if it is not empty
        if stderr_2:
            logging.warning("MUSCLE aligning... " + stderr_2)   
    # Exit system if the subprocess returned with an error
    if process_2.wait() != 0:
        return
    else:
        logging.warning(
            "MUSCLE alignment complete."
        )

        ## Print cleaned up muscle output
        # Get the titles and sequences from the generated .afa file
        titles = []
        seqs_master = []
        with open(abs_out_path) as aln_file:
            for i, line in enumerate(aln_file):
                # Recognize title lines by the '>' character
                if line[0] == ">":
                    # Record first listed identifier as title and remove the '>'
                    titles.append(line.split(" ")[0].split(">")[1])

                    # Append list containing seqs for the previous title to master list
                    if i != 0:
                        seqs_master.append(seqs)
                    # Empty the seqs list to append the sequences for this new title
                    seqs = []
                else:
                    seqs.append(line.strip())
        # Append seqs of last title to master seq list
        seqs_master.append(seqs)
        
        # Set of all possible nucleotides
        nucleotides = set("ATGC-")
        
        # zip_longest pads to the longest length 
        for seq_pair in list(itertools.zip_longest(*seqs_master)):
            print("\n")
            for idx, seq in enumerate(seq_pair):
                final_seq = []
                # If sequence is a nucleotide sequence,
                # assign nulceotide colors using the custom function n_colors
                if set(seq) <= nucleotides:
                    for letter in seq:
                        final_seq.append(n_colors(letter))

                # If sequence is not a nucleotide sequence, I assume it is an amino acid sequence
                # and assign amino acid colors using the custom function aa_colors
                else:
                    for letter in seq:
                        final_seq.append(aa_colors(letter))

                print(titles[idx], "".join(final_seq))
            
