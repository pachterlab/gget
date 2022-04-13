import sys
import time
import logging
# Add and format time stamp in logging messages
logging.basicConfig(format="%(asctime)s %(message)s", datefmt="%d %b %Y %H:%M:%S")
import os
import subprocess
# Custom functions
from .compile import (
    compile_muscle,
    MUSCLE_PATH,
    PACKAGE_PATH
)

# Path to precompiled muscle binary
PRECOMPILED_MUSCLE_PATH = os.path.join(PACKAGE_PATH, "bins/linux/muscle")

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

    # Define muscle command
    if super5:
        command = f"{muscle_path} -super5 {abs_fasta_path} -output {abs_out_path}"
    else:
        command = f"{muscle_path} -align {abs_fasta_path} -output {abs_out_path}"
     
    # Run align command  
    logging.warning(
            "MUSCLE aligning... "
        )
    start_time = time.time()

    # Assign read, write, and execute permission to muscle binary
    os.system(f"chmod 755 {muscle_path}")
    
    # Run muscle command
    os.system(command)

    # process2 = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # # Return standard output
    # print(process2.stdout)
    # # Exit system if process returned error code
    # if process2.returncode != 0:
    #     sys.exit()

    # logging.warning(
    #     f"MUSCLE alignment complete. Alignment time: {round(time.time() - start_time, 2)} seconds."
    # )
