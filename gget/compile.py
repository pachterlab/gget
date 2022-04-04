import os
import logging
# Add and format time stamp in logging messages
logging.basicConfig(format="%(asctime)s %(message)s", datefmt="%d %b %Y %H:%M:%S")

# Constants
from .constants import (
    MUSCLE_GITHUB_LINK
)

def compile_muscle():
    f"""
    Compile MUSCLE from source. 
    Link to Github repo: {MUSCLE_GITHUB_LINK}
    """
    
    logging.warning(
        "Compiling `muscle` binary from source. "
    )
    
    # Clone MUSCLE repo
    command1 = "git clone " + MUSCLE_GITHUB_LINK
    os.system(command1)
    
    # Move to src folder
    os.chdir("muscle/src/")
    
    # Run make command
    command3 = "make"
    os.system(command3)
    
    logging.warning(
        "'Muscle' binary compiled. "
    )
