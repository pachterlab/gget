import os
import platform
import logging
# Add and format time stamp in logging messages
logging.basicConfig(format="%(asctime)s %(message)s", datefmt="%d %b %Y %H:%M:%S")

# Constants
from .constants import (
    MUSCLE_GITHUB_LINK
)

# Get absolute package path
PACKAGE_PATH = os.path.abspath(os.path.dirname(__file__))
# Path to muscle binary (only exists after 'compile_muscle' was executed)
MUSCLE_PATH = os.path.join(PACKAGE_PATH, "bins/compiled/muscle/src/Linux/muscle")

def compile_muscle():
    f"""
    Compiles MUSCLE from source. 
    Link to MUSCLE Github repo: {MUSCLE_GITHUB_LINK}
    Currently does not support Windows platforms.
    """
    
    logging.warning(
        "Compiling MUSCLE binary from source... "
    )
    
    # Change path to package path
    os.chdir(PACKAGE_PATH)
    # Create folders 'bins/compiled/' inside gget package
    os.makedirs("bins/compiled/", exist_ok=True)
    # Change path to PACKAGE_PATH/bins/compiled/
    os.chdir("bins/compiled/")
    
    # Clone MUSCLE repo into PACKAGE_PATH/bins/compiled/
    command1 = "git clone " + MUSCLE_GITHUB_LINK + " -q"
    os.system(command1)
    
    # Change path to PACKAGE_PATH/bins/compiled/muscle/src/
    os.chdir("muscle/src/")
    
    # Run make command
    command2 = "make -s"
    os.system(command2)
    
    logging.warning(
        "MUSCLE compiled. "
    )
