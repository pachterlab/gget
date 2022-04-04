import os

# Constants
from .constants import (
    MUSCLE_GITHUB_LINK
)

def compile_muscle():
    f"""
    Compile MUSCLE from source. 
    Link to Github repo: {MUSCLE_GITHUB_LINK}
    """
    command1 = "wget " + MUSCLE_GITHUB_LINK
    os.system(command1)
    
    command2 = "cd src/"
    os.system(command2)
    
    command3 = "make"
    os.system(command3)
