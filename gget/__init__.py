from .gget_ref import ref
from .gget_search import search
from .gget_info import info
from .gget_seq import seq
from .gget_muscle import muscle
from .gget_blast import blast
from .gget_blat import blat
from .gget_enrichr import enrichr
from .gget_archs4 import archs4
from .gget_alphafold import alphafold
from .gget_setup import setup

import logging
logging.basicConfig(
    format="%(asctime)s %(levelname)s %(message)s",
    level=logging.INFO,
    datefmt="%c",
)
# Mute numexpr threads info
logging.getLogger("numexpr").setLevel(logging.WARNING)

__version__ = "0.3.10"
__author__ = "Laura Luebbert"
__email__ = "lauraluebbert@caltech.edu"
