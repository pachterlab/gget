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
from .gget_pdb import pdb
from .gget_gpt import gpt
from .gget_cellxgene import cellxgene
from .gget_elm import elm
from .gget_diamond import diamond
from .gget_cosmic import cosmic

import logging
logging.basicConfig(
    format="%(asctime)s %(levelname)s %(message)s",
    level=logging.INFO,
    datefmt="%c",
)
# Mute numexpr threads info
logging.getLogger("numexpr").setLevel(logging.WARNING)

__version__ = "0.28.4"
__author__ = "Laura Luebbert"
__email__ = "lauralubbert@gmail.com"
