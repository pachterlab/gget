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

import logging
# Add and format time stamp in logging messages
logger = logging.getLogger("gget")
formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
hdlr = logging.FileHandler("gget.log")
hdlr.setFormatter(formatter)
hdlr.setLevel(logging.INFO)
logger.addHandler(hdlr)
# Mute numexpr threads info
logging.getLogger("numexpr").setLevel(logging.WARNING)

__version__ = "0.2.7"
__author__ = "Laura Luebbert"
__email__ = "lauraluebbert@caltech.edu"
