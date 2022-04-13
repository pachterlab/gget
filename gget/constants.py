import uuid

# Ensembl REST API server for gget seq and info
ENSEMBL_REST_API = "http://rest.ensembl.org/"

# UniProt REST API server for gget seq
UNIPROT_REST_API = "https://www.uniprot.org/uploadlists/"

# BLAST constants
BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
# Generate a random UUID
BLAST_CLIENT = "gget_client-" + str(uuid.uuid4())

# MUSCLE constants
MUSCLE_GITHUB_LINK = "https://github.com/rcedgar/muscle.git"