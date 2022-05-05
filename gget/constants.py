import uuid

# Ensembl REST API server for gget seq and info
ENSEMBL_REST_API = "http://rest.ensembl.org/"

# UniProt REST API server for gget seq
UNIPROT_REST_API = "https://www.uniprot.org/uploadlists/"

# BLAST API endpoints
BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
# Generate a random UUID
BLAST_CLIENT = "gget_client-" + str(uuid.uuid4())

# MUSCLE Github repo
MUSCLE_GITHUB_LINK = "https://github.com/rcedgar/muscle.git"

# Enrichr API endpoints
POST_ENRICHR_URL = "https://maayanlab.cloud/Enrichr/addList"
GET_ENRICHR_URL = "https://maayanlab.cloud/Enrichr/enrich"

# Gene correlation API endpoints
GENECORR_URL = "https://maayanlab.cloud/matrixapi/coltop"