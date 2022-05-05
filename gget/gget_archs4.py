import requests
import pandas as pd
import logging

# Add and format time stamp in logging messages
logging.basicConfig(
    format="%(asctime)s %(levelname)s %(message)s",
    level=logging.INFO,
    datefmt="%c",
)

# Constants
from .constants import GENECORR_URL

def archs4(gene, which="correlation"):
    """
    Function to find the most similar genes based on co-expression
    and tissue expression using data from the ARCHS4 human and mouse 
    RNA-seq database (https://maayanlab.cloud/archs4/).

    Args:
    - gene      Short name of gene (gene symbol) to search (str).
    - which     'correlation' (default) or 'tissue'.

    If which='correlation': Returns a data frame with the gene symbols and 
    Pearson Correlation values of the 100 most similar genes.
    """

    ## Find most similar genes based on co-expression
    # Number of correlated genes to return
    gene_count = 101

     # Dictionary with arguments
    json_dict = {
        'id':gene,
        'count':gene_count
        }

    r = requests.post(url=GENECORR_URL, json=json_dict)

    if not r.ok:
        raise RuntimeError(
            f"Gene correlation API request returned with error code: {r.status_code}. "
            "Please double-check the arguments and try again.\n"
        )

    corr_data = r.json()

    # Check if the request returned an error (e.g. gene not found)
    if "error" in corr_data.keys():
        if corr_data["error"] == f'{gene} not in colids':
            logging.error(f"Search term '{gene}' did not return any results.")
            return
        else:
            logging.error(f"Gene correlation request for search term '{gene}' returned error: {corr_data['error']}")
            return

    else:
        corr_df = pd.DataFrame()
        corr_df["gene_symbol"] = corr_data["rowids"]
        corr_df["pearson_correlation"] = corr_data["values"]
        # Drop the first row (since that is the searched gene against itself)
        corr_df = corr_df.iloc[1: , :]

    return corr_df