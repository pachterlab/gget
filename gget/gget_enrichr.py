import requests
import pandas as pd
import json as json_package
import numpy as np
import logging

# Add and format time stamp in logging messages
logging.basicConfig(
    format="%(asctime)s %(levelname)s %(message)s",
    level=logging.INFO,
    datefmt="%c",
)
# Mute numexpr threads info
logging.getLogger("numexpr").setLevel(logging.WARNING)

# Plotting packages
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import textwrap

# Custom functions
from gget.gget_info import info

# Constants
from .constants import POST_ENRICHR_URL, GET_ENRICHR_URL, POST_BACKGROUND_ID_ENRICHR_URL, GET_BACKGROUND_ENRICHR_URL

def enrichr(
    genes,
    database,
    background=None,
    ensembl=False,
    plot=False,
    figsize=(10, 10),
    ax=None,
    json=False,
    save=False,
    verbose=True,
):
    """
    Perform an enrichment analysis on a list of genes using Enrichr (https://maayanlab.cloud/Enrichr/).

    Args:
    - genes       List of Entrez gene symbols to perform enrichment analysis on, passed as a list of strings, e.g. ['PHF14', 'RBM3', 'MSL1', 'PHF21A'].
                  Set 'ensembl = True' to input a list of Ensembl gene IDs, e.g. ['ENSG00000106443', 'ENSG00000102317', 'ENSG00000188895'].
    - background  List of background genes
    - database    Database to use as reference for the enrichment analysis.
                  Supported shortcuts (and their default database):
                  'pathway' (KEGG_2021_Human)
                  'transcription' (ChEA_2016)
                  'ontology' (GO_Biological_Process_2021)
                  'diseases_drugs' (GWAS_Catalog_2019)
                  'celltypes' (PanglaoDB_Augmented_2021)
                  'kinase_interactions' (KEA_2015)
                  or any database listed under Gene-set Library at: https://maayanlab.cloud/Enrichr/#libraries  
    - ensembl     Define as 'True' if 'genes' is a list of Ensembl gene IDs. (Default: False)
    - plot        True/False whether to provide a graphical overview of the first 15 results. (Default: False)
    - figsize     (width, height) of plot in inches. (Default: (10,10))
    - ax          Pass a matplotlib axes object for further customization of the plot. (Default: None)
    - json        If True, returns results in json format instead of data frame. (Default: False)
    - save        True/False whether to save the results in the local directory. (Default: False)
    - verbose     True/False whether to print progress information. Default True.

    Returns a data frame with the Enrichr results.
    """

    # Define database
    # All available libraries: https://maayanlab.cloud/Enrichr/#libraries
    db_message = f"""
    Please note that there might a more appropriate database for your application. 
    Go to https://maayanlab.cloud/Enrichr/#libraries for a full list of supported databases.
    """
    if verbose:
        logging.info(
            f"Performing Enichr analysis using database {database}. " + db_message
        )

    if database == "pathway":
        database = "KEGG_2021_Human"
      
    elif database == "transcription":
        database = "ChEA_2016"
    
    elif database == "ontology":
        database = "GO_Biological_Process_2021"
    
    elif database == "diseases_drugs":
        database = "GWAS_Catalog_2019"
  
    elif database == "celltypes":
        database = "PanglaoDB_Augmented_2021"
      
    elif database == "kinase_interactions":
        database = "KEA_2015"

    else:
        database = database
   
    # If single gene passed as string, convert to list
    if type(genes) == str:
        genes = [genes]

    ## Transform Ensembl IDs to gene symbols
    if ensembl:
        if verbose:
            logging.info("Getting gene symbols from Ensembl IDs.")

        genes_v2 = []

        for gene_id in genes:
            # Remove version number if passed
            gene_id = gene_id.split(".")[0]

            info_df = info(gene_id, pdb=False, ncbi=False, uniprot=False, verbose=False)

            # Check if Ensembl ID was found
            if isinstance(info_df, type(None)):
                logging.warning(
                    f"ID '{gene_id}' not found. Please double-check spelling/arguments."
                )
                continue

            gene_symbol = info_df.loc[gene_id]["ensembl_gene_name"]

            # If more than one gene symbol was returned, use first entry
            if isinstance(gene_symbol, list):
                genes_v2.append(str(gene_symbol[0]))
            else:
                genes_v2.append(str(gene_symbol))

        if verbose:
            logging.info(
                f"Performing Enichr analysis on the following gene symbols: {', '.join(genes_v2)}"
            )

    else:
        genes_v2 = genes

    # Remove any NaNs/Nones from the gene list
    genes_clean = []
    for gene in genes_v2:
        if not gene == np.NaN and not gene is None and not isinstance(gene, float):
            genes_clean.append(gene)

    if len(genes_clean) == 0 and ensembl:
        logging.error("No gene symbols found for given Ensembl IDs.")
        return

    # Join genes from list
    genes_clean_final = "\n".join(genes_clean)

    # Submit gene list to Enrichr API
    args_dict = {
        "list": (None, genes_clean_final),
        "description": (None, "gget client gene list"),
    }

    r1 = requests.post(POST_ENRICHR_URL, files=args_dict)

    
    # If single gene passed as string, convert to list
    if type(background) == str:
        background = [background]

    # Join background genes from list
    background_final = "\n".join(background)

    # Submit background list to Enrichr API to get background id
    if background is not None:
        # Submit gene list to Enrichr API
        args_dict_background = {
            "background": (None, background_final),
        }

        request_background_id = requests.post(POST_BACKGROUND_ID_ENRICHR_URL, files=args_dict_background)

        # Get background ID
        post_results_background= request_background_id.json()
        background_list_id = post_results_background["backgroundid"]

    if not r1.ok:
        raise RuntimeError(
            f"Enrichr HTTP POST response status code: {r1.status_code}. "
            "Please double-check arguments and try again.\n"
        )

    # Get user ID
    post_results = r1.json()
    userListId = post_results["userListId"]

    if background is None:
        query_string = f"?userListId={userListId}&backgroundType={database}"
        r2 = requests.get(GET_ENRICHR_URL + query_string)
    else:
        query_string = f"?userListId={userListId}&backgroundid={background_list_id}&backgroundType={database}"
        r2 = requests.get(GET_BACKGROUND_ENRICHR_URL + query_string)

    if not r2.ok:
        raise RuntimeError(
            f"Enrichr HTTP GET response status code: {r2.status_code}. "
            "Please double-check arguments and try again.\n"
        )


    enrichr_results = r2.json()


    # Return error if no results were found
    if len(enrichr_results) > 1:
        logging.error(
            f"No Enrichr results were found for genes {genes_clean} and database {database}. \n"
            "If the genes are Ensembl IDs, please set argument 'ensembl=True'. (For command-line, add flag [-e][--ensembl].)"
        )
        return

    ## Build data frame (standard return)
    # Define column names
    columns = [
        "rank",
        "path_name",
        "p_val",
        "z_score",
        "combined_score",
        "overlapping_genes",
        "adj_p_val",
        "Old p-value",
        "Old adjusted p-value",
    ]
    try:
        # Create data frame from Enrichr results
        df = pd.DataFrame(enrichr_results[database], columns=columns)

    except KeyError:
        logging.error(
            f"Database {database} not found. Go to https://maayanlab.cloud/Enrichr/#libraries "
            "for a full list of supported databases."
        )
        return

    # Drop last two columns ("Old p-value", "Old adjusted p-value")
    df = df.iloc[:, :-2]

    # Add database column
    df["database"] = database

    if len(df) == 0:
        logging.warning(
            f"No Enrichr results were found for genes {genes_clean} and database {database}. \n"
            "If the genes are Ensembl IDs, please set argument 'ensembl=True' (for terminal, add flag: [--ensembl])."
        )
