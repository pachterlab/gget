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

from .constants import (
    POST_ENRICHR_URL,
    GET_ENRICHR_URL,
    POST_BACKGROUND_ID_ENRICHR_URL,
    GET_BACKGROUND_ENRICHR_URL,
)
from .compile import PACKAGE_PATH
from .gget_info import info

def ensembl_to_gene_names(ensembl_ids):
    """
    Function to fetch gene names from a list of Ensembl IDs using gget info.
    """
    genes_v2 = []

    for gene_id in ensembl_ids:
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

    return genes_v2

def clean_genes_list(genes_list):
    # Remove any NaNs/Nones from the gene list
    genes_clean = []
    for gene in genes_list:
        if not isinstance(gene, float) and not gene is None and not gene == "nan":
            genes_clean.append(gene)
    return genes_clean


def enrichr(
    genes,
    database,
    background_list=None,
    background=False,
    ensembl=False,
    ensembl_bkg=False,
    plot=False,
    figsize=(10, 10),
    ax=None,
    kegg_out=None,
    kegg_rank=1,
    json=False,
    save=False,
    verbose=True,
):
    """
    Perform an enrichment analysis on a list of genes using Enrichr (https://maayanlab.cloud/Enrichr/).

    Args:
    - genes             List of Entrez gene symbols to perform enrichment analysis on, passed as a list of strings, e.g. ['PHF14', 'RBM3', 'MSL1', 'PHF21A'].
                        Set 'ensembl = True' to input a list of Ensembl gene IDs, e.g. ['ENSG00000106443', 'ENSG00000102317', 'ENSG00000188895'].
    - database          Database to use as reference for the enrichment analysis.
                        Supported shortcuts (and their default database):
                        'pathway' (KEGG_2021_Human)
                        'transcription' (ChEA_2016)
                        'ontology' (GO_Biological_Process_2021)
                        'diseases_drugs' (GWAS_Catalog_2019)
                        'celltypes' (PanglaoDB_Augmented_2021)
                        'kinase_interactions' (KEA_2015)
                        or any database listed under Gene-set Library at: https://maayanlab.cloud/Enrichr/#libraries
    - background_list   List of gene names/Ensembl IDs to be used as background genes. (Default: None)
    - background        If True, use set of > 20,000 default background genes listed here: https://github.com/pachterlab/gget/blob/main/gget/constants/enrichr_bkg_genes.txt.
                        (Default: False)
    - ensembl           Define as 'True' if 'genes' is a list of Ensembl gene IDs. (Default: False)
    - ensembl_bkg       Define as 'True' if 'background_list' is a list of Ensembl gene IDs. (Default: False)
    - plot              True/False whether to provide a graphical overview of the first 15 results. (Default: False)
    - figsize           (width, height) of plot in inches. (Default: (10,10))
    - ax                Pass a matplotlib axes object for further customization of the plot. (Default: None)
    - kegg_out          Path to file to save the highlighted KEGG pathway image, e.g. path/to/folder/kegg_pathway.png. (Default: None)
    - kegg_rank         Candidate pathway rank to be plotted in KEGG pathway image. (Default: 1)
    - json              If True, returns results in json format instead of data frame. (Default: False)
    - save              True/False whether to save the results in the local directory. (Default: False)
    - verbose           True/False whether to print progress information. (Default: True)

    Returns a data frame with the Enrichr results.
    """

    # Define database
    # All available libraries: https://maayanlab.cloud/Enrichr/#libraries
    db_message = f"""
    Please note that there might be a more appropriate database for your application. 
    Go to https://maayanlab.cloud/Enrichr/#libraries for a full list of supported databases.
    """
    if not (type(background) == bool):
        raise ValueError(
            f"Argument`background` must be a boolean True/False. If you are adding a background list, use the argument `background_list` instead."
        )

    # Handle database shortcuts
    if database == "pathway":
        database = "KEGG_2021_Human"
        if verbose:
            logging.info(
                f"Performing Enichr analysis using database {database}. " + db_message
            )

    elif database == "transcription":
        database = "ChEA_2016"
        if verbose:
            logging.info(
                f"Performing Enichr analysis using database {database}. " + db_message
            )

    elif database == "ontology":
        database = "GO_Biological_Process_2021"
        if verbose:
            logging.info(
                f"Performing Enichr analysis using database {database}. " + db_message
            )

    elif database == "diseases_drugs":
        database = "GWAS_Catalog_2019"
        if verbose:
            logging.info(
                f"Performing Enichr analysis using database {database}. " + db_message
            )

    elif database == "celltypes":
        database = "PanglaoDB_Augmented_2021"
        if verbose:
            logging.info(
                f"Performing Enichr analysis using database {database}. " + db_message
            )

    elif database == "kinase_interactions":
        database = "KEA_2015"
        if verbose:
            logging.info(
                f"Performing Enichr analysis using database {database}. " + db_message
            )

    else:
        database = database
        if verbose:
            logging.info(f"Performing Enichr analysis using database {database}.")

    # To generate a KEGG pathway image, confirm that the database is a KEGG database and pykegg is installed
    if kegg_out:
        if not database.startswith("KEGG"):
            logging.error(
                "Please specify a KEGG database when generating a KEGG pathway image."
            )
            return
        try:
            import pykegg
        except ImportError:
            logging.error(
                "Please install `pykegg` to generate a KEGG pathway image. Pykegg can be installed using pip: 'pip install pykegg'"
            )
            return

    # If single gene passed as string, convert to list
    if type(genes) == str:
        genes = [genes]

    ## Transform Ensembl IDs to gene symbols
    if ensembl:
        if verbose:
            logging.info("Getting gene symbols from Ensembl IDs.")

        genes_v2 = ensembl_to_gene_names(genes)

    else:
        genes_v2 = genes

    if len(genes_v2) == 0 and ensembl:
        logging.error("No gene symbols found for given Ensembl IDs.")
        return

    # Transform Ensembl IDs to gene symbols for background genes
    if background_list and ensembl_bkg:
        background_list = ensembl_to_gene_names(background_list)

    if not isinstance(background_list, type(None)):
        if len(background_list) == 0 and ensembl:
            logging.error("No background gene symbols found for given Ensembl IDs.")
            return

    genes_clean = clean_genes_list(genes_v2)

    if ensembl:
        if verbose:
            logging.info(
                f"Performing Enichr analysis on the following gene symbols: {', '.join(genes_clean)}"
            )

    # Join genes from list
    genes_clean_final = "\n".join(genes_clean)

    # Remove any NaNs/Nones from the background list
    if background_list:
        background_list = clean_genes_list(background_list)

    # Submit gene list to Enrichr API
    args_dict = {
        "list": (None, genes_clean_final),
        "description": (None, "gget client gene list"),
    }

    r1 = requests.post(POST_ENRICHR_URL, files=args_dict)

    if not r1.ok:
        raise RuntimeError(
            f"Enrichr HTTP POST gene list response status code: {r1.status_code}. "
            "Please double-check arguments and try again.\n"
        )

    # Get user ID
    post_results = r1.json()
    userListId = post_results["userListId"]

    # Get background genes list from user or from file of all genes
    background_final = None

    # If user gives a background list, use the user input instead of the default
    if background_list:
        if verbose:
            logging.info(
                f"Performing Enichr analysis using user-defined background gene list."
            )

        if background:
            logging.warning(
                "Since you provided a list of background genes, the 'background==True' argument to use the default background gene list is being ignored."
            )
        background_final = "\n".join(background_list)

    elif background:
        if verbose:
            logging.info(
                "Background genes set to > 20,000 default background genes listed here: https://github.com/pachterlab/gget/blob/main/gget/constants/enrichr_bkg_genes.txt."
            )
        with open(f"{PACKAGE_PATH}/constants/enrichr_bkg_genes.txt") as f:
            lines = f.read().splitlines()
        background_final = "\n".join(lines)

    # Submit background list to Enrichr API to get background id
    background_list_id = None
    if background_final:
        args_dict_background = {
            "background": (None, background_final),
        }

        request_background_id = requests.post(
            POST_BACKGROUND_ID_ENRICHR_URL, files=args_dict_background
        )

        if not request_background_id.ok:
            raise RuntimeError(
                f"""
                Enrichr HTTP POST background gene list response status code: {request_background_id.status_code}. \n
                Please double-check arguments and try again.\n
                """
            )

        # Get background ID
        post_results_background = request_background_id.json()
        background_list_id = post_results_background["backgroundid"]

    # Submit query to Enrich using gene list and background genes list
    if not background_final:
        query_string = f"?userListId={userListId}&backgroundType={database}"
        r2 = requests.get(GET_ENRICHR_URL + query_string)
    else:
        query_string = f"?userListId={userListId}&backgroundid={background_list_id}&backgroundType={database}"
        r2 = requests.post(GET_BACKGROUND_ENRICHR_URL + query_string)

    if not r2.ok:
        if background_final:
            raise RuntimeError(
                f"""
                Enrichr HTTP GET response status code: {r2.status_code} for genes {genes_clean}, background genes {background_list}, and database {database}\n
                This can be due to no results found by Enrichr.
                If the input genes are Ensembl IDs, please set argument 'ensembl=True'. (For command-line, add flag [-e][--ensembl].)\n
                If the background genes are Ensembl IDs, please set argument 'ensembl_bkg=True'. (For command-line, add flag [-e_b][--ensembl_bkg].\n
                """
            )
        else:
            raise RuntimeError(
                f"""
                Enrichr HTTP GET response status code: {r2.status_code} for genes {genes_clean}, and database {database}\n
                If the input genes are Ensembl IDs, please set argument 'ensembl=True'. (For command-line, add flag [-e][--ensembl].)\n
                """
            )

    enrichr_results = r2.json()

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
            f"""
            Database {database} not found. Go to https://maayanlab.cloud/Enrichr/#libraries 
            for a full list of supported databases.
            """
        )
        return

    # Drop last two columns ("Old p-value", "Old adjusted p-value")
    df = df.iloc[:, :-2]

    # Add database column
    df["database"] = database

    if len(df) == 0:
        logging.error(
            f"""
            No Enrichr results were found for genes {genes_clean} and database {database}. \n
            If the genes are Ensembl IDs, please set argument 'ensembl=True' (for terminal, add flag: [--ensembl]).
            """
        )

    ## Plot if plot=True
    if plot and len(df) != 0:
        if ax is None:
            fig, ax1 = plt.subplots(figsize=figsize)
        else:
            ax1 = ax

        fontsize = 12
        barcolor = "indigo"
        p_val_color = "darkorange"

        # Only plot first 15 results
        if len(df) > 15:
            overlapping_genes = df["overlapping_genes"].values[:15]
            path_names = df["path_name"].values[:15]
            adj_p_values = df["adj_p_val"].values[:15]
        else:
            overlapping_genes = df["overlapping_genes"].values
            path_names = df["path_name"].values
            adj_p_values = df["adj_p_val"].values

        # # Define bar colors by adj. p-value
        # cmap = plt.get_cmap("viridis")
        # c_values = -np.log10(adj_p_values)
        # # Plot scatter to use for colorbar legend
        # plot = ax1.scatter(c_values, c_values, c = c_values, cmap = cmap)
        # # Clear axis to remove unnecessary scatter
        # plt.cla()

        # Get gene counts
        gene_counts = []
        for gene_list in overlapping_genes:
            gene_counts.append(len(gene_list))

        # Wrap pathway labels
        labels = []
        for label in path_names:
            labels.append(
                textwrap.fill(
                    label,
                    width=40,
                    break_long_words=False,
                    max_lines=2,
                    placeholder="...",
                )
            )

        # Plot barplot
        # ax1.barh(np.arange(len(gene_counts)), gene_counts, color=cmap(c_values), align="center")
        ax1.barh(
            np.arange(len(gene_counts)), gene_counts, color=barcolor, align="center"
        )
        ax1.set_yticks(
            np.arange(len(gene_counts)), labels, linespacing=0.85, fontsize=fontsize
        )
        ax1.invert_yaxis()
        # Set x-limit to be gene count + 1
        ax1.set_xlim(0, ax1.get_xlim()[1] + 1)

        # # Add colorbar legend
        # cb = plt.colorbar(plot)
        # cb.set_label("$-log_{10}$(adjusted P value)", fontsize=fontsize)
        # cb.ax1.tick_params(labelsize=fontsize)

        # Add adj. P value secondary x-axis
        ax2 = ax1.twiny()
        ax2.scatter(
            -np.log10(adj_p_values),
            np.arange(len(gene_counts)),
            color=p_val_color,
            s=20,
        )
        # Change label and color of p-value axis
        ax2.set_xlabel(
            "$-log_{10}$(adjusted P value)", fontsize=fontsize, color=p_val_color
        )
        ax2.spines["top"].set_color(p_val_color)
        ax2.tick_params(axis="x", colors=p_val_color, labelsize=fontsize)

        # Add alpha=0.05 p-value cutoff
        ax2.axvline(-np.log10(0.05), color=p_val_color, ls="--", lw=2, alpha=0.5)
        ax2.text(
            -np.log10(0.05) + 0.02,
            -0.3,
            "p = 0.05",
            ha="left",
            va="top",
            rotation="vertical",
            color=p_val_color,
            fontsize=fontsize,
            alpha=0.5,
        )

        # Set label and color of count axis
        ax1.set_xlabel(
            f"Number of overlapping genes (query size: {len(genes_clean)})",
            color=barcolor,
            fontsize=fontsize,
        )
        ax2.spines["bottom"].set_color(barcolor)
        ax1.tick_params(axis="x", labelsize=fontsize, colors=barcolor)
        # Set bottom x axis to keep only integers since counts cannot be floats
        ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
        # Change fontsize of y-tick labels
        ax1.tick_params(axis="y", labelsize=fontsize)

        # Set title
        ax1.set_title(
            f"Enrichr results from database {database}", fontsize=fontsize + 2
        )

        # Set axis margins
        ax1.margins(y=0, x=0)

        # Remove grids
        ax1.grid(False)
        ax2.grid(False)

        plt.tight_layout()

        if save:
            fig.savefig(
                "gget_enrichr_results.png",
                dpi=300,
                bbox_inches="tight",
                transparent=True,
            )

    # Generate KEGG pathway image
    if kegg_out:
        candidate_rank = df[df["rank"] == kegg_rank].iloc[0, :]
        kegg_img = pykegg.visualize(
            candidate_rank["path_name"],
            candidate_rank["overlapping_genes"],
            db=database,
            output=kegg_out,
        )

    if json:
        results_dict = json_package.loads(df.to_json(orient="records"))
        if save:
            with open("gget_enrichr_results.json", "w", encoding="utf-8") as f:
                json_package.dump(results_dict, f, ensure_ascii=False, indent=4)

        # Return results in json format
        return results_dict

    else:
        if save:
            df.to_csv("gget_enrichr_results.csv", index=False)

        # Return data frame
        return df
