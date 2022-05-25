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
# Plotting packages
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import textwrap

# Constants
from .constants import POST_ENRICHR_URL, GET_ENRICHR_URL


def enrichr(genes, database, plot=False, json=False, save=False):
    """
    Perform an enrichment analysis on a list of genes using Enrichr (https://maayanlab.cloud/Enrichr/).

    Args:
    - genes       Genes to perform enrichment analysis on, passed as a list of strings,
                  e.g. ['PHF14', 'RBM3', 'MSL1', 'PHF21A'].
    - database    Database to use as reference for the enrichment analysis.
                  Supported shortcuts (and their default database):
                  'pathway' (KEGG_2021_Human)
                  'transcription' (ChEA_2016)
                  'ontology' (GO_Biological_Process_2021)
                  'diseases_drugs' (GWAS_Catalog_2019)
                  'celltypes' (PanglaoDB_Augmented_2021)
                  'kinase_interactions' (KEA_2015)
                  or any database listed under Gene-set Library at: https://maayanlab.cloud/Enrichr/#libraries
    - plot        True/False whether to provide a graphical overview of the first 15 results.
    - json        If True, returns results in json format instead of data frame. Default: False.
    - save        True/False whether to save the results in the local directory.

    Returns a data frame with the Enrichr results.
    """

    # Define database
    # All available libraries: https://maayanlab.cloud/Enrichr/#libraries
    db_message = f"""
    Please note that there might a more appropriate database for your application. 
    Go to https://maayanlab.cloud/Enrichr/#libraries for a full list of supported databases.
    """

    if database == "pathway":
        database = "KEGG_2021_Human"
        logging.info(
            f"Performing Enichr analysis using database {database}. " + db_message
        )
    elif database == "transcription":
        database = "ChEA_2016"
        logging.info(
            f"Performing Enichr analysis using database {database}. " + db_message
        )
    elif database == "ontology":
        database = "GO_Biological_Process_2021"
        logging.info(
            f"Performing Enichr analysis using database {database}. " + db_message
        )
    elif database == "diseases_drugs":
        database = "GWAS_Catalog_2019"
        logging.info(
            f"Performing Enichr analysis using database {database}. " + db_message
        )
    elif database == "celltypes":
        database = "PanglaoDB_Augmented_2021"
        logging.info(
            f"Performing Enichr analysis using database {database}. " + db_message
        )
    elif database == "kinase_interactions":
        database = "KEA_2015"
        logging.info(
            f"Performing Enichr analysis using database {database}. " + db_message
        )
    else:
        database = database

    logging.info(f"Performing Enichr analysis using database {database}.")

    # If single gene passed as string, convert to list
    if type(genes) == str:
        genes = [genes]

    # Remove any NaNs/Nones from the gene list
    genes_clean = []
    for gene in genes:
        if not gene == np.NaN and not gene is None:
            genes_clean.append(gene)
    # Join genes from list
    genes_clean = "\n".join(genes_clean)

    ## Submit gene list to Enrichr API
    args_dict = {
        "list": (None, genes_clean),
        "description": (None, "gget client gene list"),
    }

    r1 = requests.post(POST_ENRICHR_URL, files=args_dict)

    if not r1.ok:
        raise RuntimeError(
            f"Enrichr HTTP POST response status code: {r1.status_code}. "
            "Please double-check arguments and try again.\n"
        )

    # Get user ID
    post_results = r1.json()
    userListId = post_results["userListId"]

    ## Fetch results from Enrichr API
    # Build query with user ID
    query_string = f"?userListId={userListId}&backgroundType={database}"

    r2 = requests.get(GET_ENRICHR_URL + query_string)
    if not r2.ok:
        raise RuntimeError(
            f"Enrichr HTTP GET response status code: {r2.status_code}. "
            "Please double-check arguments and try again.\n"
        )

    enrichr_results = r2.json()

    # Return error if no results were found
    if len(enrichr_results) > 1:
        logging.error(
            f"No Enrichr results were found for your genes in database {database}. "
            "Please double-check the arguments and try again."
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
            f"No Enrichr results were found for genes {genes_clean} and database {database}."
        )

    ## Plot if plot=True
    if plot and len(df) != 0:
        fig, ax = plt.subplots(figsize=(10, 10))

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
        # plot = ax.scatter(c_values, c_values, c = c_values, cmap = cmap)
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
        # bar = ax.barh(labels, gene_counts, color=cmap(c_values), align="center")
        bar = ax.barh(labels, gene_counts, color=barcolor, align="center")
        ax.invert_yaxis()
        # Set x-limit to be gene count + 1
        ax.set_xlim(0, ax.get_xlim()[1] + 1)

        # # Add colorbar legend
        # cb = plt.colorbar(plot)
        # cb.set_label("$-log_{10}$(adjusted P value)", fontsize=fontsize)
        # cb.ax.tick_params(labelsize=fontsize)

        # Add adj. P value secondary x-axis
        ax2 = ax.twiny()
        ax2.scatter(-np.log10(adj_p_values), labels, color=p_val_color, s=20)
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
        ax.set_xlabel(
            f"Number of overlapping genes (query size: {len(genes)})",
            color=barcolor,
            fontsize=fontsize,
        )
        ax2.spines["bottom"].set_color(barcolor)
        ax.tick_params(axis="x", labelsize=fontsize, colors=barcolor)
        # Set bottom x axis to keep only integers since counts cannot be floats
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        # Change fontsize of y-tick labels
        ax.tick_params(axis="y", labelsize=fontsize)

        # Set title
        ax.set_title(f"Enrichr results from database {database}", fontsize=fontsize + 2)

        # Set axis margins
        ax.margins(y=0, x=0)

        plt.tight_layout()

        if save:
            fig.savefig("gget_enrichr_results.png", dpi=300, bbox_inches="tight")

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
