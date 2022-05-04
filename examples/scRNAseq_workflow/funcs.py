## Custom plotting functions for the scRNAseq workflow

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from scipy import stats
import scanpy as sc
from IPython.display import display, HTML
sc.set_figure_params(figsize=(4, 4), frameon=False)
sc.settings.n_jobs=2

# Plotting defaults
color = "#e66101"
filter_color = "#f7dbc6"
fontsize = 18

def yex(ax, col="k", alpha=0.75):
    """
    Function to add linear graph to plot.
    Call after defining x and y scales
    """
    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    ]
    # Plot both limits against each other
    ax.plot(lims, lims, col, alpha=alpha, zorder=4)
    ax.set(**{
        "aspect": "equal",
        "xlim": lims,
        "ylim": lims
    })
    return ax

def nd(arr):
    """
    Function to transform numpy matrix to nd array.
    """
    return np.asarray(arr).reshape(-1)

titles = ["Knee Plot", "Library Saturation Plot"]
def knee_library_sat_plot(adata, expected_num_cells):
    fig, axs = plt.subplots(figsize=(30, 14), ncols=1, nrows=2, sharex='col')

    ## Plot knee plots
    ax = axs[0]

    knee = np.sort((np.array(adata.X.sum(axis=1))).flatten())[::-1]
    # Create masks for cells that will be kept/filtered out based on minimum UMI count
    min_umi = knee[expected_num_cells]
    keep_mask = np.ma.masked_where(knee > min_umi, knee).mask
    filter_mask = np.invert(keep_mask)

    # Plot all cells in filter out color
    ax.plot(knee, range(len(knee)), linewidth=5, color=filter_color)
    # Overlay with cells to keep in dark color
    ax.plot(knee[keep_mask], range(len(knee[keep_mask])), linewidth=5, color=color)
    ax.set_xscale('log')
    ax.set_yscale('log')

    # Add lines to show filtering cut offs
    ax.axvline(x=knee[expected_num_cells], linewidth=1.5, color="k", ls="--")
    ax.axhline(y=expected_num_cells, linewidth=1.5, color="k", ls="--")

    # Define axis labels and title
    #     ax.set_xlabel("UMI Counts", fontsize=fontsize)
    ax.set_ylabel("Set of Barcodes", fontsize=fontsize)
    ax.set_title(titles[0], fontsize=fontsize+2)

    # Change fontsize of tick labels
    ax.tick_params(axis="both", labelsize=fontsize)

    # Set x and y lim to be the same as for the library sat plot
    ax.set_xlim(1, 10**6)
    ax.set_ylim(1, 10**6)

    # Add grid and set below graph
    ax.grid(True, which="both", color="lightgray")
    ax.set_axisbelow(True)

    # Add invisible linear graph to match square library sat plot
    yex(ax, alpha=0)

    ## Plot library saturation plot
    ax = axs[1]

    # Number of UMI counts per cell
    x = np.asarray(adata.X.sum(axis=1))[:,0]
    # Number of different genes detected per cell
    y = np.asarray(np.sum(adata.X>0, axis=1))[:,0]

    # Create masks for cells that will be kept/filtered out based on minimum UMI count from knee plot
    knee = np.sort((np.array(adata.X.sum(axis=1))).flatten())[::-1]
    min_umi = knee[expected_num_cells]
    keep_mask = np.ma.masked_where(x > min_umi, x).mask
    filter_mask = np.invert(keep_mask)

    ax.scatter(x[keep_mask], y[keep_mask], alpha=0.5, color=color, zorder=2)
    ax.scatter(x[filter_mask], y[filter_mask], alpha=0.5, color=filter_color, zorder=2)

    # Add lines to show filtering cut off
    ax.axvline(x=min_umi, linewidth=1.5, color="k", ls="--", zorder=3)

    ax.set_xscale('log')
    ax.set_yscale('log')

    # Define axis labels and title
    ax.set_xlabel("UMI Counts", fontsize=fontsize)
    ax.set_ylabel("Genes Detected", fontsize=fontsize)
    ax.set_title(titles[1], fontsize=fontsize+2)

    # Change fontsize of tick labels
    ax.tick_params(axis='both', labelsize=fontsize)

    # Set x and y lim
    ax.set_xlim(1, 10**6)
    ax.set_ylim(1, 10**6)

    # Add grid and set below graph
    ax.grid(True, which="both", color="lightgray")
    ax.set_axisbelow(True)

    # Add linear graph
    yex(ax)

    fig.show()

def library_sat_gene_fraction(adata, genes, expected_num_cells, gene_type="QC"):
    """
    Plot a library saturation with the fraction of 
    gene counts of passed genes indicated in a colormap.
    """
    fig, ax = plt.subplots(figsize=(10,7))

    # Calculate the sum counts of all mito genes per cell
    total_exp = []
    for idx, gene in enumerate(genes):
        gene_exp = nd(adata[:, adata.var.index.str.contains(gene)].X.todense())
        if idx == 0:
            total_exp = gene_exp
        else:
            # Add the expression of this mito gene to the list
            total_exp = total_exp + gene_exp

#     print(f"Maximum total UMI count for all apoptosis/stress genes in one cell: {total_exp.max()}.")

    # Calculate the fraction of the total sum per cell
    total_mito_fraction = total_exp / np.sum(total_exp)

    x = nd(adata.X.sum(1))
    y = nd((adata.X>0).sum(1))
    c = total_mito_fraction
    idx = c.argsort()
    x2, y2, c2 = x[idx], y[idx], c[idx]
    scatter = ax.scatter(x2, y2, c=c2, cmap="Reds")
    # Set range of colorbar
    scatter.set_clim(0, 0.001) 
    
    # Add line to show filtering cut off
    knee = np.sort((np.array(adata.X.sum(axis=1))).flatten())[::-1]
    min_umi = knee[expected_num_cells]
    ax.axvline(x=min_umi, linewidth=1.5, color="k", ls="--", zorder=3)

    # Add a colorbar legend to the last graph
    cb = fig.colorbar(scatter, ax=ax)
    cb.set_label(label=f"Fraction of {gene_type} gene counts", fontsize=fontsize)
    # Set ticks below colorbar
    cb.ax.xaxis.set_ticks_position("bottom")

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel("UMI Counts", fontsize=fontsize)
    ax.set_ylabel("Genes Detected", fontsize=fontsize)
    ax.set_xlim(1, 10**6)
    ax.set_ylim(1, 10**6)
    # ax.set_title(titles[1], fontsize=fontsize+2)
    
    # Change fontsize of tick labels
    ax.tick_params(axis="both", labelsize=fontsize)
    
    yex(ax)

    # Add grid and set below graph
    # ax.grid(True, which="both", color="lightgray")
    ax.set_axisbelow(True)

    plt.tight_layout()
    
    fig.show()
    
def find_ids(adata, mito_ids):
    """
    Find Ensembl IDs in an AnnData object.
    """
    mito_gene_ids = []
    for i, gene in enumerate(mito_ids):
        gene_name_id = []

        gni = adata.var.iloc[np.where(adata.var.index.str.contains(gene))]

        if len(gni) > 0:
            gene_name_id = gni.index[0]
        else:
            gene_name_id = np.nan

        mito_gene_ids.append(gene_name_id)
        
    return mito_gene_ids

def sc_heatmap(adata, marker_genes):
    """
    Use scanpy to plot heatmap of normalized gene expression.
    """
    sc.pl.heatmap(
        adata,
        marker_genes,
        groupby="leiden",
        use_raw=False, 
        cmap="inferno", 
        standard_scale="var",
        swap_axes=True,
        figsize=(15, 7)
    )
    
def rank_marker_genes(adata, marker_genes, marker_gene_dict):
    """
    Use scanpy to rank marker genes.
    """
    # Create a copy of adata
    adata_test = adata.copy()

    # Find indices of all marker genes in adata
    ens_idx = np.isin(adata_test.var_names, marker_genes)

    # Remove all genes from adata_test except marker genes
    adata_test = adata_test[:,ens_idx].copy()

    # Rank marker genes
    sc.tl.rank_genes_groups(
        adata_test, 
        groupby="leiden", 
        method="wilcoxon", 
        corr_method="bonferroni", 
        use_raw=False
    )
    
    # Plot ranked genes per cluster
    sc.pl.rank_genes_groups(adata_test, n_genes=5, sharey=False)

def celltype_heatmap(adata, colors, figsize=(8, 8)):
    # Rank gene group based on celltype
    sc.tl.rank_genes_groups(
        adata,
        groupby="celltype",
        use_raw=False,
        method="wilcoxon",
        corr_method="bonferroni",
    )
    
    # Create dendrorgam
    sc.tl.dendrogram(
        adata,
        groupby="celltype",
        use_raw=False,
        cor_method="pearson"
    )
    
    # Run sc.pl.rank_genes_groups_heatmap once to create adata.uns["celltype_colors"] object
    sc.pl.rank_genes_groups_heatmap(
        adata,
        show_gene_labels=False,
        use_raw=False,
        show=False
    )
    plt.close()

    # Relabel celltype colors
    adata.uns["celltype_colors"] = colors
    
    # Plot heatmap with dendrogram
    sc.pl.rank_genes_groups_heatmap(
        adata,
        show_gene_labels=False,
        use_raw=False,
        cmap="inferno",
        standard_scale="var",
        figsize=figsize
    )

def volcano_df(adata, control_mask, experiment_mask):
    control = np.array(adata[control_mask].X.mean(axis=0))[0]
    experiment = np.array(adata[experiment_mask].X.mean(axis=0))[0]
    
    # Create df
    df_volcano = pd.DataFrame()

    ## Add columns with gene names and mean counts for each gene for each batch
    df_volcano["GeneNames"] = adata.var.index
    df_volcano["control"] = control
    df_volcano["experiment"] = experiment

    ## Compute log2 fold change
    # Since we already normalized and logged the values, we just subtract them from each other in order to get the log fold change.
    df_volcano["logFC"] = df_volcano["control"] - df_volcano["experiment"]
    df_volcano["logFC"] = df_volcano["logFC"].fillna(0)

    ## Compute p-value
    df_volcano["p-value"] = 0
    # Dense sparse matrices
    matrix_1 = adata[control_mask].X.todense()
    matrix_2 = adata[experiment_mask].X.todense()
    # Compute and save p-value
    _, df_volcano["p-value"] = stats.ttest_ind(matrix_1, matrix_2, equal_var=False)

    return df_volcano

def volcano_plot(df_volcano, min_fold_change=2, alpha=0.05, figsize=(7, 7)):
    fig, ax = plt.subplots(figsize=figsize)
    
    # Define dotsize and translucence
    s = 20
    a = 0.5
    
    # Get fold change and p-value cutoffs
    xline = np.log(min_fold_change)
    yline = -np.log10(alpha)

    # Plot scatter plot
    x = df_volcano["logFC"].values
    y = -np.log10(df_volcano["p-value"].values.astype(float))
    labels = df_volcano["GeneNames"]
    
    ax.scatter(x, y, color="grey", s=s, alpha=a)

    mask1 = np.logical_and(x > xline, y > yline)
    ax.scatter(x[mask1], y[mask1], color="r", s=s, alpha=a)

    mask2 = np.logical_and(x < -xline, y > yline)
    ax.scatter(x[mask2], y[mask2], color="r", s=s, alpha=a)

    mask = np.logical_or(mask1, mask2)

    # Add lines to show cutoffs
    ax.axvline(x=-xline, color="grey", linestyle="--")
    ax.axhline(y=yline, color="grey", linestyle="--")
    ax.axvline(x=xline, color="grey", linestyle="--")

#     ax.text(ax.get_xlim()[1] - 1.1, yline + 1, f"p-value = {alpha}", fontsize=fontsize - 2)
#     ax.text(
#         -xline - 0.14,
#         ax.get_ylim()[1] - 25,
#         f"Fold Change = {min_fold_change}",
#         rotation="vertical",
#         fontsize=fontsize - 2,
#     )

    # Define axis labels
    ax.set_xlabel("$log_2$ Fold Change", fontsize=fontsize)
    ax.set_ylabel("$-log_{10}$ p-value", fontsize=fontsize)
    # Change fontsize of tick labels
    ax.tick_params(axis="both", labelsize=fontsize)

    ax.set_axisbelow(True)

    fig.show()
    
def pretty_print(df, cols):
    """
    Function to wrap columns cols of 
    a data frame df for easier reading.
    """
    for col in cols:
        df.loc[:, col] = df[col].str.wrap(30)
    
    return display(HTML(df.to_html().replace("\\n","<br>")))

def chr_locations(df_blat):
    fig, ax = plt.subplots(figsize=(10,5))

    fontsize = 13
    color = "#154A78"
    edgecolor = color

    x = df_blat["chromosome"].value_counts().index
    y = df_blat["chromosome"].value_counts().values
    ax.bar(x, y, color=color, edgecolor=edgecolor)

    ax.set_xlabel("Chromosome", fontsize=fontsize)
    ax.set_ylabel("Differentially Expressed Gene Counts", fontsize=fontsize)
    # Change fontsize of tick labels
    ax.tick_params(axis="both", labelsize=fontsize-1)
    
    # Set y axis to keep only integers since counts cannot be floats
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    # Add grid and set below graph
    ax.grid(True, which="both", color="lightgray")
    ax.set_axisbelow(True)

    fig.show()