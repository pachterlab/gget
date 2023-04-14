import cellxgene_census
import logging

# Add and format time stamp in logging messages
logging.basicConfig(
    format="%(asctime)s %(levelname)s %(message)s",
    level=logging.INFO,
    datefmt="%c",
)
# Mute numexpr threads info
logging.getLogger("numexpr").setLevel(logging.WARNING)


def cellxgene(
    tissue=None,
    cell_type=None,
    gene=None,
    ensembl_id=False,
    sex=None,
    development_stage=None,
    disease=None,
    species="homo_sapiens",
    adata=True,
    column_names=[
        "dataset_id",
        "assay",
        "suspension_type",
        "sex",
        "tissue_general",
        "tissue",
        "cell_type",
    ],
):
    """
    Query data from CZ CELLxGENE Discover (https://cellxgene.cziscience.com/) using the
    CZ CELLxGENE Discover Census (https://github.com/chanzuckerberg/cellxgene-census).

    NOTE: Querying large datasets requires a large amount of RAM. Use the 'tissue', 'cell_type' and 'gene' 
    arguments to define the (sub)dataset of interest.
    The CZ CELLxGENE Discover Census recommends >16 GB of memory and a >5 Mbps internet connection.

    Args:
        - tissue              Str or list of tissues, e.g. ['lung', 'blood']. Default: None.
                              See https://cellxgene.cziscience.com/gene-expression for available tissues. Default: None.
        - cell_type           Str or list of celltypes, e.g. ['mucus secreting cell', 'neuroendocrine cell']. Default: None.
                              See https://cellxgene.cziscience.com/gene-expression and select a tissue for available celltypes.
        - gene                Str or list of gene names or Ensembl IDs, e.g. ['ACE2', 'SLC5A1'] or ['ENSG00000130234', 'ENSG00000100170']. Default: None.
                              NOTE: Set ensembl_id=True when providing Ensembl IDs instead of gene names.
                              See https://cellxgene.cziscience.com/gene-expression for available genes.
        - ensembl_id          True/False whether provided genes are Ensembl IDs or gene names. Default: False.
        - sex                 Str or list of sexes, e.g. 'female'. Default: None.
        - development_stage   Str or list of development stage(s). Default: None.
        - disease             Str or list of disease(s). Default: None.
        - species             Choice of 'homo_sapiens' or 'mus_musculus'. Default: 'homo_sapiens'.
        - column_names        List of columns to return (stored in .obs when adata=True).
                              Default: ["dataset_id", "assay", "suspension_type", "sex", "tissue_general", "tissue", "cell_type"]
                              For more options see: https://api.cellxgene.cziscience.com/curation/ui/#/ -> dataset
        - adata               True/False. Returns AnnData object. Default: True.
                              If False, only returns dataset metadata (corresponds to adata.obs).

    Returns AnnData object (when adata=True) or dataframe (when adata=False).
    """
    # Clean up arguments
    if isinstance(cell_type, str):
        cell_type = [cell_type]
    if isinstance(tissue, str):
        tissue = [tissue]
    if isinstance(development_stage, str):
        development_stage = [development_stage]
    if isinstance(disease, str):
        disease = [disease]
    if isinstance(gene, str):
        gene = [gene]
    if isinstance(sex, str):
        sex = [sex]

    # Define value filter
    args = [cell_type, tissue, development_stage, disease, sex]
    arg_names = ["cell_type", "tissue", "development_stage", "disease", "sex"]
    for i, (arg_name, arg) in enumerate(zip(arg_names, args)):
        if arg:
            if i == 0:
                obs_value_filter = f"{arg_name} in {str(arg)}"
            else:
                obs_value_filter = obs_value_filter + f" and {arg_name} in {str(arg)}"

    # Fetch AnnData object
    if adata:
        logging.info(
            "Fetching AnnData object from CZ CELLxGENE Discover. This might take a few minutes..."
        )
        with cellxgene_census.open_soma() as census:
            adata = cellxgene_census.get_anndata(
                census=census,
                organism=species,
                var_value_filter=f"{'feature_id' if ensembl_id else 'feature_name'} in {gene}",
                obs_value_filter=obs_value_filter,
                column_names={"obs": column_names},
            )

            return adata

    # Fetch metadata
    else:
        logging.info("Fetching metadata from CZ CELLxGENE Discover...")
        with cellxgene_census.open_soma() as census:

            # Reads SOMADataFrame as a slice
            cell_metadata = census["census_data"][species].obs.read(
                value_filter=obs_value_filter, column_names=column_names
            )

            # Concatenates results to pyarrow.Table
            cell_metadata = cell_metadata.concat()

            # Converts to pandas.DataFrame
            cell_metadata = cell_metadata.to_pandas()

            return cell_metadata
