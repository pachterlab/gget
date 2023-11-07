import logging

# Add and format time stamp in logging messages
logging.basicConfig(
    format="%(asctime)s %(levelname)s %(message)s",
    level=logging.INFO,
    datefmt="%c",
)
# Mute numexpr threads info
logging.getLogger("numexpr").setLevel(logging.WARNING)


def convert_to_list(lst):
    """
    Function to convert all non-list instances in a list to list.
    Returns list of lists.
    """
    temp = []
    for el in lst:
        if isinstance(el, str):
            temp.append([el])
        else:
            temp.append(el)
    return temp


def cellxgene(
    species="homo_sapiens",
    gene=None,
    ensembl=False,
    column_names=[
        "dataset_id",
        "assay",
        "suspension_type",
        "sex",
        "tissue_general",
        "tissue",
        "cell_type",
    ],
    meta_only=False,
    tissue=None,
    cell_type=None,
    development_stage=None,
    disease=None,
    sex=None,
    is_primary_data=True,
    dataset_id=None,
    tissue_general_ontology_term_id=None,
    tissue_general=None,
    assay_ontology_term_id=None,
    assay=None,
    cell_type_ontology_term_id=None,
    development_stage_ontology_term_id=None,
    disease_ontology_term_id=None,
    donor_id=None,
    self_reported_ethnicity_ontology_term_id=None,
    self_reported_ethnicity=None,
    sex_ontology_term_id=None,
    suspension_type=None,
    tissue_ontology_term_id=None,
    census_version="stable",
    verbose=True,
    out=None,
):
    """
    Query data from CZ CELLxGENE Discover (https://cellxgene.cziscience.com/) using the
    CZ CELLxGENE Discover Census (https://github.com/chanzuckerberg/cellxgene-census).

    NOTE: Querying large datasets requires a large amount of RAM. Use the cell metadata attributes
    to define the (sub)dataset of interest.
    The CZ CELLxGENE Discover Census recommends >16 GB of memory and a >5 Mbps internet connection.

    General args:
        - species        Choice of 'homo_sapiens' or 'mus_musculus'. Default: 'homo_sapiens'.
        - gene           Str or list of gene name(s) or Ensembl ID(s), e.g. ['ACE2', 'SLC5A1'] or ['ENSG00000130234', 'ENSG00000100170']. Default: None.
                         NOTE: Set ensembl=True when providing Ensembl ID(s) instead of gene name(s).
                         See https://cellxgene.cziscience.com/gene-expression for examples of available genes.
        - ensembl        True/False (default: False). Set to True when genes are provided as Ensembl IDs.
        - column_names   List of metadata columns to return (stored in AnnData.obs when meta_only=False).
                         Default: ["dataset_id", "assay", "suspension_type", "sex", "tissue_general", "tissue", "cell_type"]
                         For more options see: https://api.cellxgene.cziscience.com/curation/ui/#/ -> Schemas -> dataset
        - meta_only      True/False (default: False). If True, returns only metadata dataframe (corresponds to AnnData.obs).
        - census_version Str defining version of Census, e.g. "2023-05-15" or "latest" or "stable". Default: "stable".
        - verbose        True/False whether to print progress information. Default True.
        - out            If provided, saves the generated AnnData h5ad (or csv when meta_only=True) file with the specified path. Default: None.

    Cell metadata attributes:
        - tissue                          Str or list of tissue(s), e.g. ['lung', 'blood']. Default: None.
                                          See https://cellxgene.cziscience.com/gene-expression for examples of available tissues.
        - cell_type                       Str or list of celltype(s), e.g. ['mucus secreting cell', 'neuroendocrine cell']. Default: None.
                                          See https://cellxgene.cziscience.com/gene-expression and select a tissue to see examples of available celltypes.
        - development_stage               Str or list of development stage(s). Default: None.
        - disease                         Str or list of disease(s). Default: None.
        - sex                             Str or list of sex(es), e.g. 'female'. Default: None.
        - is_primary_data                 True/False (default: True). If True, returns only the canonical instance of the cellular observation.
                                          This is commonly set to False for meta-analyses reusing data or for secondary views of data.
        - dataset_id                      Str or list of CELLxGENE dataset ID(s). Default: None.
        - tissue_general_ontology_term_id Str or list of high-level tissue UBERON ID(s). Default: None.
                                          Also see: https://github.com/chanzuckerberg/single-cell-data-portal/blob/9b94ccb0a2e0a8f6182b213aa4852c491f6f6aff/backend/wmg/data/tissue_mapper.py
        - tissue_general                  Str or list of high-level tissue label(s). Default: None.
                                          Also see: https://github.com/chanzuckerberg/single-cell-data-portal/blob/9b94ccb0a2e0a8f6182b213aa4852c491f6f6aff/backend/wmg/data/tissue_mapper.py
        - tissue_ontology_term_id         Str or list of tissue ontology term ID(s) as defined in the CELLxGENE dataset schema. Default: None.
        - assay_ontology_term_id          Str or list of assay ontology term ID(s) as defined in the CELLxGENE dataset schema. Default: None.
        - assay                           Str or list of assay(s) as defined in the CELLxGENE dataset schema. Default: None.
        - cell_type_ontology_term_id      Str or list of celltype ontology term ID(s) as defined in the CELLxGENE dataset schema. Default: None.
        - development_stage_ontology_term_id        Str or list of development stage ontology term ID(s) as defined in the CELLxGENE dataset schema. Default: None.
        - disease_ontology_term_id        Str or list of disease ontology term ID(s) as defined in the CELLxGENE dataset schema. Default: None.
        - donor_id                        Str or list of donor ID(s) as defined in the CELLxGENE dataset schema. Default: None.
        - self_reported_ethnicity_ontology_term_id  Str or list of self reported ethnicity ontology ID(s) as defined in the CELLxGENE dataset schema. Default: None.
        - self_reported_ethnicity         Str or list of self reported ethnicity as defined in the CELLxGENE dataset schema. Default: None.
        - sex_ontology_term_id            Str or list of sex ontology ID(s) as defined in the CELLxGENE dataset schema. Default: None.
        - suspension_type                 Str or list of suspension type(s) as defined in the CELLxGENE dataset schema. Default: None.

    Returns AnnData object (when meta_only=False) or dataframe (when meta_only=True).
    """
    # Check if cellxgene_census is installed
    try:
        import cellxgene_census
    except ImportError:
        logging.error(
            """
            Some third-party dependencies are missing. Please run the following command: 
            >>> gget.setup('cellxgene') or $ gget setup cellxgene

            Alternative: Install the cellxgene-census package using pip (https://pypi.org/project/cellxgene-census).
            """
        )
        return

    # List of metadata arguments
    args = [
        dataset_id,
        tissue_general_ontology_term_id,
        tissue_general,
        assay_ontology_term_id,
        assay,
        cell_type_ontology_term_id,
        cell_type,
        development_stage_ontology_term_id,
        development_stage,
        disease_ontology_term_id,
        disease,
        donor_id,
        self_reported_ethnicity_ontology_term_id,
        self_reported_ethnicity,
        sex_ontology_term_id,
        sex,
        suspension_type,
        tissue_ontology_term_id,
        tissue,
    ]

    if all(el is None for el in args):
        logging.warning(
            """
            You are attempting to query the entire Census dataset which requires a large amount of RAM (100's of GBs) and high network bandwidth. 
            Use the cell metadata arguments (e.g. 'tissue', 'cell_type', 'disease', etc...) to define the (sub)dataset of interest.
            """
        )

    # Convert args to string to get argument names
    arg_names = []
    for arg in args:
        arg_names.append([i for i, j in locals().items() if j == arg][0])

    # Convert all arguments to list
    args = convert_to_list(args)

    # Define metadata filter
    if is_primary_data:
        obs_value_filter = f"is_primary_data == True"
    else:
        obs_value_filter = None
    for arg_name, arg in zip(arg_names, args):
        if arg:
            if obs_value_filter is None:
                obs_value_filter = f"{arg_name} in {str(arg)}"
            else:
                obs_value_filter = obs_value_filter + f" and {arg_name} in {str(arg)}"

    # Fetch AnnData object
    if not meta_only:
        if verbose:
            logging.info(
                "Fetching AnnData object from CZ CELLxGENE Discover. This might take a few minutes..."
            )

        if gene:
            var_value_filter = f"{'feature_id' if ensembl else 'feature_name'} in {gene}"
        else:
            var_value_filter = None

        with cellxgene_census.open_soma(census_version=census_version) as census:
            adata = cellxgene_census.get_anndata(
                census=census,
                organism=species,
                var_value_filter=var_value_filter,
                obs_value_filter=obs_value_filter,
                column_names={"obs": column_names},
            )

            if out:
                adata.write(out)

            return adata

    # Fetch metadata
    else:
        if verbose:
            logging.info("Fetching metadata from CZ CELLxGENE Discover...")
        with cellxgene_census.open_soma(census_version=census_version) as census:
            # Reads SOMADataFrame as a slice
            cell_metadata = census["census_data"][species].obs.read(
                value_filter=obs_value_filter, column_names=column_names
            )

            # Concatenates results to pyarrow.Table
            cell_metadata = cell_metadata.concat()

            # Converts to pandas.DataFrame
            cell_metadata = cell_metadata.to_pandas()

            if out:
                cell_metadata.to_csv(out, index=False)

            return cell_metadata
