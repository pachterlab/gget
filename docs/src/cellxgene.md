> Python arguments are equivalent to long-option arguments (`--arg`), unless otherwise specified. Flags are True/False arguments in Python.  The manual for any gget tool can be called from the command-line using the `-h` `--help` flag.  
## gget cellxgene 🍱  
Query data from [CZ CELLxGENE Discover](https://cellxgene.cziscience.com/) using the [CZ CELLxGENE Discover Census](https://github.com/chanzuckerberg/cellxgene-census).  

Returns: An AnnData object containing the count matrix and metadata of single-cell RNA sequencing data from the defined tissues/genes/etc.  

Before using `gget cellxgene` for the first time, run `gget setup cellxgene` / `gget.setup("cellxgene")` once (also see [`gget setup`](setup.md)).  

**Optional arguments**  
`-s` `--species`  
Choice of 'homo_sapiens' or 'mus_musculus'. Default: 'homo_sapiens'.  

`-g` `--gene`  
 Str or list of gene name(s) or Ensembl ID(s). Default: None.  
 NOTE: Use `-e / --ensembl` (Python: `ensembl=True`) when providing Ensembl ID(s) instead of gene name(s).  
 See https://cellxgene.cziscience.com/gene-expression for examples of available genes.  

 `-cv` `--census_version`  
 Str defining version of Census, e.g. "2023-05-08", or "latest" or "stable". Default: "stable".  

`-cn` `--column_names`  
List of metadata columns to return (stored in AnnData.obs).  
Default: ['dataset_id', 'assay', 'suspension_type', 'sex', 'tissue_general', 'tissue', 'cell_type']  
For more options see: https://api.cellxgene.cziscience.com/curation/ui/#/ -> Schemas -> dataset  

`-o` `--out`   
Path to file to save generated AnnData .h5ad file (or .csv with `-mo / --meta_only` (`anndata=False`)).  
Required when using from command line!  

**Flags**  
`-e` `--ensembl`  
Use when genes are provided as Ensembl IDs instead of gene names.  

`-mo` `--meta_only`  
Only returns metadata dataframe (corresponds to AnnData.obs).  

`-q` `--quiet`   
Command-line only. Prevents progress information from being displayed.  
Python: Use `verbose=False` to prevent progress information from being displayed.  

**Optional arguments corresponding to CZ CELLxGENE Discover metadata attributes**  
`--tissue`  
Str or list of tissue(s), e.g. ['lung', 'blood']. Default: None.  
See https://cellxgene.cziscience.com/gene-expression for examples of available tissues.  

`--cell_type`  
Str or list of celltype(s), e.g. ['mucus secreting cell', 'neuroendocrine cell']. Default: None.  
See https://cellxgene.cziscience.com/gene-expression and select a tissue to see examples of available celltypes.  

`--development_stage`  
Str or list of development stage(s). Default: None.  

`--disease`  
Str or list of disease(s). Default: None.  

`--sex`  
Str or list of sex(es), e.g. 'female'. Default: None.  

`--dataset_id`  
Str or list of CELLxGENE dataset ID(s). Default: None.  

`--tissue_general_ontology_term_id`  
Str or list of high-level tissue UBERON ID(s). Default: None.  
Tissue labels and their corresponding UBERON IDs are listed [here](https://github.com/chanzuckerberg/single-cell-data-portal/blob/9b94ccb0a2e0a8f6182b213aa4852c491f6f6aff/backend/wmg/data/tissue_mapper.py).  

`--tissue_general`  
Str or list of high-level tissue label(s). Default: None.  
Tissue labels and their corresponding UBERON IDs are listed [here](https://github.com/chanzuckerberg/single-cell-data-portal/blob/9b94ccb0a2e0a8f6182b213aa4852c491f6f6aff/backend/wmg/data/tissue_mapper.py).  

`--tissue_ontology_term_id`  
Str or list of tissue ontology term ID(s) as defined in the CELLxGENE dataset schema. Default: None.  

`--assay_ontology_term_id`  
Str or list of assay ontology term ID(s) as defined in the CELLxGENE dataset schema. Default: None.  

`--assay`  
Str or list of assay(s) as defined in the CELLxGENE dataset schema. Default: None.  

`--cell_type_ontology_term_id`  
Str or list of celltype ontology term ID(s) as defined in the CELLxGENE dataset schema. Default: None.  

`--development_stage_ontology_term_id`  
Str or list of development stage ontology term ID(s) as defined in the CELLxGENE dataset schema. Default: None.  

`--disease_ontology_term_id`  
Str or list of disease ontology term ID(s) as defined in the CELLxGENE dataset schema. Default: None.  

`--donor_id`  
Str or list of donor ID(s) as defined in the CELLxGENE dataset schema. Default: None.  

`--self_reported_ethnicity_ontology_term_id`  
Str or list of self reported ethnicity ontology ID(s) as defined in the CELLxGENE dataset schema. Default: None.  

`--self_reported_ethnicity`  
Str or list of self reported ethnicity as defined in the CELLxGENE dataset schema. Default: None.  

`--sex_ontology_term_id`  
Str or list of sex ontology ID(s) as defined in the CELLxGENE dataset schema. Default: None.  

`--suspension_type`  
Str or list of suspension type(s) as defined in the CELLxGENE dataset schema. Default: None.  

  
### Examples
```bash
gget cellxgene --gene ACE2 ABCA1 SLC5A1 --tissue lung --cell_type 'mucus secreting cell' 'neuroendocrine cell' -o example_adata.h5ad
```
```python
# Python
adata = gget.cellxgene(
    gene = ["ACE2", "ABCA1", "SLC5A1"],
    tissue = "lung",
    cell_type = ["mucus secreting cell", "neuroendocrine cell"]
)
adata
```
&rarr; Returns an AnnData object containing the scRNAseq ACE2, ABCA1, and SLC5A1 count matrix of 3322 human lung mucus secreting and neuroendocrine cells from CZ CELLxGENE Discover and their corresponding metadata.  

Fetch metadata (corresponds to AnnData.obs) only:  
```bash
gget cellxgene --meta_only --gene ENSMUSG00000015405 --ensembl --tissue lung --species mus_musculus -o example_meta.csv
```
```python
# Python
df = gget.cellxgene(
    meta_only = True,
    gene = "ENSMUSG00000015405",
    ensembl = True,
    tissue = "lung",  
    species = "mus_musculus"
)
df
```
&rarr; Returns only the metadata from ENSMUSG00000015405 (ACE2) expression datasets corresponding to mouse lung cells.  
