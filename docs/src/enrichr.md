> Python arguments are equivalent to long-option arguments (`--arg`), unless otherwise specified. Flags are True/False arguments in Python. The manual for any gget tool can be called from the command-line using the `-h` `--help` flag.  
## gget enrichr 💰
Perform an enrichment analysis on a list of genes using [Enrichr](https://maayanlab.cloud/Enrichr/).  
Return format: JSON (command-line) or data frame/CSV (Python).
  
**Positional argument**  
`genes`  
Short names (gene symbols) of genes to perform enrichment analysis on, e.g. PHF14 RBM3 MSL1 PHF21A.  
Alternatively: use flag `--ensembl` to input a list of Ensembl gene IDs, e.g. ENSG00000106443 ENSG00000102317 ENSG00000188895.

**Other required arguments**  
`-db` `--database`  
Database to use as reference for the enrichment analysis.  
Supports any database listed [here](https://maayanlab.cloud/Enrichr/#libraries) under 'Gene-set Library' or one of the following shortcuts:  
'pathway'       (KEGG_2021_Human)  
'transcription'     (ChEA_2016)  
'ontology'      (GO_Biological_Process_2021)  
'diseases_drugs'   (GWAS_Catalog_2019)   
'celltypes'      (PanglaoDB_Augmented_2021)  
'kinase_interactions'   (KEA_2015)  
  
**Optional arguments**  
`-o` `--out`   
Path to the file the results will be saved in, e.g. path/to/directory/results.csv (or .json). Default: Standard out.   
Python: `save=True` will save the output in the current working directory.

`figsize`  
Python only. (width, height) of plot in inches. (Default: (10,10))

`ax`  
Python only. Pass a matplotlib axes object for plot customization. (Default: None)
  
**Flags**  
`-e` `--ensembl`  
Add this flag if `genes` are given as Ensembl gene IDs.  
 
`-csv` `--csv`  
Command-line only. Returns results in CSV format.  
Python: Use `json=True` to return output in JSON format.
  
`plot`  
Python only. `plot=True` provides a graphical overview of the first 15 results (default: False).  
  
  
### Example
```bash
gget enrichr -db ontology ACE2 AGT AGTR1
```
```python
# Python
gget.enrichr(["ACE2", "AGT", "AGTR1"], database="ontology", plot=True)
```
&rarr; Returns pathways/functions involving genes ACE2, AGT, and AGTR1 from the *GO Biological Process 2021* database. In Python, `plot=True` returns a graphical overview of the results:

![alt text](https://github.com/pachterlab/gget/blob/main/figures/gget_enrichr_results.png?raw=true)

#### [More examples](https://github.com/pachterlab/gget_examples)
