> Python arguments are equivalent to long-option arguments (`--arg`), unless otherwise specified. Flags are True/False arguments in Python.  The manual for any gget tool can be called from the command-line using the `-h` `--help` flag.  
## gget cosmic 🪐
Search for genes, mutations, and other factors associated with cancer using the [COSMIC](https://cancer.sanger.ac.uk/cosmic) (Catalogue Of Somatic Mutations In Cancer) database.  
Return format: JSON (command-line) or data frame/CSV (Python).  
This module was written by [@AubakirovArman](https://github.com/AubakirovArman).

NOTE: License fees apply for the commercial use of COSMIC. You can read more about licensing COSMIC data [here](https://cancer.sanger.ac.uk/cosmic/license).

**Positional argument**  
`searchterm`   
Search term, which can be a mutation, or gene name (or Ensembl ID), or sample, etc.  
Examples for the searchterm and entitity arguments:   

| searchterm   | entitity    |
|--------------|-------------|
| EGFR         | mutations   | -> Find mutations in the EGFR gene that are associated with cancer
| v600e        | mutations   | -> Find genes for which a v600e mutation is associated with cancer
| COSV57014428 | mutations   | -> Find mutations associated with this COSMIC mutations ID
| EGFR         | genes       | -> Get the number of samples, coding/simple mutations, and fusions observed in COSMIC for EGFR
| prostate     | cancer      | -> Get number of tested samples and mutations for prostate cancer
| prostate     | tumour_site | -> Get number of tested samples, genes, mutations, fusions, etc. with 'prostate' as primary tissue site
| ICGC         | studies     | -> Get project code and descriptions for all studies from the ICGC (International Cancer Genome Consortium)
| EGFR         | pubmed      | -> Find PubMed publications on EGFR and cancer
| ICGC         | samples     | -> Get metadata on all samples from the ICGC (International Cancer Genome Consortium)
| COSS2907494  | samples     | -> Get metadata on this COSMIC sample ID (cancer type, tissue, # analyzed genes, # mutations, etc.)

NOTE: (Python only) Set to `None` when downloading COSMIC databases with `download_cosmic=True`.

**Optional arguments**  
`-e` `--entity`  
'mutations' (default), 'genes', 'cancer', 'tumour site', 'studies', 'pubmed', or 'samples'.  
Defines the type of the results to return. 

`-l` `--limit`  
Limits number of hits to return. Default: 100.  

`-o` `--out`   
Path to the file the results will be saved in, e.g. path/to/directory/results.csv (or .json). Default: Standard out.   
Python: `save=True` will save the output in the current working directory.  

**Flags**  
`-csv` `--csv`  
Command-line only. Returns results in CSV format.  
Python: Use `json=True` to return output in JSON format.

`-q` `--quiet`   
Command-line only. Prevents progress information from being displayed.  
Python: Use `verbose=False` to prevent progress information from being displayed.  

  
### Example
```bash
gget cosmic -e genes EGFR
```
```python
# Python
gget.cosmic("EGFR", entity="genes")
```
&rarr; Returns the COSMIC hits for gene 'EGFR' in the format:  

| Gene     | Alternate IDs     | Tested samples     | Simple Mutations        | Fusions | Coding Mutations | ... |
| -------------- |-------------------------| ------------------------| -------------- | ----------|-----|---|
| EGFR| EGFR,ENST00000275493.6,... | 210280 | 31900 | 0 | 31900 | ... |
| . . . | . . . | . . . | . . . | . . . | . . . | . . . | ... | 

