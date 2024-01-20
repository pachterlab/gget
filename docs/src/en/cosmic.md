> Python arguments are equivalent to long-option arguments (`--arg`), unless otherwise specified. Flags are True/False arguments in Python.  The manual for any gget tool can be called from the command-line using the `-h` `--help` flag.  
## gget cosmic 🪐
Search for genes, mutations, etc associated with cancers using the [COSMIC](https://cancer.sanger.ac.uk/cosmic) (Catalogue Of Somatic Mutations In Cancer) database.  
Return format: JSON (command-line) or data frame/CSV (Python).

**Positional argument**  
`searchterm`   
Search term, which can be a mutation, or gene, or sample, etc. as defined using the `entity` argument. Example: 'EGFR'  

**Optional arguments**  
`-e` `--entity`  
'mutations' (default), 'genes', 'cancer', 'tumour site', 'studies', 'pubmed', or 'samples'.  
Defines the type of the supplied search term.  

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

