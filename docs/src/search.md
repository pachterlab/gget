> Python arguments are equivalent to long-option arguments (`--arg`), unless otherwise specified. Flags are True/False arguments in Python. The manual for any gget tool can be called from the command-line using the `-h` `--help` flag.  
## gget search 🔎
Fetch genes and transcripts from [Ensembl](https://www.ensembl.org/) using free-form search terms.   
Note: Only returns results based on matches in the "gene name" or "description" sections in the Ensembl database.  
Return format: JSON (command-line) or data frame/CSV (Python).

**Positional argument**
`searchwords`   
One or more free form search words, e.g. gaba nmda. (Note: Search is not case-sensitive.)

**Other required arguments**   
`-s` `--species`  
Species or database to be searched.  
A species can be passed in the format 'genus_species', e.g. 'homo_sapiens' or 'arabidopsis_thaliana'.  
To pass a specific database, pass the name of the CORE database, e.g. 'mus_musculus_dba2j_core_105_1'.  
All available databases for each Ensembl release can be found [here](http://ftp.ensembl.org/pub/).  
  
Supported shortcuts: 'human', 'mouse'. 

**Optional arguments**  
`-r` `--release`   
Defines the Ensembl release number from which the files are fetched, e.g. 104. Default: None -> latest Ensembl release is used.  
Note: *Does not apply to plant species* (you can pass a specific plant core database (which includes a release number) to the `species` argument instead).  
This argument is overwritten if a specific database (which includes a release number) is passed to the species argument.   

`-t` `--id_type`  
'gene' (default) or 'transcript'  
Returns genes or transcripts, respectively.

`-ao` `--andor`  
'or' (default) or 'and'  
'or': Returns all genes that INCLUDE AT LEAST ONE of the searchwords in their name/description.  
'and': Returns only genes that INCLUDE ALL of the searchwords in their name/description.

`-l` `--limit`   
Limits the number of search results, e.g. 10. Default: None.  

`-o` `--out`  
Path to the csv the results will be saved in, e.g. path/to/directory/results.csv (or .json). Default: Standard out.   
Python: `save=True` will save the output in the current working directory.

**Flags**  
`-csv` `--csv`  
Command-line only. Returns results in CSV format.  
Python: Use `json=True` to return output in JSON format.

`-q` `--quiet`   
Command-line only. Prevents progress information from being displayed.  
Python: Use `verbose=False` to prevent progress information from being displayed. 

`wrap_text`  
Python only. `wrap_text=True` displays data frame with wrapped text for easy reading (default: False).  
 
    
    
### Example
```bash
gget search -s human gaba gamma-aminobutyric
```
```python
# Python
gget.search(["gaba", "gamma-aminobutyric"], "homo_sapiens")
```
&rarr; Returns all genes that contain at least one of the search words in their name or Ensembl/external reference description:

| ensembl_id     | gene_name     | ensembl_description     | ext_ref_description        | biotype | url |
| -------------- |-------------------------| ------------------------| -------------- | ----------|-----|
| ENSG00000034713| GABARAPL2 | 	GABA type A receptor associated protein like 2 [Source:HGNC Symbol;Acc:HGNC:13291] | GABA type A receptor associated protein like 2 | protein_coding | https://uswest.ensembl.org/homo_sapiens/Gene/Summary?g=ENSG00000034713 |
| . . .            | . . .                     | . . .                     | . . .            | . . .       | . . . |
    
#### [More examples](https://github.com/pachterlab/gget_examples)
