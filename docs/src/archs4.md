> Python arguments are equivalent to long-option arguments (`--arg`), unless otherwise specified. Flags are True/False arguments in Python. The manual for any gget tool can be called from the command-line using the `-h` `--help` flag.  
## gget archs4 🐁
Find the most correlated genes to a gene of interest or find the gene's tissue expression atlas using [ARCHS4](https://maayanlab.cloud/archs4/).  
Return format: JSON (command-line) or data frame/CSV (Python).

**Positional argument**  
`gene`  
Short name (gene symbol) of gene of interest, e.g. STAT4.  
Alternatively: use flag `--ensembl` to input an Ensembl gene IDs, e.g. ENSG00000138378.

**Optional arguments**  
 `-w` `--which`  
'correlation' (default) or 'tissue'.  
'correlation' returns a gene correlation table that contains the 100 most correlated genes to the gene of interest. The Pearson correlation is calculated over all samples and tissues in [ARCHS4](https://maayanlab.cloud/archs4/).  
'tissue' returns a tissue expression atlas calculated from human or mouse samples (as defined by 'species') in [ARCHS4](https://maayanlab.cloud/archs4/).  

`-s` `--species`  
'human' (default) or 'mouse'.   
Defines whether to use human or mouse samples from [ARCHS4](https://maayanlab.cloud/archs4/).  
(Only for tissue expression atlas.)

`-o` `--out`   
Path to the file the results will be saved in, e.g. path/to/directory/results.csv (or .json). Default: Standard out.   
Python: `save=True` will save the output in the current working directory.  
  
**Flags**   
`-e` `--ensembl`  
Add this flag if `gene` is given as an Ensembl gene ID.  

`-csv` `--csv`  
Command-line only. Returns results in CSV format.  
Python: Use `json=True` to return output in JSON format.
  
  
### Examples
```bash
gget archs4 ACE2
```
```python
# Python
gget.archs4("ACE2")
```
&rarr; Returns the 100 most correlated genes to ACE2:  

| gene_symbol     | pearson_correlation     |
| -------------- |-------------------------| 
| SLC5A1 | 0.579634 | 	
| CYP2C18 | 0.576577 | 	
| . . . | . . . | 	
  

```bash
gget archs4 -w tissue ACE2
```
```python
# Python
gget.archs4("ACE2", which="tissue")
```
&rarr; Returns the tissue expression of ACE2 (by default, human data is used):

| id     | min     | q1 |  median | q3 | max |
| ------ |--------| ------ |--------| ------ |--------| 
| System.Urogenital/Reproductive System.Kidney.RENAL CORTEX | 0.113644 | 8.274060 | 9.695840 | 10.51670 | 11.21970 |
| System.Digestive System.Intestine.INTESTINAL EPITHELIAL CELL | 0.113644 | 	5.905560 | 9.570450 | 13.26470 | 13.83590 | 
| . . . | . . . | . . . | . . . | . . . | . . . |

#### [More examples](https://github.com/pachterlab/gget_examples)
