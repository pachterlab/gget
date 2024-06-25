> Python arguments are equivalent to long-option arguments (`--arg`), unless otherwise specified. Flags are True/False arguments in Python. The manual for any gget tool can be called from the command-line using the `-h` `--help` flag.  
## gget opentargets 📖
Fetch associated diseases from [OpenTargets](https://platform.opentargets.org/) using Ensembl IDs.   
Return format: JSON/CSV (command-line) or data frame (Python).

**Positional argument**  
`ens_id`  
Ensembl gene ID, e.g ENSG00000169194.

**Optional arguments**  
`-r` `--resource`   
Defines the type of information to return in the output. Default: 'diseases'.   
Possible resources are:     
'diseases' - Returns associated diseases.

`-l` `--limit`  
Limit the number of results, e.g 10. Default: No limit.

`-o` `--out`    
Path to the JSON file the results will be saved in, e.g. path/to/directory/results.json. Default: Standard out.  
Python: `save=True` will save the output in the current working directory.

**Flags**   
`-csv` `--csv`  
Command-line only. Returns the output in CSV format, instead of JSON format.

`-q` `--quiet`   
Command-line only. Prevents progress information from being displayed.  
Python: Use `verbose=False` to prevent progress information from being displayed. 

`wrap_text`  
Python only. `wrap_text=True` displays data frame with wrapped text for easy reading (default: False).  
  
  
### Example
```bash
gget opentargets ENSG00000169194 -l 5
```
```python
# Python
import gget
gget.opentargets('ENSG00000169194', limit=5)
```
&rarr; Returns the top 5 diseases associated with the gene ENSG00000169194.

| id            | name          | description                                                           | score         |
|---------------|---------------|-----------------------------------------------------------------------|---------------|
| EFO_0000274   | atopic eczema | A chronic inflammatory genetically determined disease of the skin ... | 0.66364347241831 |
    
#### [More examples](https://github.com/pachterlab/gget_examples)
