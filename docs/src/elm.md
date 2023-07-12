> Python arguments are equivalent to long-option arguments (`--arg`), unless otherwise specified. Flags are True/False arguments in Python. The manual for any gget tool can be called from the command-line using the `-h` `--help` flag.  
## gget elm 💡
Fetch functional sites information from [ELM](http://elm.eu.org/) using an amino acid sequence or Uniprot ID.  
NOTE:  
Please limit your searches to a maximum of 1 per minute for amino acid sequences (1 per 3 minutes for Uniprot IDs). If you exceed this limit, you will recieve a "429 Too many requests" error. Also please note that this does not always work for sequences longer than 2000 amino acids: URLs may be truncated beyond this length.

Return format: JSON (command-line) or data frame/CSV (Python).

**Positional argument**  
`sequence`  
Amino acid sequence or UniProt ID.  
Use `--uniprot` (Python: `uniprot=True`) if a UniProt ID is supplied.

**Optional arguments**  
`-o` `--out`   
Path to the file the results will be saved in, e.g. path/to/directory/results.csv (or .json). Default: Standard out.   
Python: `save=True` will save the output in the current working directory.  

**Flags**   

`-u` `--uniprot`  
Search using Uniprot ID instead of amino acid sequence  

`-csv` `--csv`  
Command-line only. Returns results in CSV format.  
Python: Use `json=True` to return output in JSON format.
  
`-q` `--quiet`   
Command-line only. Prevents progress information from being displayed.  
Python: Use `verbose=False` to prevent progress information from being displayed. 


### Example
```bash
gget elm DRVYVHPFHL
```
```python
# Python
gget.info("DRVYVHPFHL")
```
&rarr; Returns functional site information about each elm identifier:  

|      | elm_identifier     | Accession     | Functional site class | Functional site description | ELM Description | ELMs with same func. site | Pattern | Pattern Probability | Present in taxons | Interaction Domain | ... |
| -------------- |-------------------------| ------------------------| -------------- | ----------|-----|----|----|----|----|----|----|
| LIG_Arc_Nlobe_1| ELME000534 | Arc N-lobe binding ligand | The activity-regulated cytoskeleton-associated... | The motif peptide binds Arc in an unusual conf... | NaN| [^P][P]G{0,1}[^P][YFH][^P] | 0.0043852 |Mammalia Tetrapoda| Activity-regulated cytoskeleton-associated pro... |... |
| . . .            | . . .                     | . . .                     | . . .            | . . .       | . . . | . . . | . . . | . . . | . . . | . . . | ... |
  

#### [More examples](https://github.com/pachterlab/gget_examples)
