> Python arguments are equivalent to long-option arguments (`--arg`), unless otherwise specified. Flags are True/False arguments in Python. The manual for any gget tool can be called from the command-line using the `-h` `--help` flag.  
## gget elm 💡
Fetch functional sites information from [ELM](http://elm.eu.org/) using an amino acid sequence or Uniprot ID.  
Return format: JSON (command-line) or data frame/CSV (Python).

**Positional argument**  
`-seq` `--sequence`  
Amino acid sequence or Uniprot ID.

**Optional arguments**  

**Flags**   

`-u` `--uniprot`  
Search using Uniprot ID instead of amino acid sequence 
Python: `uniprot=True` turns on search using Uniprot ID (default: False).   
  

### Example
```bash
gget elm DRVYVHPFHL
```
```python
# Python
gget.info("DRVYVHPFHL")
```
&rarr; Returns extensive information about each requested Ensembl ID:  

|      | elm_identifier     | Accession     | Functional site class | Functional site description | ELM Description | ELMs with same func. site | Pattern | Pattern Probability | Present in taxons | Interaction Domain | ... |
| -------------- |-------------------------| ------------------------| -------------- | ----------|-----|----|----|----|----|----|----|
| LIG_Arc_Nlobe_1| ELME000534 | Arc N-lobe binding ligand | The activity-regulated cytoskeleton-associated... | The motif peptide binds Arc in an unusual conf... | NaN| [^P][P]G{0,1}[^P][YFH][^P] | 0.0043852 |Mammalia Tetrapoda| Activity-regulated cytoskeleton-associated pro... |... |
| . . .            | . . .                     | . . .                     | . . .            | . . .       | . . . | . . . | . . . | . . . | . . . | . . . | ... |
  
#### [More examples](https://github.com/pachterlab/gget_examples)
