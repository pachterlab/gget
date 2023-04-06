> Python arguments are equivalent to long-option arguments (`--arg`), unless otherwise specified. Flags are True/False arguments in Python. The manual for any gget tool can be called from the command-line using the `-h` `--help` flag.  
## gget info 💡
Fetch extensive gene and transcript metadata from [Ensembl](https://www.ensembl.org/), [UniProt](https://www.uniprot.org/), and [NCBI](https://www.ncbi.nlm.nih.gov/) using Ensembl IDs.  
Return format: JSON (command-line) or data frame/CSV (Python).

**Positional argument**  
`ens_ids`   
One or more Ensembl IDs.

**Optional arguments**  
`-o` `--out`   
Path to the file the results will be saved in, e.g. path/to/directory/results.csv (or .json). Default: Standard out.    
Python: `save=True` will save the output in the current working directory.

**Flags**  
`-n` `--ncbi`  
TURN OFF results from [NCBI](https://www.ncbi.nlm.nih.gov/).  
Python: `ncbi=False` prevents data retrieval from NCBI (default: True).    

`-u` `--uniprot`  
TURN OFF results from [UniProt](https://www.uniprot.org/).  
Python: `uniprot=False` prevents data retrieval from UniProt (default: True).   

`-pdb` `--pdb`  
INCLUDE [PDB](https://www.ebi.ac.uk/pdbe/) IDs in output (might increase runtime).  
Python: `pdb=False` prevents data retrieval from PDB (default: False).   

`-q` `--quiet`   
Command-line only. Prevents progress information from being displayed.  
Python: Use `verbose=False` to prevent progress information from being displayed.  

`-csv` `--csv`  
Command-line only. Returns results in CSV format.  
Python: Use `json=True` to return output in JSON format.

`wrap_text`  
Python only. `wrap_text=True` displays data frame with wrapped text for easy reading (default: False).  


### Example
```bash
gget info ENSG00000034713 ENSG00000104853 ENSG00000170296
```
```python
# Python
gget.info(["ENSG00000034713", "ENSG00000104853", "ENSG00000170296"])
```
&rarr; Returns extensive information about each requested Ensembl ID:  

|      | uniprot_id     | ncbi_gene_id     | primary_gene_name | synonyms | protein_names | ensembl_description | uniprot_description | ncbi_description | biotype | canonical_transcript | ... |
| -------------- |-------------------------| ------------------------| -------------- | ----------|-----|----|----|----|----|----|----|
| ENSG00000034713| P60520 | 11345 | GABARAPL2 | [ATG8, ATG8C, FLC3A, GABARAPL2, GATE-16, GATE16, GEF-2, GEF2] | Gamma-aminobutyric acid receptor-associated protein like 2 (GABA(A) receptor-associated protein-like 2)... | GABA type A receptor associated protein like 2 [Source:HGNC Symbol;Acc:HGNC:13291] | FUNCTION: Ubiquitin-like modifier involved in intra- Golgi traffic (By similarity). Modulates intra-Golgi transport through coupling between NSF activity and ... | Enables ubiquitin protein ligase binding activity. Involved in negative regulation of proteasomal protein catabolic process and protein... | protein_coding | ENST00000037243.7 |... |
| . . .            | . . .                     | . . .                     | . . .            | . . .       | . . . | . . . | . . . | . . . | . . . | . . . | ... |
  
#### [More examples](https://github.com/pachterlab/gget_examples)
