> Python arguments are equivalent to long-option arguments (`--arg`), unless otherwise specified. Flags are True/False arguments in Python. The manual for any gget tool can be called from the command-line using the `-h` `--help` flag.  
## gget seq 🧬
Fetch nucleotide or amino acid sequence of a gene (and all its isoforms) or a transcript by Ensembl ID.   
Return format: FASTA.

**Positional argument**  
`ens_ids`   
One or more Ensembl IDs.

**Optional arguments**  
`-o` `--out`   
Path to the file the results will be saved in, e.g. path/to/directory/results.fa. Default: Standard out.   
Python: `save=True` will save the output in the current working directory.

**Flags**  
`-t` `--translate`  
Returns amino acid (instead of nucleotide) sequences.  
Nucleotide sequences are fetched from [Ensembl](https://www.ensembl.org/).  
Amino acid sequences are fetched from [UniProt](https://www.uniprot.org/).

`-iso` `--isoforms`   
Returns the sequences of all known transcripts.  
(Only for gene IDs.)


### Examples  
```bash
gget seq ENSG00000034713 ENSG00000104853 ENSG00000170296
```
```python
# Python
gget.seq(["ENSG00000034713", "ENSG00000104853", "ENSG00000170296"])
```
&rarr; Returns the nucleotide sequences of ENSG00000034713, ENSG00000104853, and ENSG00000170296 in FASTA format.  


```bash
gget seq -t -iso ENSG00000034713
```
```python
# Python
gget.seq("ENSG00000034713", translate=True, isoforms=True)
```
&rarr; Returns the amino acid sequences of all known transcripts of ENSG00000034713 in FASTA format.

#### [More examples](https://github.com/pachterlab/gget_examples)
