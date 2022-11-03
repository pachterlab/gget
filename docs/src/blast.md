Python arguments are equivalent to long-option arguments (`--arg`), unless otherwise specified. Flags are True/False arguments in Python.  
The manual for any gget tool can be called from the command-line using the `-h` `--help` flag.  
  
## gget blast 💥
BLAST a nucleotide or amino acid sequence to any [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) database.  
Return format: JSON (command-line) or data frame/CSV (Python).

**Positional argument**  
`sequence`   
Nucleotide or amino acid sequence, or path to FASTA or .txt file.

**Optional arguments**  
`-p` `--program`  
'blastn', 'blastp', 'blastx', 'tblastn', or 'tblastx'.  
Default: 'blastn' for nucleotide sequences; 'blastp' for amino acid sequences.

`-db` `--database`  
'nt', 'nr', 'refseq_rna', 'refseq_protein', 'swissprot', 'pdbaa', or 'pdbnt'.  
Default: 'nt' for nucleotide sequences; 'nr' for amino acid sequences.  
[More info on BLAST databases](https://ncbi.github.io/blast-cloud/blastdb/available-blastdbs.html)

`-l` `--limit`  
Limits number of hits to return. Default: 50.  

`-e` `--expect`  
Defines the [expect value](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=FAQ#expect) cutoff. Default: 10.0.  

`-o` `--out`   
Path to the file the results will be saved in, e.g. path/to/directory/results.csv (or .json). Default: Standard out.   
Python: `save=True` will save the output in the current working directory.

**Flags**  
`-lcf` `--low_comp_filt`  
Turns on [low complexity filter](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=FAQ#LCR).  

`-mbo` `--megablast_off`  
Turns off MegaBLAST algorithm. Default: MegaBLAST on (blastn only).  

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
gget blast MKWMFKEDHSLEHRCVESAKIRAKYPDRVPVIVEKVSGSQIVDIDKRKYLVPSDITVAQFMWIIRKRIQLPSEKAIFLFVDKTVPQSR
```
```python
# Python
gget.blast("MKWMFKEDHSLEHRCVESAKIRAKYPDRVPVIVEKVSGSQIVDIDKRKYLVPSDITVAQFMWIIRKRIQLPSEKAIFLFVDKTVPQSR")
```
&rarr; Returns the BLAST result of the sequence of interest. `gget blast` automatically detects this sequence as an amino acid sequence and therefore sets the BLAST program to *blastp* with database *nr*.  

| Description     | Scientific Name	     | Common Name     | Taxid        | Max Score | Total Score | Query Cover | ... |
| -------------- |-------------------------| ------------------------| -------------- | ----------|-----|---|---|
| PREDICTED: gamma-aminobutyric acid receptor-as...| Colobus angolensis palliatus	 | 	NaN | 336983 | 180	 | 180 | 100% | ... |
| . . . | . . . | . . . | . . . | . . . | . . . | . . . | ... | 

BLAST from .fa or .txt file:  
```bash
gget blast fasta.fa
```
```python
# Python
gget.blast("fasta.fa")
```
&rarr; Returns the BLAST results of the first sequence contained in the fasta.fa file. 

#### [More examples](https://github.com/pachterlab/gget_examples)
