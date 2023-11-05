> Python arguments are equivalent to long-option arguments (`--arg`), unless otherwise specified. Flags are True/False arguments in Python. The manual for any gget tool can be called from the command-line using the `-h` `--help` flag.  
## gget elm 🎭
Locally predict Eukaryotic Linear Motifs from an amino acid sequence or UniProt ID using data from the [ELM database](http://elm.eu.org/).    
Return format: JSON (command-line) or data frame/CSV (Python). This module resturns two data frames (or JSON formatted files) (see examples).     

Before using `gget elm` for the first time, run `gget setup elm` / `gget.setup("elm")` once (also see [`gget setup`](setup.md)).   

**Positional argument**  
`sequence`  
Amino acid sequence or Uniprot ID (str).  
When providing a Uniprot ID, use flag `--uniprot` (Python: `uniprot==True`).  

**Optional arguments**  
`-s` `--sensitivity`  
Sensitivity of DIAMOND alignment (str). Default: "very-sensitive".   
One of the following: fast, mid-sensitive, sensitive, more-sensitive, very-sensitive, or ultra-sensitive.  

`-t` `--threads`  
Number of threads used in DIAMOND alignment (int). Default: 1.  

`-db` `--diamond_binary`  
Path to DIAMOND binary (str). Default: None -> Uses DIAMOND binary installed with `gget`.  

`-o` `--out`   
Path to the folder to save results in, e.g. path/to/directory. Default: Standard out; temporary files are deleted.   

**Flags**  
`-u` `--uniprot`  
Set to True if `sequence` is a Uniprot ID instead of an amino acid sequence. Default: False.  

`-csv` `--csv`  
Command-line only. Returns results in CSV format.  
Python: Use `json=True` to return output in JSON format.

`-q` `--quiet`   
Command-line only. Prevents progress information from being displayed.  
Python: Use `verbose=False` to prevent progress information from being displayed.   

### Example
```bash
gget setup elm          # Downloads/updates local ELM database
gget elm -o gget_elm_results LIAQSIGQASFV
gget elm -o gget_elm_results --uniprot Q02410
```
```python
# Python
gget.setup(“elm”)      # Downloads/updates local ELM database
ortholog_df, regex_df = gget.elm(“LIAQSIGQASFV”)
ortholog_df, regex_df = gget.elm(“Q02410”, uniprot=True)
```
&rarr; Returns two data frames containing extensive information about linear motifs associated with orthologous proteins and motifs found in the input sequence directly based on their regex expressions.    

#### [More examples](https://github.com/pachterlab/gget_examples)
