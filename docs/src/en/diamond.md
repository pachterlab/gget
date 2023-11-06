> Python arguments are equivalent to long-option arguments (`--arg`), unless otherwise specified. Flags are True/False arguments in Python. The manual for any gget tool can be called from the command-line using the `-h` `--help` flag.  
## gget diamond 💎
Align multiple protein or translated DNA sequences using [DIAMOND](https://www.nature.com/articles/nmeth.3176) (DIAMOND is similar to BLAST, but this is a local computation).     
Return format: JSON (command-line) or data frame/CSV (Python).  

**Positional argument**  
`query`  
Sequences (str or list) or path to FASTA file containing sequences to be aligned against the reference.  

**Required arguments**  
`-ref` `--reference`  
Reference sequences (str or list) or path to FASTA file containing reference sequences.  

**Optional arguments**  
`-db` `--diamond_db`  
Path to save DIAMOND database created from `reference` (str).  
Default: None -> Temporary db file will be deleted after alignment or saved in `out` if `out` is provided.  

`-s` `--sensitivity`  
Sensitivity of alignment (str). Default: "very-sensitive".   
One of the following: fast, mid-sensitive, sensitive, more-sensitive, very-sensitive, or ultra-sensitive.  

`-t` `--threads`  
Number of threads used (int). Default: 1.  

`-db` `--diamond_binary`  
Path to DIAMOND binary (str). Default: None -> Uses DIAMOND binary installed with `gget`.  

`-o` `--out`   
Path to the folder to save results in (str), e.g. "path/to/directory". Default: Standard out; temporary files are deleted.   

**Flags**  
`-csv` `--csv`  
Command-line only. Returns results in CSV format.  
Python: Use `json=True` to return output in JSON format.  

`-q` `--quiet`   
Command-line only. Prevents progress information from being displayed.  
Python: Use `verbose=False` to prevent progress information from being displayed.   

### Example
```bash
# !!! Make sure to list the positional argument first here so it is not added as a reference sequence
gget diamond GGETISAWESQME ELVISISALIVE LQVEFRANKLIN PACHTERLABRQCKS -ref GGETISAWESQMEELVISISALIVELQVEFRANKLIN PACHTERLABRQCKS
```
```python
# Python
gget.diamond(["GGETISAWESQME", "ELVISISALIVE", "LQVEFRANKLIN", "PACHTERLABRQCKS"], reference=["GGETISAWESQMEELVISISALIVELQVEFRANKLIN", "PACHTERLABRQCKS"])
```
&rarr;    

#### [More examples](https://github.com/pachterlab/gget_examples)