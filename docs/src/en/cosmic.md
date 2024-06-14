> Python arguments are equivalent to long-option arguments (`--arg`), unless otherwise specified. Flags are True/False arguments in Python.  The manual for any gget tool can be called from the command-line using the `-h` `--help` flag.  
## gget cosmic 🪐
Search for genes, mutations, and other factors associated with cancer using the [COSMIC](https://cancer.sanger.ac.uk/cosmic) (Catalogue Of Somatic Mutations In Cancer) database.  
Return format: JSON (command-line) or data frame/CSV (Python) when `download_cosmic=False`. When `download_cosmic=True`, downloads the requested database into the specified folder.    

This module was written in part by [@AubakirovArman](https://github.com/AubakirovArman) (information querying) and [@josephrich98](https://github.com/josephrich98) (database download).  

NOTE: License fees apply for the commercial use of COSMIC. You can read more about licensing COSMIC data [here](https://cancer.sanger.ac.uk/cosmic/license).

**Positional argument (for querying information)**  
`searchterm`   
Search term, which can be a mutation, or gene name (or Ensembl ID), or sample, etc.  
Examples for the searchterm and entitity arguments:   

| searchterm   | entitity    | |
|--------------|-------------| ---|
| EGFR         | mutations   | -> Find mutations in the EGFR gene that are associated with cancer |
| v600e        | mutations   | -> Find genes for which a v600e mutation is associated with cancer |
| COSV57014428 | mutations   | -> Find mutations associated with this COSMIC mutations ID |
| EGFR         | genes       | -> Get the number of samples, coding/simple mutations, and fusions observed in COSMIC for EGFR |
| prostate     | cancer      | -> Get number of tested samples and mutations for prostate cancer |
| prostate     | tumour_site | -> Get number of tested samples, genes, mutations, fusions, etc. with 'prostate' as primary tissue site |
| ICGC         | studies     | -> Get project code and descriptions for all studies from the ICGC (International Cancer Genome Consortium) |
| EGFR         | pubmed      | -> Find PubMed publications on EGFR and cancer |
| ICGC         | samples     | -> Get metadata on all samples from the ICGC (International Cancer Genome Consortium) |
| COSS2907494  | samples     | -> Get metadata on this COSMIC sample ID (cancer type, tissue, # analyzed genes, # mutations, etc.) |

NOTE: (Python only) Set to `None` when downloading COSMIC databases with `download_cosmic=True`.

**Optional arguments (for querying information)**  
`-e` `--entity`  
'mutations' (default), 'genes', 'cancer', 'tumour site', 'studies', 'pubmed', or 'samples'.  
Defines the type of the results to return. 

`-l` `--limit`  
Limits number of hits to return. Default: 100.  

**Flags (for downloading COSMIC databases)**  
`-d` `--download_cosmic`  
Switches into database download mode.  

`-gm` `--gget_mutate`  
TURNS OFF creation of a modified version of the database for use with gget mutate.  
Python: `gget_mutate` is True by default. Set `gget_mutate=False` to disable.  

**Optional arguments (for downloading COSMIC databases)**  
`-mc` `--mutation_class`
'cancer' (default), 'cell_line', 'census', 'resistance', 'screen', or 'cancer_example'  
Type of COSMIC database to download.  

`-cv` `--cosmic_version`  
Version of the COSMIC database. Default: None -> Defaults to latest version.  

`-gv` `--grch_version`  
Version of the human GRCh reference genome the COSMIC database was based on (37 or 38). Default: 37  

**Optional arguments (general)**  
`-o` `--out`   
Path to the file (or folder when downloading databases with the `download_cosmic` flag) the results will be saved in, e.g. 'path/to/results.json'.  
Default: None  
-> When download_cosmic=False: Results will be returned to standard out  
-> When download_cosmic=True: Database will be downloaded into current working directory  

**Flags (general)**  
`-csv` `--csv`  
Command-line only. Returns results in CSV format.  
Python: Use `json=True` to return output in JSON format.

`-q` `--quiet`   
Command-line only. Prevents progress information from being displayed.  
Python: Use `verbose=False` to prevent progress information from being displayed.  

  
### Examples
#### Query information
```bash
gget cosmic EGFR
```
```python
# Python
gget.cosmic("EGFR")
```
&rarr; Returns mutations in the EGFR gene that are associated with cancer in the format:

| Gene     | Syntax     | Alternate IDs                  | Canonical  |
| -------- |------------| -------------------------------| ---------- |
| EGFR     | c.*2446A>G | EGFR c.*2446A>G, EGFR p.?, ... | y          |
| EGFR     | c.(2185_2283)ins(18) | EGFR c.(2185_2283)ins(18), EGFR p.?, ... | y          |
| . . .    | . . .      | . . .                          | . . .      | 


### Downloading COSMIC databases
```bash
gget cosmic --download_cosmic
```
```python
# Python
gget.cosmic(searchterm=None, download_cosmic=True)
```
&rarr; Downloads the COSMIC cancer database of the latest COSMIC release into the current working directory.

