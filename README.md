# gget (gene-get)
[![pypi version](https://img.shields.io/pypi/v/gget)](https://pypi.org/project/gget/0.0.5/)
[![license](https://img.shields.io/pypi/l/gget)](LICENSE)
[![DOI](https://zenodo.org/badge/458943224.svg)](https://zenodo.org/badge/latestdoi/458943224)

gget features 3 main functionalities:  
- [**gget ref**](#gget-ref)  
Fetch genome references (GTF and FASTAs) from the [Ensembl FTP site](http://ftp.ensembl.org/pub/) for a specific species and release.
- [**gget search**](#gget-search)   
Query [Ensembl](https://www.ensembl.org/) for genes from a defined species using free form search words.
- [**gget lookup**](#gget-lookup)  
Look up gene or transcript Ensembl IDs for their common name, description, sequence, homologs, synonyms, corresponding transcript/gene and more from the Ensembl database as well as external references.

## Installation
```
pip install gget
```

For use in Jupyter Lab / Google Colab:
```python
from gget import ref, search, lookup
```

## gget ref
Function to fetch GTF and FASTA (cDNA and DNA) URLs from the [Ensembl FTP site](http://ftp.ensembl.org/pub/). Returns a dictionary/json containing the requested URLs with their respective Ensembl version and release date and time.

### Options
`-s` `--species`  
Species for which the FTPs will be fetched in the format genus_species, e.g. homo_sapiens.

`-w` `--which`  
Defines which results to return. Possible entries are:
'all' - Returns GTF, cDNA, and DNA links and associated info (default). 
Or one or a combination of the following:
'gtf' - Returns the GTF FTP link and associated info.
'cdna' - Returns the cDNA FTP link and associated info.
'dna' - Returns the DNA FTP link and associated info.

`-r` `--release`  
Ensemble release the FTPs will be fetched from, e.g. 104 (default: None &rarr; uses latest Ensembl release).

`-ftp` `--ftp`  
If True: returns only a list containing the requested FTP links (default: False).

`-o` `--out`  
Path to the file the results will be saved in, e.g. path/to/directory/results.json (default: None &rarr; just prints results).  
For Jupyter Lab / Google Colab: `save=True` will save the output to the current working directory.

### gget ref EXAMPLES
#### Fetch GTF, DNA, and cDNA FTP links for a specific species
Jupyter Lab / Google Colab:
```python
ref("homo_sapiens")
```

Terminal:
```
$ gget ref -s homo_sapiens
```
&rarr; Returns a json with the latest human reference genome GTF, DNA, and cDNA links, their respective release dates and time, and the Ensembl release from which the links were fetched, in the format:
```
{
            species: {
                "transcriptome_cdna": {
                    "ftp": cDNA FTP download URL,
                    "ensembl_release": Ensembl release,
                    "release_date": Day-Month-Year,
                    "release_time": HH:MM,
                    "bytes": cDNA FTP file size in bytes
                },
                "genome_dna": {
                    "ftp": DNA FTP download URL,
                    "ensembl_release": Ensembl release,
                    "release_date": Day-Month-Year,
                    "release_time": HH:MM,
                    "bytes": DNA FTP file size in bytes
                },
                "annotation_gtf": {
                    "ftp": GTF FTP download URL,
                    "ensembl_release": Ensembl release,
                    "release_date": Day-Month-Year,
                    "release_time": HH:MM,
                    "bytes": GTF FTP file size in bytes
                }
            }
        }
```

#### Fetch GTF, DNA, and cDNA FTP links for a specific species from a specific Ensembl release
For example, for Ensembl release 104:  
Jupyter Lab / Google Colab:
```python
ref("homo_sapiens", release=104)
```

Terminal:
```
$ gget ref -s homo_sapiens -r 104
```
&rarr; Returns a json with the human reference genome GTF, DNA, and cDNA links, and their respective release dates and time, from Ensembl release 104

#### Save the results
Jupyter Lab / Google Colab:
```python
ref("homo_sapiens", save=True)
```

Terminal:
```
$ gget ref -s homo_sapiens -o path/to/directory/ref_results.json
```
&rarr; Saves the results in path/to/directory/ref_results.json.  
For Jupyter Lab / Google Colab: Saves the results in a json file named ref_results.json in the current working directory.

#### Fetch only certain types of links for a specific species 
Jupyter Lab / Google Colab:
```python
ref("homo_sapiens", which=["gtf", "dna"])
```

Terminal:
```
$ gget ref -s homo_sapiens -w gtf, dna
```
&rarr; Returns only the links to the latest human reference GTF and DNA files, in this order, in a space-separated list (terminal), or comma-separated list (Jupyter Lab / Google Colab).    
For Jupyter Lab / Google Colab: Combining this command with `save=True`, will save the results in a text file named ref_results.txt in the current working directory.

## gget search 
> :warning: **gget search currently only supports genes listed in the Ensembl core API, which includes limited external references.** Manually searching the [Ensembl website](https://uswest.ensembl.org/) might yield more results.

### Options
`-sw` `--searchwords`  
One or more free form searchwords for the query, e.g. gaba, nmda. Searchwords are not case-sensitive.

`-s` `--species`  
Species or database to be searched.  
Species can be passed in the format 'genus_species', e.g. 'homo_sapiens'. To pass a specific CORE database (e.g. a specific mouse strain), enter the name of the CORE database, e.g. 'mus_musculus_dba2j_core_105_1'. All availabale species databases can be found here: http://ftp.ensembl.org/pub/release-105/mysql/

`-l` `--limit`  
Limits the number of search results to the top `[limit]` genes found.

`-o` `--out`  
Path to the file the results will be saved in, e.g. path/to/directory/results.csv (default: None &rarr; just prints results).  
For Jupyter Lab / Google Colab: `save=True` will save the output to the current working directory.

### gget search EXAMPLES
#### Query Ensembl for genes from a specific species using multiple searchwords
Jupyter Lab / Google Colab:
```python
search(["gaba", "gamma-aminobutyric acid"], "homo_sapiens")
```

Terminal:
```
$ gget search -sw gaba, gamma-aminobutyric acid -s homo_sapiens
```
&rarr; 

### Query Ensembl for genes from a specific species using a single searchword
Jupyter Lab:
```python
gget("searchword1", "genus_species")
```

Terminal:
```
gget search -sw searchword1 -sp genus_species
```

### Query Ensembl for genes from a specific species using multiple searchwords while limiting the number of returned search results 
For example, limiting the number of results to 10:  
Jupyter Lab:
```python
gget(["searchword1", searchword2, searchword3"], "genus_species", limit=10)
```

Terminal:
```
gget search -sw searchword1 searchword2 searchword3 -sp genus_species -l 10
```

### Query Ensembl for genes from any of the 236 species databases found [here](http://ftp.ensembl.org/pub/release-105/mysql/).   
For example, for the database "nothobranchius_furzeri_core_105_2":  
Jupyter Lab:
```python
gget("searchword1", "nothobranchius_furzeri_core_105_2")
```

Terminal:
```
gget search -sw searchword1 -sp nothobranchius_furzeri_core_105_2 
```

**Note:**  
`gget search` supports the following species abbreviations:  
"homo_sapiens" -> "human"  
"mus_musculus" -> "mouse"  
"caenorhabditis_elegans" -> "roundworm"  
All other species can be called using their specific [database](http://ftp.ensembl.org/pub/release-105/mysql/), as shown in the example above.

### [Click here for examples](https://github.com/lauraluebbert/gget/tree/main/examples)

Author: Laura Luebbert  
[![DOI](https://zenodo.org/badge/458943224.svg)](https://zenodo.org/badge/latestdoi/458943224)
