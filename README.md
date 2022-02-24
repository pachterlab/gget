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
Look up gene or transcript Ensembl IDs for their common name, description, sequence, homologs, synonyms, corresponding transcript/gene and more.

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

#### [Click here for more examples](https://github.com/lauraluebbert/gget/tree/main/examples)

___
## gget search 
Query [Ensembl](https://www.ensembl.org/) for genes from a defined species using free form search words.  
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
&rarr; Returns all genes that contain at least one of the searchwords in their Ensembl or external reference description, in the format:

| Ensembl_ID     | Ensembl_description     | Ext_ref_description     | Biotype        | Gene_name | URL |
| -------------- |-------------------------| ------------------------| -------------- | ----------|-----|
| ENSG00000034713| GABA type A receptor associated protein like 2 [Source:HGNC Symbol;Acc:HGNC:13291] | 	GABA type A receptor associated protein like 2 | protein_coding | GABARAPL2 | https://uswest.ensembl.org/homo_sapiens/Gene/Summary?g=ENSG00000034713 |
| . . .            | . . .                     | . . .                     | . . .            | . . .       | . . . |


#### Query Ensembl for genes from a specific species using a single searchword
Jupyter Lab:
```python
search("gaba", "homo_sapiens")
```

Terminal:
```
$ gget search -sw gaba -s homo_sapiens
```
&rarr; Returns all genes that contain the searchword in their Ensembl or external reference description.

#### Query Ensembl for genes from a specific species using multiple searchwords while limiting the number of returned search results 
Jupyter Lab:
```python
search(["gaba", "gamma-aminobutyric acid"], "homo_sapiens", limit=10)
```

Terminal:
```
$ gget search -sw gaba, gamma-aminobutyric acid -s homo_sapiens -l 10
```
&rarr; Returns the first 10 genes that contain at least one of the searchwords in their Ensembl or external reference description.

#### Query Ensembl for genes from any of the 236 species databases found [here](http://ftp.ensembl.org/pub/release-105/mysql/), e.g. a specific mouse strain.   
Jupyter Lab:
```python
search("brain", "mus_musculus_cbaj_core_105_1")
```

Terminal:
```
$ gget search -sw brain -s mus_musculus_cbaj_core_105_1 
```
&rarr; Returns genes from the CBA/J mouse strain that contain the searchword in their Ensembl or external reference description.

#### [Click here for more examples](https://github.com/lauraluebbert/gget/tree/main/examples)

___
## gget lookup  
Look up gene or transcript Ensembl IDs for their common name, description, sequence, homologs, synonyms, corresponding transcript/gene and more from the Ensembl database as well as external references.

### Options
`-id` `--ens_ids`  
One or more Ensembl IDs.

`-seq` `--seq`  
Returns basepair sequence of gene (or parent gene if transcript ID is passed) (default: False).

`-H` `--homology`  
Returns homology information of ID (default: False).

`-x` `--xref`  
Returns information from external references (default: False).

`-o` `--out`  
Path to the file the results will be saved in, e.g. path/to/directory/results.json (default: None &rarr; just prints results).  
For Jupyter Lab / Google Colab: `save=True` will save the output to the current working directory.

### gget lookup EXAMPLES
#### Look up a list of Ensembl IDs
Jupyter Lab / Google Colab:
```python
lookup(["ENSG00000034713, ENSG00000104853, ENSG00000170296"])
```

Terminal:
```
$ gget lookup -id ENSG00000034713, ENSG00000104853, ENSG00000170296
```
&rarr; Returns a json containing information about each ID, amongst others the common name, description, and corresponding transcript/gene, in the format:

#ADD LOOKUP OUTPUT!!!

#### Look up one Ensembl ID and include sequence, homology information and external reference description
Jupyter Lab / Google Colab:
```python
lookup("ENSG00000034713", seq=True, homology=True, xref=True)
```
Terminal:
```
$ gget lookup -id ENSG00000034713 -seq -H -x
```
&rarr; Returns a json containing the sequence, homology information, and external reference description of each ID in addition to the standard information mentioned above.

#### [Click here for more examples](https://github.com/lauraluebbert/gget/tree/main/examples)

___ 
Author: Laura Luebbert  
