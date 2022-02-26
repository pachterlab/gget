# gget (gene-get)
[![pypi version](https://img.shields.io/pypi/v/gget)](https://pypi.org/project/gget/0.0.6/)
![PyPI - Downloads](https://img.shields.io/pypi/dm/gget)
[![license](https://img.shields.io/pypi/l/gget)](LICENSE)
[![DOI](https://zenodo.org/badge/458943224.svg)](https://zenodo.org/badge/latestdoi/458943224)

gget has three main commands:  
- [**gget ref**](#gget-ref)  
Fetch links to GTF and FASTA files from the [Ensembl FTP site](http://ftp.ensembl.org/pub/).
- [**gget search**](#gget-search)   
Query [Ensembl](https://www.ensembl.org/) for genes using free form search words.
- [**gget spy**](#gget-spy)  
Look up genes or transcripts by their Ensembl ID.

## Installation
```bash
pip install gget
```

For use in Jupyter Lab / Google Colab:
```python
from gget import ref, search, spy
```

## Getting started
```bash
# Fetch Homo sapiens  GTF, DNA, and cDNA FTPs from Ensembl release 104
$ gget ref -s homo_sapiens -r 104

# Search zebra finch genes with "mito" in their description and limit to the top 10 genes
$ gget search -sw mito -s taeniopygia_guttata -l 10

# Look up Ensembl ID ENSSCUG00000017183 and return its sequence
$ gget spy -id ENSSCUG00000017183 -seq
```

# Manual
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

### gget ref examples
#### Fetch GTF, DNA, and cDNA FTP links for a specific species

```python
# Jupyter Lab / Google Colab:
ref("homo_sapiens")

# Terminal:
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

```python
# Jupyter Lab / Google Colab:
ref("homo_sapiens", release=104)

# Terminal
$ gget ref -s homo_sapiens -r 104
```
&rarr; Returns a json with the human reference genome GTF, DNA, and cDNA links, and their respective release dates and time, from Ensembl release 104.

#### Save the results

```python
# Jupyter Lab / Google Colab:
ref("homo_sapiens", save=True)

# Terminal 
$ gget ref -s homo_sapiens -o path/to/directory/ref_results.json
```
&rarr; Saves the results in path/to/directory/ref_results.json.  
For Jupyter Lab / Google Colab: Saves the results in a json file named ref_results.json in the current working directory.

#### Fetch only certain types of links for a specific species 

```python
# Jupyter Lab / Google Colab:
ref("homo_sapiens", which=["gtf", "dna"])

# Terminal 
$ gget ref -s homo_sapiens -w gtf,dna
```
&rarr; Returns a dictionary/json containing the latest human reference GTF and DNA files, in this order, and their respective release dates and time.    

#### Fetch only certain types of links for a specific species and return only the links

```python
# Jupyter Lab / Google Colab:
ref("homo_sapiens", which=["gtf", "dna"], ftp=True)

# Terminal 
$ gget ref -s homo_sapiens -w gtf,dna -ftp
```
&rarr; Returns only the links (wihtout additional information) to the latest human reference GTF and DNA files, in this order, in a space-separated list (terminal), or comma-separated list (Jupyter Lab / Google Colab).    
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

```python
# Jupyter Lab / Google Colab:
search(["gaba", "gamma-aminobutyric"], "homo_sapiens")

# Terminal 
$ gget search -sw gaba,gamma-aminobutyric -s homo_sapiens
```
&rarr; Returns all genes that contain at least one of the searchwords in their Ensembl or external reference description, in the format:

| Ensembl_ID     | Ensembl_description     | Ext_ref_description     | Biotype        | Gene_name | URL |
| -------------- |-------------------------| ------------------------| -------------- | ----------|-----|
| ENSG00000034713| GABA type A receptor associated protein like 2 [Source:HGNC Symbol;Acc:HGNC:13291] | 	GABA type A receptor associated protein like 2 | protein_coding | GABARAPL2 | https://uswest.ensembl.org/homo_sapiens/Gene/Summary?g=ENSG00000034713 |
| . . .            | . . .                     | . . .                     | . . .            | . . .       | . . . |


#### Query Ensembl for genes from a specific species using a single searchword

```python
# Jupyter Lab / Google Colab:
search("gaba", "homo_sapiens")

# Terminal 
$ gget search -sw gaba -s homo_sapiens
```
&rarr; Returns all genes that contain the searchword in their Ensembl or external reference description.

#### Query Ensembl for genes from a specific species using multiple searchwords while limiting the number of returned search results 

```python
# Jupyter Lab / Google Colab:
search("gaba", "homo_sapiens", limit=10)

# Terminal 
$ gget search -sw gaba -s homo_sapiens -l 10
```
&rarr; Returns the first 10 genes that contain the searchword in their Ensembl or external reference description. If more than one searchword is passed, `limit` will limit the number of genes per searchword. 

#### Query Ensembl for genes from any of the 236 species databases found [here](http://ftp.ensembl.org/pub/release-105/mysql/), e.g. a specific mouse strain.   

```python
# Jupyter Lab / Google Colab:
search("brain", "mus_musculus_cbaj_core_105_1")

# Terminal 
$ gget search -sw brain -s mus_musculus_cbaj_core_105_1 
```
&rarr; Returns genes from the CBA/J mouse strain that contain the searchword in their Ensembl or external reference description.

#### [Click here for more examples](https://github.com/lauraluebbert/gget/tree/main/examples)

___
## gget spy  
"Spy" on gene or transcript Ensembl IDs for their common name, description, sequence, homologs, synonyms, corresponding transcript/gene and more from the Ensembl database as well as external references.

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

### gget spy EXAMPLES
#### Look up a list of Ensembl IDs

```python
# Jupyter Lab / Google Colab:
spy(["ENSG00000034713", "ENSG00000104853", "ENSG00000170296"])

# Terminal 
$ gget spy -id ENSG00000034713,ENSG00000104853,ENSG00000170296
```
&rarr; Returns a json containing information about each ID, amongst others the common name, description, and corresponding transcript/gene, in the format:
```
{
            "Ensembl ID": {
                        "species": genus_species,
                        "object_type": e.g. Gene,
                        "biotype": e.g. protein_coding,
                        "canonical_transcript": Transcript ID,
                        "display_name": Common gene name,
                        "description": Ensemble description,
                        "assembly_name": Name of species assmebly,
                        "seq_region_name": Sequence region,
                        "start": Sequence start position,
                        "end": Sequence end position,
                        "strand": Strand
                        },
}
```
#### Look up one Ensembl ID and include sequence, homology information and external reference description

```python
# Jupyter Lab / Google Colab:
spy("ENSG00000034713", seq=True, homology=True, xref=True)

# Terminal 
$ gget spy -id ENSG00000034713 -seq -H -x
```
&rarr; Returns a json containing the sequence, homology information, and external reference description of each ID in addition to the standard information mentioned above.

#### [Click here for more examples](https://github.com/lauraluebbert/gget/tree/main/examples)

___ 
Author: Laura Luebbert 
