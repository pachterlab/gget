# gget

Query ENSEMBL for genes using free form search words (gget **search**) or fetch FTP download links by species (gget **FetchTP**).

## Installation
```
pip install gget
```

For use in Jupyter Lab:
```python
from gget import gget, fetchtp
```

## gget FetchTP usage

### Fetch GTF, DNA, and cDNA FTP links for a specific species

Jupyter Lab:
```python
fetchtp("genus_species")
```

Terminal:
```
gget fetchtp -sp genus_species
```
where `genus_species` defines the species for which the FTPs are fetched, e.g. `homo_sapiens`.

This returns a json with the GTF, DNA, and cDNA links, their respective release dates and time, and the Ensembl release from which the links were fetched in the format:
```python
{
            species: {
                "transcriptome": {
                    "ftp": cDNA FTP download URL,
                    "ensembl_release": Ensembl release #,
                    "release_date": Day-Month-Year,
                    "release_time": HH:MM,
                    "bytes": cDNA FTP file size in bytes
                },
                "genome": {
                    "ftp": DNA FTP download URL,
                    "ensembl_release": Ensembl release #,
                    "release_date": Day-Month-Year,
                    "release_time": HH:MM,
                    "bytes": DNA FTP file size in bytes
                },
                "annotation": {
                    "ftp": GTF FTP download URL,
                    "ensembl_release": Ensembl release #,
                    "release_date": Day-Month-Year,
                    "release_time": HH:MM,
                    "bytes": GTF FTP file size in bytes
                }
            }
        }
```

### Fetch GTF, DNA, and cDNA FTP links for a specific species from a specific Ensembl release, e.g. release 104

Jupyter Lab:
```python
fetchtp("genus_species", release=104)
```

Terminal:
```
gget fetchtp -sp genus_species -r 104
```
where `release` / `-r` defines the Ensembl release from which the FTPs are fetched (by default the latest release is used).

### Fetch only the GTF link for a specific species 

Jupyter Lab:
```python
fetchtp("genus_species", return_val="gtf")
```

Terminal:
```
gget fetchtp -sp genus_species -rv gtf
```
where `gtf` alters the return value from the default `json` such that only the annotation (GTF) download link for the defined species is returned. Alternative entries for `return_val` / `-rv` are `dna` or `cdna`, which return only the genome (DNA) or the transcriptome (cDNA) download links, respectively.

This functionality can be combined with single-cell RNA-seq data pre-processing tools such as [kallisto bustools](https://pachterlab.github.io/kallistobustools/kb_usage/kb_ref/) or [cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references) to build a transcriptome index by automatically fetching the latest FTP links from Ensembl:
```
# kb ref
kb ref \
-i INDEX \
-g T2G \
-f1 FASTA \
$(gget fetchtp -sp homo_sapiens -rv dna) \
$(gget fetchtp -sp homo_sapiens -rv gtf)

# cellranger mkref
cellranger mkref \
--genome=output_genome \
--fasta=$(gget fetchtp -sp genus_species -rv dna)
--genes=$(gget fetchtp -sp genus_species -rv gtf)
```

## gget search usage

### Query Ensembl for genes from a specific species using multiple searchwords
Jupyter Lab:
```python
gget(["searchword1", "searchword2", "searchword3"], "genus_species")
```

Terminal:
```
gget search -w searchword1 searchword2 searchword3 -sp genus_species
```

### Query Ensembl for genes from a specific species using a single searchword
Jupyter Lab:
```python
gget("searchword1", "genus_species")
```

Terminal:
```
gget search -w searchword1 -sp genus_species
```

### Query Ensembl for genes from a specific species using multiple searchwords, limiting the number of results to 10
Jupyter Lab:
```python
gget(["searchword1", searchword2, searchword3"], "genus_species", limit=10)
```

Terminal:
```
gget search -w searchword1 searchword2 searchword3 -sp genus_species -l 10
```

### Query Ensembl for genes from any of the 236 species databases found [here](http://ftp.ensembl.org/pub/release-105/mysql/). For example, for the database "nothobranchius_furzeri_core_105_2":
Jupyter Lab:
```python
gget("searchword1", "nothobranchius_furzeri_core_105_2")
```

Terminal:
```
gget search -w searchword1 -sp nothobranchius_furzeri_core_105_2 
```

