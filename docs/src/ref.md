Python arguments are equivalent to long-option arguments (`--arg`), unless otherwise specified. Flags are True/False arguments in Python.  
The manual for any gget tool can be called from the command-line using the `-h` `--help` flag.  

## gget ref 📖
Fetch FTPs and their respective metadata (or use flag `ftp` to only return the links) for reference genomes and annotations from [Ensembl](https://www.ensembl.org/) by species.  
Return format: dictionary/JSON.

**Positional argument**  
`species`  
Species for which the FTPs will be fetched in the format genus_species, e.g. homo_sapiens.  
Note: Not required when calling flag [--list_species].   
Supported shortcuts: 'human', 'mouse'

**Optional arguments**  
`-w` `--which`  
Defines which results to return. Default: 'all' -> Returns all available results.  
Possible entries are one or a combination (as comma-separated list) of the following:  
'gtf' - Returns the annotation (GTF).  
'cdna' - Returns the trancriptome (cDNA).  
'dna' - Returns the genome (DNA).  
'cds' - Returns the coding sequences corresponding to Ensembl genes. (Does not contain UTR or intronic sequence.)  
'cdrna' - Returns transcript sequences corresponding to non-coding RNA genes (ncRNA).  
'pep' - Returns the protein translations of Ensembl genes.  

`-r` `--release`  
Defines the Ensembl release number from which the files are fetched, e.g. 104. Default: latest Ensembl release.  

`-o` `--out`    
Path to the JSON file the results will be saved in, e.g. path/to/directory/results.json. Default: Standard out.  
Python: `save=True` will save the output in the current working directory.

**Flags**  
`-l` `--list_species`   
Lists all available species. (Python: combine with `species=None`.)  

`-ftp` `--ftp`   
Returns only the requested FTP links.  

`-d` `--download`   
Command-line only. Downloads the requested FTPs to the current directory (requires [curl](https://curl.se/docs/) to be installed).

  
### Examples
**Use `gget ref` in combination with [kallisto | bustools](https://www.kallistobus.tools/kb_usage/kb_ref/) to build a reference index:**
```bash
kb ref -i INDEX -g T2G -f1 FASTA $(gget ref --ftp -w dna,gtf homo_sapiens)
```
&rarr; kb ref builds a reference index using the latest DNA and GTF files of species **Homo sapiens** passed to it by `gget ref`.  
  
List all available genomes from Ensembl release 103:  
```bash
gget ref --list_species -r 103
```
```python
# Python
gget.ref(species=None, list_species=True, release=103)
```
&rarr; Returns a list with all available genomes (checks if GTF and FASTAs are available) from Ensembl release 103.   
(If no release is specified, `gget ref` will always return information from the latest Ensembl release.)  
  
Get the genome reference for a specific species:   
```bash
gget ref -w gtf,dna homo_sapiens
```
```python
# Python
gget.ref("homo_sapiens", which=["gtf", "dna"])
```
&rarr; Returns a JSON with the latest human GTF and FASTA FTPs, and their respective metadata, in the format:
```
{
    "homo_sapiens": {
        "annotation_gtf": {
            "ftp": "http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz",
            "ensembl_release": 106,
            "release_date": "28-Feb-2022",
            "release_time": "23:27",
            "bytes": "51379459"
        },
        "genome_dna": {
            "ftp": "http://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
            "ensembl_release": 106,
            "release_date": "21-Feb-2022",
            "release_time": "09:35",
            "bytes": "881211416"
        }
    }
}
```

#### [More examples](https://github.com/pachterlab/gget_examples)
