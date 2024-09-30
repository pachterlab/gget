> Python arguments are equivalent to long-option arguments (`--arg`), unless otherwise specified. Flags are True/False arguments in Python. The manual for any gget tool can be called from the command-line using the `-h` `--help` flag.  
## gget ref 📖
Fetch FTPs and their respective metadata (or use flag `ftp` to only return the links) for reference genomes and annotations from [Ensembl](https://www.ensembl.org/) by species.  
Return format: dictionary/JSON.

**Positional argument**  
`species`  
Species for which the FTPs will be fetched in the format genus_species, e.g. homo_sapiens.  
Supports all available vertebrate and invertebrate (plants, fungi, protists, and invertebrate metazoa) genomes from Ensembl, except bacteria.  
Note: Not required when using flags `--list_species` or `--list_iv_species`.   
Supported shortcuts: 'human', 'mouse', 'human_grch37' (accesses the GRCh37 genome assembly)

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

`-od` `--out_dir`   
Path to the directory where the FTPs will be saved, e.g. path/to/directory/. Default: Current working directory.

`-o` `--out`    
Path to the JSON file the results will be saved in, e.g. path/to/directory/results.json. Default: Standard out.  
Python: `save=True` will save the output in the current working directory.

**Flags**  
`-l` `--list_species`   
Lists all available vertebrate species. (Python: combine with `species=None`.)  

`-liv` `--list_iv_species`   
Lists all available invertebrate species. (Python: combine with `species=None`.)  

`-ftp` `--ftp`   
Returns only the requested FTP links.  

`-d` `--download`   
Command-line only. Downloads the requested FTPs to the directory specified by `out_dir` (requires [curl](https://curl.se/docs/) to be installed).

`-q` `--quiet`   
Command-line only. Prevents progress information from being displayed.  
Python: Use `verbose=False` to prevent progress information from being displayed. 
  
  
### Examples
**Use `gget ref` in combination with [kallisto | bustools](https://www.kallistobus.tools/kb_usage/kb_ref/) to build a reference index:**
```bash
kb ref -i INDEX -g T2G -f1 FASTA $(gget ref --ftp -w dna,gtf homo_sapiens)
```
&rarr; kb ref builds a reference index using the latest DNA and GTF files of species **Homo sapiens** passed to it by `gget ref`.  

<br/><br/>

**List all available genomes from Ensembl release 103:**  
```bash
gget ref --list_species -r 103
```
```python
# Python
gget.ref(species=None, list_species=True, release=103)
```
&rarr; Returns a list with all available genomes (checks if GTF and FASTAs are available) from Ensembl release 103.   
(If no release is specified, `gget ref` will always return information from the latest Ensembl release.)  

<br/><br/>

**Get the genome reference for a specific species:**   
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

# References
If you use `gget ref` in a publication, please cite the following articles:   

- Luebbert, L., & Pachter, L. (2023). Efficient querying of genomic reference databases with gget. Bioinformatics. [https://doi.org/10.1093/bioinformatics/btac836](https://doi.org/10.1093/bioinformatics/btac836)

- Martin FJ, Amode MR, Aneja A, Austine-Orimoloye O, Azov AG, Barnes I, Becker A, Bennett R, Berry A, Bhai J, Bhurji SK, Bignell A, Boddu S, Branco Lins PR, Brooks L, Ramaraju SB, Charkhchi M, Cockburn A, Da Rin Fiorretto L, Davidson C, Dodiya K, Donaldson S, El Houdaigui B, El Naboulsi T, Fatima R, Giron CG, Genez T, Ghattaoraya GS, Martinez JG, Guijarro C, Hardy M, Hollis Z, Hourlier T, Hunt T, Kay M, Kaykala V, Le T, Lemos D, Marques-Coelho D, Marugán JC, Merino GA, Mirabueno LP, Mushtaq A, Hossain SN, Ogeh DN, Sakthivel MP, Parker A, Perry M, Piližota I, Prosovetskaia I, Pérez-Silva JG, Salam AIA, Saraiva-Agostinho N, Schuilenburg H, Sheppard D, Sinha S, Sipos B, Stark W, Steed E, Sukumaran R, Sumathipala D, Suner MM, Surapaneni L, Sutinen K, Szpak M, Tricomi FF, Urbina-Gómez D, Veidenberg A, Walsh TA, Walts B, Wass E, Willhoft N, Allen J, Alvarez-Jarreta J, Chakiachvili M, Flint B, Giorgetti S, Haggerty L, Ilsley GR, Loveland JE, Moore B, Mudge JM, Tate J, Thybert D, Trevanion SJ, Winterbottom A, Frankish A, Hunt SE, Ruffier M, Cunningham F, Dyer S, Finn RD, Howe KL, Harrison PW, Yates AD, Flicek P. Ensembl 2023. Nucleic Acids Res. 2023 Jan 6;51(D1):D933-D941. doi: [10.1093/nar/gkac958](https://doi.org/10.1093/nar/gkac958). PMID: 36318249; PMCID: PMC9825606.
