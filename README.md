# gget
[![pypi version](https://img.shields.io/pypi/v/gget)](https://pypi.org/project/gget)
[![Downloads](https://pepy.tech/badge/gget)](https://pepy.tech/project/gget)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/gget/README.html)
[![license](https://img.shields.io/pypi/l/gget)](LICENSE)
![status](https://github.com/pachterlab/gget/workflows/CI/badge.svg)
![Code Coverage](https://img.shields.io/badge/Coverage-83%25-green.svg)  

## ✨ Version ≥ 0.3.0: [`gget alphafold`](#gget-alphafold-) 

## ✨ What's new in version ≥ 0.2.0
- JSON is now the default output format for the command-line interface for modules that previously returned data frame (CSV) format by default (the output can be converted to data frame/CSV using flag `[-csv][--csv]`). Data frame/CSV remains the default output for Jupyter Lab / Google Colab (and can be converted to JSON with `json=True`).
- For all modules, the first required argument was converted to a positional argument and should not be named anymore in the command-line, e.g. `gget ref -s human` &rarr; `gget ref human`.
- `gget info`: `[--expand]` is deprecated. The module will now always return all of the available information.
- Slight changes to the output returned by `gget info`, including the return of versioned Ensembl IDs.
- `gget info` and `gget seq` now support 🪱 WormBase and 🪰 FlyBase IDs.
- `gget archs4` and `gget enrichr` now also take Ensembl IDs as input with added flag `[-e][--ensembl]` (`ensembl=True` in Jupyter Lab / Google Colab).
- `gget seq` argument `seqtype` was replaced by flag `[-t][--translate]` (`translate=True/False` in Jupyter Lab / Google Colab) which will return either nucleotide (`False`) or amino acid (`True`) sequences.
- `gget search` argument `seqtype` was renamed to `id_type` for clarity (still taking the same arguments 'gene' or 'transcript').
- Version ≥ 0.2.6: `gget ref` supports plant genomes! 🌱  
  
Note: [UniProt](https://www.uniprot.org/) changed the structure of their API on June 28, 2022. Please upgrade to `gget` version ≥ 0.2.5 if you use any of the modules querying data from UniProt (`gget info` and `gget seq`).  
___

`gget` is a free, open-source command-line tool and Python package that enables efficient querying of genomic databases. `gget`  consists of a collection of separate but interoperable modules, each designed to facilitate one type of database querying in a single line of code.  

If you use `gget` in a publication, please [cite*](#cite):    
```
Luebbert, L. & Pachter, L. (2022). Efficient querying of genomic reference databases with gget. bioRxiv 2022.05.17.492392
```
Read the manuscript here: https://doi.org/10.1101/2022.05.17.492392  
  
![alt text](https://github.com/pachterlab/gget/blob/main/figures/gget_overview.png?raw=true)
  
`gget` currently consists of the following nine modules:
- [`gget ref`](#gget-ref-)  
Fetch File Transfer Protocols (FTPs) and metadata for reference genomes and annotations from [Ensembl](https://www.ensembl.org/) by species.
- [`gget search`](#gget-search-)   
Fetch genes and transcripts from [Ensembl](https://www.ensembl.org/) using free-form search terms.
- [`gget info`](#gget-info-)  
Fetch extensive gene and transcript metadata from [Ensembl](https://www.ensembl.org/), [UniProt](https://www.uniprot.org/), and [NCBI](https://www.ncbi.nlm.nih.gov/) using Ensembl IDs.  
- [`gget seq`](#gget-seq-)  
Fetch nucleotide or amino acid sequences of genes or transcripts from [Ensembl](https://www.ensembl.org/) or [UniProt](https://www.uniprot.org/), respectively.  
- [`gget blast`](#gget-blast-)  
BLAST a nucleotide or amino acid sequence to any [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) database.
- [`gget blat`](#gget-blat-)  
Find the genomic location of a nucleotide or amino acid sequence using [BLAT](https://genome.ucsc.edu/cgi-bin/hgBlat).
- [`gget muscle`](#gget-muscle-)  
Align multiple nucleotide or amino acid sequences to each other using [Muscle5](https://www.drive5.com/muscle/).
- [`gget enrichr`](#gget-enrichr-)  
Perform an enrichment analysis on a list of genes using [Enrichr](https://maayanlab.cloud/Enrichr/).
- [`gget archs4`](#gget-archs4-)  
Find the most correlated genes to a gene of interest or find the gene's tissue expression atlas using [ARCHS4](https://maayanlab.cloud/archs4/).
- [`gget alphafold`](#gget-alphafold-)  
Predict the 3D structure of a protein from its amino acid sequence using a simplified version of [DeepMind](https://www.deepmind.com/)’s [AlphaFold2](https://github.com/deepmind/alphafold).  


## Installation
```bash
pip install --upgrade gget
```
or
```bash
conda install -c bioconda gget
```

For use in Jupyter Lab / Google Colab:
```python
import gget
```

## 🪄 Quick start guide
```bash
# Fetch all Homo sapiens reference and annotation FTPs from the latest Ensembl release
$ gget ref homo_sapiens

# Search human genes with "ace2" or "angiotensin converting enzyme 2" in their name/description
$ gget search -s homo_sapiens 'ace2' 'angiotensin converting enzyme 2'

# Look up gene ENSG00000130234 (ACE2) and its transcript ENST00000252519
$ gget info ENSG00000130234 ENST00000252519

# Fetch the amino acid sequence of the canonical transcript of gene ENSG00000130234
$ gget seq --translate ENSG00000130234

# Quickly find the genomic location of (the start of) that amino acid sequence
$ gget blat MSSSSWLLLSLVAVTAAQSTIEEQAKTFLDKFNHEAEDLFYQSSLAS

# Blast (the start of) that amino acid sequence
$ gget blast MSSSSWLLLSLVAVTAAQSTIEEQAKTFLDKFNHEAEDLFYQSSLAS

# Align nucleotide or amino acid sequences stored in a FASTA file
$ gget muscle path/to/file.fa

# Use Enrichr for an ontology analysis of a list of genes
$ gget enrichr -db ontology ACE2 AGT AGTR1 ACE AGTRAP AGTR2 ACE3P

# Get the human tissue expression of gene ACE2
$ gget archs4 -w tissue ACE2

# Predict the protein structure of GFP from its amino acid sequence
$ gget setup alphafold # setup only needs to be run once
$ gget alphafold MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK

```
Jupyter Lab / Google Colab:
```python  
import gget
gget.ref("homo_sapiens")
gget.search(["ace2", "angiotensin converting enzyme 2"], "homo_sapiens")
gget.info(["ENSG00000130234", "ENST00000252519"])
gget.seq("ENSG00000130234", translate=True)
gget.blat("MSSSSWLLLSLVAVTAAQSTIEEQAKTFLDKFNHEAEDLFYQSSLAS")
gget.blast("MSSSSWLLLSLVAVTAAQSTIEEQAKTFLDKFNHEAEDLFYQSSLAS")
gget.muscle("path/to/file.fa")
gget.enrichr(["ACE2", "AGT", "AGTR1", "ACE", "AGTRAP", "AGTR2", "ACE3P"], database="ontology", plot=True)
gget.archs4("ACE2", which="tissue")

gget.setup("alphafold") # setup only needs to be run once
gget.alphafold("MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK")
```
#### [More examples](https://github.com/pachterlab/gget_examples)
___

# Manual
Jupyter Lab / Google Colab arguments are equivalent to long-option arguments (`--arg`), unless otherwise specified. Flags are True/False arguments in Jupyter Lab / Google Colab.  
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
Jupyter Lab / Google Colab: `save=True` will save the output in the current working directory.

**Flags**  
`-l` `--list_species`   
Lists all available species. (Jupyter Lab / Google Colab: combine with `species=None`.)  

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
# Jupyter Lab / Google Colab:
gget.ref(species=None, list_species=True, release=103)
```
&rarr; Returns a list with all available genomes (checks if GTF and FASTAs are available) from Ensembl release 103.   
(If no release is specified, `gget ref` will always return information from the latest Ensembl release.)  
  
Get the genome reference for a specific species:   
```bash
gget ref -w gtf,dna homo_sapiens
```
```python
# Jupyter Lab / Google Colab:
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
___

## gget search 🔎
Fetch genes and transcripts from [Ensembl](https://www.ensembl.org/) using free-form search terms.   
Return format: JSON (command-line) or data frame/CSV (Jupyter Lab / Google Colab).

**Positional argument**
`searchwords`   
One or more free form search words, e.g. gaba nmda. (Note: Search is not case-sensitive.)

**Other required arguments**   
`-s` `--species`  
Species or database to be searched.  
A species can be passed in the format 'genus_species', e.g. 'homo_sapiens'.  
To pass a specific database, pass the name of the CORE database, e.g. 'mus_musculus_dba2j_core_105_1'.  
All availabale databases can be found [here](http://ftp.ensembl.org/pub/release-106/mysql/).  
Supported shortcuts: 'human', 'mouse'. 

**Optional arguments**  
`-t` `--id_type`  
'gene' (default) or 'transcript'  
Returns genes or transcripts, respectively.

`-ao` `--andor`  
'or' (default) or 'and'  
'or': Returns all genes that INCLUDE AT LEAST ONE of the searchwords in their name/description.  
'and': Returns only genes that INCLUDE ALL of the searchwords in their name/description.

`-l` `--limit`   
Limits the number of search results, e.g. 10. Default: None.  

`-o` `--out`  
Path to the csv the results will be saved in, e.g. path/to/directory/results.csv (or .json). Default: Standard out.   
Jupyter Lab / Google Colab: `save=True` will save the output in the current working directory.

**Flags**  
`-csv` `--csv`  
Command-line only. Returns results in CSV format.  
Jupyter Lab / Google Colab: Use `json=True` to return output in JSON format.

`wrap_text`  
Jupyter Lab / Google Colab only. `wrap_text=True` displays data frame with wrapped text for easy reading (default: False).  
  
    
### Example
```bash
gget search -s human gaba gamma-aminobutyric
```
```python
# Jupyter Lab / Google Colab:
gget.search(["gaba", "gamma-aminobutyric"], "homo_sapiens")
```
&rarr; Returns all genes that contain at least one of the search words in their name or Ensembl/external reference description:

| ensembl_id     | gene_name     | ensembl_description     | ext_ref_description        | biotype | url |
| -------------- |-------------------------| ------------------------| -------------- | ----------|-----|
| ENSG00000034713| GABARAPL2 | 	GABA type A receptor associated protein like 2 [Source:HGNC Symbol;Acc:HGNC:13291] | GABA type A receptor associated protein like 2 | protein_coding | https://uswest.ensembl.org/homo_sapiens/Gene/Summary?g=ENSG00000034713 |
| . . .            | . . .                     | . . .                     | . . .            | . . .       | . . . |
    
#### [More examples](https://github.com/pachterlab/gget_examples)
___

## gget info 💡
Fetch extensive gene and transcript metadata from [Ensembl](https://www.ensembl.org/), [UniProt](https://www.uniprot.org/), and [NCBI](https://www.ncbi.nlm.nih.gov/) using Ensembl IDs.  
Return format: JSON (command-line) or data frame/CSV (Jupyter Lab / Google Colab).

**Positional argument**  
`ens_ids`   
One or more Ensembl IDs.

**Optional arguments**  
`-o` `--out`   
Path to the file the results will be saved in, e.g. path/to/directory/results.csv (or .json). Default: Standard out.    
Jupyter Lab / Google Colab: `save=True` will save the output in the current working directory.

**Flags**  
`-q` `--quiet`   
Command-line only. Prevents progress information from being displayed.  
Jupyter Lab / Google Colab: Use `verbose=False` to prevent progress information from being displayed.  

`-csv` `--csv`  
Command-line only. Returns results in CSV format.  
Jupyter Lab / Google Colab: Use `json=True` to return output in JSON format.

`wrap_text`  
Jupyter Lab / Google Colab only. `wrap_text=True` displays data frame with wrapped text for easy reading (default: False).  

  
### Example
```bash
gget info ENSG00000034713 ENSG00000104853 ENSG00000170296
```
```python
# Jupyter Lab / Google Colab:
gget.info(["ENSG00000034713", "ENSG00000104853", "ENSG00000170296"])
```
&rarr; Returns extensive information about each requested Ensembl ID:  

|      | uniprot_id     | ncbi_gene_id     | primary_gene_name | synonyms | protein_names | ensembl_description | uniprot_description | ncbi_description | biotype | canonical_transcript | ... |
| -------------- |-------------------------| ------------------------| -------------- | ----------|-----|----|----|----|----|----|----|
| ENSG00000034713| P60520 | 11345 | GABARAPL2 | [ATG8, ATG8C, FLC3A, GABARAPL2, GATE-16, GATE16, GEF-2, GEF2] | Gamma-aminobutyric acid receptor-associated protein like 2 (GABA(A) receptor-associated protein-like 2)... | GABA type A receptor associated protein like 2 [Source:HGNC Symbol;Acc:HGNC:13291] | FUNCTION: Ubiquitin-like modifier involved in intra- Golgi traffic (By similarity). Modulates intra-Golgi transport through coupling between NSF activity and ... | Enables ubiquitin protein ligase binding activity. Involved in negative regulation of proteasomal protein catabolic process and protein... | protein_coding | ENST00000037243.7 |... |
| . . .            | . . .                     | . . .                     | . . .            | . . .       | . . . | . . . | . . . | . . . | . . . | . . . | ... |
  
#### [More examples](https://github.com/pachterlab/gget_examples)
___

## gget seq 🧬
Fetch nucleotide or amino acid sequence of a gene (and all its isoforms) or a transcript by Ensembl ID.   
Return format: FASTA.

**Positional argument**  
`ens_ids`   
One or more Ensembl IDs.

**Optional arguments**  
`-o` `--out`   
Path to the file the results will be saved in, e.g. path/to/directory/results.fa. Default: Standard out.   
Jupyter Lab / Google Colab: `save=True` will save the output in the current working directory.

**Flags**  
`-t` `--translate`  
Returns amino acid (instead of nucleotide) sequences.  
Nucleotide sequences are fetched from [Ensembl](https://www.ensembl.org/).  
Amino acid sequences are fetched from [UniProt](https://www.uniprot.org/).

`-iso` `--isoforms`   
Returns the sequences of all known transcripts.  
(Only for gene IDs.)
   
  
### Examples  
```bash
gget seq ENSG00000034713 ENSG00000104853 ENSG00000170296
```
```python
# Jupyter Lab / Google Colab:
gget.seq(["ENSG00000034713", "ENSG00000104853", "ENSG00000170296"])
```
&rarr; Returns the nucleotide sequences of ENSG00000034713, ENSG00000104853, and ENSG00000170296 in FASTA format.  
  

```bash
gget seq -t -iso ENSG00000034713
```
```python
# Jupyter Lab / Google Colab:
gget.seq("ENSG00000034713", translate=True, isoforms=True)
```
&rarr; Returns the amino acid sequences of all known transcripts of ENSG00000034713 in FASTA format.

#### [More examples](https://github.com/pachterlab/gget_examples)
___

## gget blast 💥
BLAST a nucleotide or amino acid sequence to any [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) database.  
Return format: JSON (command-line) or data frame/CSV (Jupyter Lab / Google Colab).

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
Jupyter Lab / Google Colab: `save=True` will save the output in the current working directory.

**Flags**  
`-lcf` `--low_comp_filt`  
Turns on [low complexity filter](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=FAQ#LCR).  

`-mbo` `--megablast_off`  
Turns off MegaBLAST algorithm. Default: MegaBLAST on (blastn only).  

`-q` `--quiet`   
Command-line only. Prevents progress information from being displayed.  
Jupyter Lab / Google Colab: Use `verbose=False` to prevent progress information from being displayed.  

`-csv` `--csv`  
Command-line only. Returns results in CSV format.  
Jupyter Lab / Google Colab: Use `json=True` to return output in JSON format.

`wrap_text`  
Jupyter Lab / Google Colab only. `wrap_text=True` displays data frame with wrapped text for easy reading (default: False).   
  
### Example
```bash
gget blast MKWMFKEDHSLEHRCVESAKIRAKYPDRVPVIVEKVSGSQIVDIDKRKYLVPSDITVAQFMWIIRKRIQLPSEKAIFLFVDKTVPQSR
```
```python
# Jupyter Lab / Google Colab:
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
# Jupyter Lab / Google Colab:
gget.blast("fasta.fa")
```
&rarr; Returns the BLAST results of the first sequence contained in the fasta.fa file. 

#### [More examples](https://github.com/pachterlab/gget_examples)
___

## gget blat 🎯
Find the genomic location of a nucleotide or amino acid sequence using [BLAT](https://genome.ucsc.edu/cgi-bin/hgBlat).   
Return format: JSON (command-line) or data frame/CSV (Jupyter Lab / Google Colab).

**Positional argument**  
`sequence`   
Nucleotide or amino acid sequence, or path to FASTA or .txt file.

**Optional arguments**  
`-st` `--seqtype`    
'DNA', 'protein', 'translated%20RNA', or 'translated%20DNA'.   
Default: 'DNA' for nucleotide sequences; 'protein' for amino acid sequences.  

`-a` `--assembly`  
'human' (hg38) (default), 'mouse' (mm39), 'zebrafinch' (taeGut2),   
or any of the species assemblies available [here](https://genome.ucsc.edu/cgi-bin/hgBlat) (use short assembly name).

`-o` `--out`   
Path to the file the results will be saved in, e.g. path/to/directory/results.csv (or .json). Default: Standard out.   
Jupyter Lab / Google Colab: `save=True` will save the output in the current working directory.  
  
**Flags**  
`-csv` `--csv`  
Command-line only. Returns results in CSV format.  
Jupyter Lab / Google Colab: Use `json=True` to return output in JSON format.
  

### Example
```bash
gget blat -a taeGut2 MKWMFKEDHSLEHRCVESAKIRAKYPDRVPVIVEKVSGSQIVDIDKRKYLVPSDITVAQFMWIIRKRIQLPSEKAIFLFVDKTVPQSR
```
```python
# Jupyter Lab / Google Colab:
gget.blat("MKWMFKEDHSLEHRCVESAKIRAKYPDRVPVIVEKVSGSQIVDIDKRKYLVPSDITVAQFMWIIRKRIQLPSEKAIFLFVDKTVPQSR", assembly="taeGut2")
```
&rarr; Returns BLAT results for assembly taeGut2 (zebra finch). In the above example, `gget blat` automatically detects this sequence as an amino acid sequence and therefore sets the BLAT seqtype to *protein*.

| genome     | query_size     | aligned_start     | aligned_end        | matches | mismatches | %_aligned | ... |
| -------------- |-------------------------| ------------------------| -------------- | ----------|-----|---|---|
| taeGut2| 88 | 	12 | 88 | 77 | 0 | 87.5 | ... |

#### [More examples](https://github.com/pachterlab/gget_examples)
___

## gget muscle 🦾
Align multiple nucleotide or amino acid sequences to each other using [Muscle5](https://www.drive5.com/muscle/).  
Return format: ClustalW formatted standard out or aligned FASTA (.afa).  

**Positional argument**  
`fasta`   
Path to FASTA or .txt file containing the nucleotide or amino acid sequences to be aligned.  

**Optional arguments**  
`-o` `--out`   
Path to the aligned FASTA file the results will be saved in, e.g. path/to/directory/results.afa. Default: Standard out.   
Jupyter Lab / Google Colab: `save=True` will save the output in the current working directory.

**Flags**  
`-s5` `--super5`  
Aligns input using the [Super5 algorithm](https://drive5.com/muscle5/Muscle5_SuppMat.pdf) instead of the [Parallel Perturbed Probcons (PPP) algorithm](https://drive5.com/muscle5/Muscle5_SuppMat.pdf) to decrease time and memory.  
Use for large inputs (a few hundred sequences).
  
  
### Example
```bash
gget muscle fasta.fa
```
```python
# Jupyter Lab / Google Colab:
gget.muscle("fasta.fa")
```
&rarr; Returns an overview of the aligned sequences with ClustalW coloring. (To return an aligned FASTA (.afa) file, use `--out` argument (or `save=True` in Jupyter Lab/Google Colab).) In the above example, the 'fasta.fa' includes several sequences to be aligned (e.g. isoforms returned from `gget seq`). 

![alt text](https://github.com/pachterlab/gget/blob/main/figures/example_muscle_return.png?raw=true)

#### [More examples](https://github.com/pachterlab/gget_examples)
___

## gget enrichr 💰
Perform an enrichment analysis on a list of genes using [Enrichr](https://maayanlab.cloud/Enrichr/).  
Return format: JSON (command-line) or data frame/CSV (Jupyter Lab / Google Colab).
  
**Positional argument**  
`genes`  
Short names (gene symbols) of genes to perform enrichment analysis on, e.g. PHF14 RBM3 MSL1 PHF21A.  
Alternatively: use flag `--ensembl` to input a list of Ensembl gene IDs, e.g. ENSG00000106443 ENSG00000102317 ENSG00000188895.

**Other required arguments**  
`-db` `--database`  
Database to use as reference for the enrichment analysis.  
Supports any database listed [here](https://maayanlab.cloud/Enrichr/#libraries) under 'Gene-set Library' or one of the following shortcuts:  
'pathway'       (KEGG_2021_Human)  
'transcription'     (ChEA_2016)  
'ontology'      (GO_Biological_Process_2021)  
'diseases_drugs'   (GWAS_Catalog_2019)   
'celltypes'      (PanglaoDB_Augmented_2021)  
'kinase_interactions'   (KEA_2015)  
  
**Optional arguments**  
`-o` `--out`   
Path to the file the results will be saved in, e.g. path/to/directory/results.csv (or .json). Default: Standard out.   
Jupyter Lab / Google Colab: `save=True` will save the output in the current working directory.

`figsize`  
Jupyter Lab / Google Colab only. (width, height) of plot in inches. (Default: (10,10))

`ax`  
Jupyter Lab / Google Colab only. Pass a matplotlib axes object for plot customization. (Default: None)
  
**Flags**  
`-e` `--ensembl`  
Add this flag if `genes` are given as Ensembl gene IDs.  
 
`-csv` `--csv`  
Command-line only. Returns results in CSV format.  
Jupyter Lab / Google Colab: Use `json=True` to return output in JSON format.
  
`plot`  
Jupyter Lab / Google Colab only. `plot=True` provides a graphical overview of the first 15 results (default: False).  
  
  
### Example
```bash
gget enrichr -db ontology ACE2 AGT AGTR1
```
```python
# Jupyter Lab / Google Colab:
gget.enrichr(["ACE2", "AGT", "AGTR1"], database="ontology", plot=True)
```
&rarr; Returns pathways/functions involving genes ACE2, AGT, and AGTR1 from the *GO Biological Process 2021* database. In Jupyter Lab / Google Colab, `plot=True` returns a graphical overview of the results:

![alt text](https://github.com/pachterlab/gget/blob/main/figures/gget_enrichr_results.png?raw=true)

#### [More examples](https://github.com/pachterlab/gget_examples)
___

## gget archs4 🐁
Find the most correlated genes to a gene of interest or find the gene's tissue expression atlas using [ARCHS4](https://maayanlab.cloud/archs4/).  
Return format: JSON (command-line) or data frame/CSV (Jupyter Lab / Google Colab).

**Positional argument**  
`gene`  
Short name (gene symbol) of gene of interest, e.g. STAT4.  
Alternatively: use flag `--ensembl` to input an Ensembl gene IDs, e.g. ENSG00000138378.

**Optional arguments**  
 `-w` `--which`  
'correlation' (default) or 'tissue'.  
'correlation' returns a gene correlation table that contains the 100 most correlated genes to the gene of interest. The Pearson correlation is calculated over all samples and tissues in [ARCHS4](https://maayanlab.cloud/archs4/).  
'tissue' returns a tissue expression atlas calculated from human or mouse samples (as defined by 'species') in [ARCHS4](https://maayanlab.cloud/archs4/).  

`-s` `--species`  
'human' (default) or 'mouse'.   
Defines whether to use human or mouse samples from [ARCHS4](https://maayanlab.cloud/archs4/).  
(Only for tissue expression atlas.)

`-o` `--out`   
Path to the file the results will be saved in, e.g. path/to/directory/results.csv (or .json). Default: Standard out.   
Jupyter Lab / Google Colab: `save=True` will save the output in the current working directory.  
  
**Flags**   
`-e` `--ensembl`  
Add this flag if `gene` is given as an Ensembl gene ID.  

`-csv` `--csv`  
Command-line only. Returns results in CSV format.  
Jupyter Lab / Google Colab: Use `json=True` to return output in JSON format.
  
  
### Examples
```bash
gget archs4 ACE2
```
```python
# Jupyter Lab / Google Colab:
gget.archs4("ACE2")
```
&rarr; Returns the 100 most correlated genes to ACE2:  

| gene_symbol     | pearson_correlation     |
| -------------- |-------------------------| 
| SLC5A1 | 0.579634 | 	
| CYP2C18 | 0.576577 | 	
| . . . | . . . | 	

```bash
gget archs4 -w tissue ACE2
```
```python
# Jupyter Lab / Google Colab:
gget.archs4("ACE2", which="tissue")
```
&rarr; Returns the tissue expression of ACE2 (by default, human data is used):

| id     | min     | q1 |  median | q3 | max |
| ------ |--------| ------ |--------| ------ |--------| 
| System.Urogenital/Reproductive System.Kidney.RENAL CORTEX | 0.113644 | 8.274060 | 9.695840 | 10.51670 | 11.21970 |
| System.Digestive System.Intestine.INTESTINAL EPITHELIAL CELL | 0.113644 | 	5.905560 | 9.570450 | 13.26470 | 13.83590 | 
| . . . | . . . | . . . | . . . | . . . | . . . |

#### [More examples](https://github.com/pachterlab/gget_examples)

___

## gget setup ⚙️
Function to install/download third-party dependencies for a specified gget module.

**Positional argument**  
`module`  
gget module for which dependencies should be installed.  

### Example
```bash
gget setup alphafold
```
```python
# Jupyter Lab / Google Colab:
gget.setup("alphafold")
```
&rarr; Installs all (modified) third-party dependencies and downloads model parameters (~4GB) required to run `gget alphafold`.  

___

## gget alphafold 🪢
Predict the 3D structure of a protein from its amino acid sequence using a simplified version of [DeepMind](https://www.deepmind.com/)’s [AlphaFold2](https://github.com/deepmind/alphafold) originally released and benchmarked for [AlphaFold Colab](https://colab.research.google.com/github/deepmind/alphafold/blob/main/notebooks/AlphaFold.ipynb).  
Returns: Predicted structure (PDB) and alignment error (json).  

Before using `gget alphafold` for the first time, run `gget setup alphafold` / `gget.setup("alphafold")` once (also see `gget setup` above).  

**Positional argument**  
`sequence`  
Amino acid sequence (str), list of sequences (for multimers), or path to FASTA file.

**Optional arguments**  
`-o` `--out`   
Path to folder to save prediction results in (str). Default: "./[date_time]_gget_alphafold_prediction".  
  
**Flags**   
`-r` `--relax`   
AMBER relax the best model. 

`plot`  
Jupyter Lab / Google Colab only. `plot=True` provides an interactive, 3D graphical overview of the predicted structure and alignment quality using [py3Dmol](https://pypi.org/project/py3Dmol/) and [matplotlib](https://matplotlib.org/) (default: True).  

`show_sidechains`  
Jupyter Lab / Google Colab only. `show_sidechains=True` includes side chains in the plot (default: True).  
  
### Example
```bash
gget alphafold MAAHKGAEHHHKAAEHHEQAAKHHHAAAEHHEKGEHEQAAHHADTAYAHHKHAEEHAAQAAKHDAEHHAPKPH
```
```python
# Jupyter Lab / Google Colab:
gget.alphafold("MAAHKGAEHHHKAAEHHEQAAKHHHAAAEHHEKGEHEQAAHHADTAYAHHKHAEEHAAQAAKHDAEHHAPKPH")
```
&rarr; Returns the predicted structure (PDB) and predicted alignment error (.json) in a new folder ("./[date_time]_gget_alphafold_prediction"). The Python interface also returns the following plots:

https://user-images.githubusercontent.com/56094636/182939299-36ac2a8f-0560-4a64-b1f8-6cff39ef2a75.mp4

#### [More examples](https://github.com/pachterlab/gget_examples)
___

# Cite 

If you use `gget` in a publication, please cite:   
Luebbert L. & Pachter L. (2022). Efficient querying of genomic reference databases with _gget_. bioRxiv 2022.05.17.492392; doi: https://doi.org/10.1101/2022.05.17.492392

- If using `gget archs4`, please also cite:   
  Lachmann A, Torre D, Keenan AB, Jagodnik KM, Lee HJ, Wang L, Silverstein MC, Ma’ayan A. Massive mining of publicly available RNA-seq data from human and mouse. Nature Communications 9. Article number: 1366 (2018), doi:10.1038/s41467-018-03751-6

  Bray NL, Pimentel H, Melsted P and Pachter L, Near optimal probabilistic RNA-seq quantification, Nature Biotechnology 34, p 525--527 (2016). https://doi.org/10.1038/nbt.3519

- If using `gget enrichr`, please also cite:     
  Chen EY, Tan CM, Kou Y, Duan Q, Wang Z, Meirelles GV, Clark NR, Ma'ayan A.
Enrichr: interactive and collaborative HTML5 gene list enrichment analysis tool. BMC Bioinformatics. 2013; 128(14).  

  Kuleshov MV, Jones MR, Rouillard AD, Fernandez NF, Duan Q, Wang Z, Koplev S, Jenkins SL, Jagodnik KM, Lachmann A, McDermott MG, Monteiro CD, Gundersen GW, Ma'ayan A.
Enrichr: a comprehensive gene set enrichment analysis web server 2016 update. Nucleic Acids Research. 2016; gkw377.  

  Xie Z, Bailey A, Kuleshov MV, Clarke DJB., Evangelista JE, Jenkins SL, Lachmann A, Wojciechowicz ML, Kropiwnicki E, Jagodnik KM, Jeon M, & Ma’ayan A.
Gene set knowledge discovery with Enrichr. Current Protocols, 1, e90. 2021. doi: 10.1002/cpz1.90. 

- If using `gget blast`, please also cite:  
  Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. Basic local alignment search tool. J Mol Biol. 1990 Oct 5;215(3):403-10. doi: 10.1016/S0022-2836(05)80360-2. PMID: 2231712.

- If using `gget blat`, please also cite:   
  Kent WJ. BLAT--the BLAST-like alignment tool. Genome Res. 2002 Apr;12(4):656-64. doi: 10.1101/gr.229202. PMID: 11932250; PMCID: PMC187518.

- If using `gget muscle`, please also cite:  
  Edgar RC (2021), MUSCLE v5 enables improved estimates of phylogenetic tree confidence by ensemble bootstrapping, bioRxiv 2021.06.20.449169. https://doi.org/10.1101/2021.06.20.449169.
  
- If using `gget alphafold`, please also cite:  
  Jumper, J., Evans, R., Pritzel, A. et al. Highly accurate protein structure prediction with AlphaFold. Nature 596, 583–589 (2021). https://doi.org/10.1038/s41586-021-03819-2
  
  And, if applicable:  
  Evans, R. et al. Protein complex prediction with AlphaFold-Multimer. bioRxiv 2021.10.04.463034; https://doi.org/10.1101/2021.10.04.463034
  
___
# Disclaimer  
`gget` is only as accurate as the databases/servers/APIs it queries from. The accuracy or reliability of the data is not guaranteed or warranted in any way and the providers disclaim liability of any kind whatsoever, including, without limitation, liability for quality, performance, merchantability and fitness for a particular purpose arising out of the use, or inability to use the data.
