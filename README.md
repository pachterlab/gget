# gget
[![pypi version](https://img.shields.io/pypi/v/gget)](https://pypi.org/project/gget)
![PyPI - Downloads](https://img.shields.io/pypi/dm/gget)
[![license](https://img.shields.io/pypi/l/gget)](LICENSE)
[![DOI](https://zenodo.org/badge/458943224.svg)](https://zenodo.org/badge/latestdoi/458943224)
![Code Coverage](https://img.shields.io/badge/Coverage-87%25-green.svg)

**gget** consists of nine tools:
- [**gget ref**](#gget-ref)  
Fetch FTPs and metadata for reference genomes and annotations from [Ensembl](https://www.ensembl.org/) by species.
- [**gget search**](#gget-search)   
Fetch genes and transcripts from [Ensembl](https://www.ensembl.org/) using free-form search terms.
- [**gget info**](#gget-info)  
Fetch extensive gene and transcript metadata from [Ensembl](https://www.ensembl.org/), [UniProt](https://www.uniprot.org/), and [NCBI](https://www.ncbi.nlm.nih.gov/) using Ensembl IDs.
- [**gget seq**](#gget-seq)  
Fetch nucleotide or amino acid sequences of genes or transcripts from [Ensembl](https://www.ensembl.org/) or [UniProt](https://www.uniprot.org/), respectively.
- [**gget blast**](#gget-blast)  
BLAST a nucleotide or amino acid sequence against any [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) database.
- [**gget blat**](#gget-blat)  
Find the genomic location of a nucleotide or amino acid sequence using [BLAT](https://genome.ucsc.edu/cgi-bin/hgBlat).
- [**gget muscle**](#gget-muscle)  
Align multiple nucleotide or amino acid sequences against each other using [Muscle5](https://www.drive5.com/muscle/).
- [**gget enrichr**](#gget-enrichr)  
Perform an enrichment analysis on a list of genes using [Enrichr](https://maayanlab.cloud/Enrichr/).
- [**gget archs4**](#gget-archs4)  
Find the most correlated genes or the tissue expression atlas of a gene of interest using [ARCHS4](https://maayanlab.cloud/archs4/).

## Installation
```bash
pip install gget
```

For use in Jupyter Lab / Google Colab:
```python
import gget
```

## Quick start guide
```bash
# Fetch all Homo sapiens reference and annotation FTPs from the latest Ensembl release
$ gget ref -s homo_sapiens

# Search human genes with "ace2" AND "angiotensin" in their name/description
$ gget search -sw ace2,angiotensin -s homo_sapiens -ao and 

# Look up gene ENSG00000130234 (ACE2) with expanded info (returns all transcript isoforms for genes)
$ gget info -id ENSG00000130234 -e

# Fetch the amino acid sequence of the canonical transcript of gene ENSG00000130234
$ gget seq -id ENSG00000130234 --seqtype transcript

# Quickly find the genomic location of (the start of) that amino acid sequence
$ gget blat -seq MSSSSWLLLSLVAVTAAQSTIEEQAKTFLDKFNHEAEDLFYQSSLAS

# Blast (the start of) that amino acid sequence
$ gget blast -seq MSSSSWLLLSLVAVTAAQSTIEEQAKTFLDKFNHEAEDLFYQSSLAS

# Align nucleotide or amino acid sequences stored in a FASTA file
$ gget muscle -fa path/to/file.fa

# Use Enrichr to find the ontology of a list of genes
$ gget enrichr -g ACE2 AGT AGTR1 ACE AGTRAP AGTR2 ACE3P -db ontology

# Get the human tissue expression atlas of gene ACE2
$ gget archs4 -g ACE2 -w tissue
```
Jupyter Lab / Google Colab:
```python  
gget.ref("homo_sapiens")
gget.search(["ace2", "angiotensin"], "homo_sapiens", andor="and")
gget.info("ENSG00000130234", expand=True)
gget.seq("ENSG00000130234", seqtype="transcript")
gget.blat("MSSSSWLLLSLVAVTAAQSTIEEQAKTFLDKFNHEAEDLFYQSSLAS")
gget.blast("MSSSSWLLLSLVAVTAAQSTIEEQAKTFLDKFNHEAEDLFYQSSLAS")
gget.muscle("path/to/file.fa")
gget.enrichr(["ACE2", "AGT", "AGTR1", "ACE", "AGTRAP", "AGTR2", "ACE3P"], database="ontology", plot=True)
gget.archs4("ACE2", which="tissue")
```
___

# Instruction Manual
Jupyter Lab / Google Colab arguments are equivalent to long-option arguments (`--arg`).
The manual to any gget tool can be called from terminal using a `-h` `--help` flag.

## gget ref
Fetch FTPs and their respective metadata (or use flag `ftp` to only return the links) for reference genomes and annotations from [Ensembl](https://www.ensembl.org/) by species.
Return format: dictionary/json.

#### Required arguments
`-s` `--species`  
Species for which the FTPs will be fetched in the format genus_species, e.g. homo_sapiens.
Note: Not required when calling flag `--list_species`.  

#### Optional arguments
`-w` `--which`
Defines which results to return. Default: 'all' -> Returns all available results.
Possible entries are one or a combination of the following:
'gtf'   - Returns the annotation (GTF).
'cdna'  - Returns the trancriptome (cDNA).
'dna'   - Returns the genome (DNA).
'cds    - Returns the coding sequences corresponding to Ensembl genes. (Does not contain UTR or intronic sequence.)
'cdrna' - Returns transcript sequences corresponding to non-coding RNA genes (ncRNA).
'pep'   - Returns the protein translations of Ensembl genes.

`-r` `--release`
Defines the Ensembl release number from which the files are fetched, e.g. 104. Default: latest Ensembl release.

`-o` `--out`  
Path to the json file the results will be saved in, e.g. path/to/directory/results.json. Default: Standard out.
Jupyter Lab / Google Colab: `save=True` will save the output in the current working directory.

#### Flags
`-l` `--list_species`  
Lists all available species. (Jupyter Lab / Google Colab: combine with 'species=None'.)

`-ftp` `--ftp`  
Returns only the requested FTP links.

`-d` `--download`  
Downloads the requested FTPs to the current directory (requires [curl](https://curl.se/docs/) to be installed).

#### [Examples](https://github.com/pachterlab/gget_examples)
___

## gget search 
Fetch genes and transcripts from [Ensembl](https://www.ensembl.org/) using free-form search terms. 
Return format: data frame.

#### Required arguments
`-sw` `--searchwords`  
One or more free form search terms, e.g. gaba, nmda. (Note: Search is not case-sensitive.)

`-s` `--species`  
Species or database to be searched. 
A species can be passed in the format 'genus_species', e.g. 'homo_sapiens'. 
To pass a specific database, pass the name of the CORE database, e.g. 'mus_musculus_dba2j_core_105_1'. 
All availabale databases can be found [here](http://ftp.ensembl.org/pub/release-106/mysql/).

#### Optional arguments
`-st` `--seqtype`  
'gene' (default) or 'transcript'
Returns genes or transcripts, respectively.

`-ao` `--andor`  
'or' (default) or 'and'
'or': Returns all genes that INCLUDE AT LEAST ONE of the searchwords in their name/description.
'and': Returns only genes that INCLUDE ALL of the searchwords in their name/description.

`-l` `--limit`  
Limits the number of search results, e.g. 10. Default: None.

`-o` `--out`  
Path to the csv the results will be saved in, e.g. path/to/directory/results.csv. Default: Standard out.  
Jupyter Lab / Google Colab: `save=True` will save the output in the current working directory.
    
#### [Examples](https://github.com/pachterlab/gget_examples)
___

## gget info  
Fetch extensive gene and transcript metadata from [Ensembl](https://www.ensembl.org/), [UniProt](https://www.uniprot.org/), and [NCBI](https://www.ncbi.nlm.nih.gov/) using Ensembl IDs.
Return format: data frame.

#### Required arguments
`-id` `--ens_ids`  
One or more Ensembl IDs.

#### Optional arguments
`-o` `--out`  
Path to the csv the results will be saved in, e.g. path/to/directory/results.csv. Default: Standard out.  
Jupyter Lab / Google Colab: `save=True` will save the output in the current working directory.

#### Flags
`-e` `--expand`  
Expands returned information (only for gene and transcript IDs). 
For genes, adds information on all known transcripts.
For transcripts, adds information on all known translations and exons.

`wrap_text`
Jupyter Lab / Google Colab only. `wrap_text=True` displays data frame with wrapped text for easy reading (default: False). 

#### [Examples](https://github.com/pachterlab/gget_examples)
___

## gget seq  
Fetch nucleotide or amino acid sequence of a gene (and all its isoforms) or a transcript by Ensembl ID. 
Return format: FASTA.

#### Required arguments
`-id` `--ens_ids`  
One or more Ensembl IDs.

#### Optional arguments
`-st` `--seqtype`
'gene' (default) or 'transcript'.
Defines whether nucleotide or amino acid sequences are returned.
Nucleotide sequences are fetched from [Ensembl](https://www.ensembl.org/).
Amino acid sequences are fetched from [UniProt](https://www.uniprot.org/).

`-o` `--out`  
Path to the file the results will be saved in, e.g. path/to/directory/results.fa. Default: Standard out. 
Jupyter Lab / Google Colab: `save=True` will save the output in the current working directory.

#### Flags
`-i` `--isoforms`  
Returns the sequences of all known transcripts (for `seqtype=gene` only).

#### [Examples](https://github.com/pachterlab/gget_examples)
___

## gget blast
BLAST a nucleotide or amino acid sequence against any [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) database.
Return format: data frame.

#### Required arguments
`-seq` `--sequence`  
Nucleotise or amino acid sequence, or path to FASTA file.

#### Optional arguments
`-p` `--program`
'blastn', 'blastp', 'blastx', 'tblastn', or 'tblastx'.
Default: 'blastn' for nucleotide sequences; 'blastp' for amino acid sequences.

`-db` `--database`
'nt', 'nr', 'refseq_rna', 'refseq_protein', 'swissprot', 'pdbaa', or 'pdbnt'.
Default: 'nt' for nucleotide sequences; 'nr' for amino acid sequences.
[More info on BLAST databases](https://ncbi.github.io/blast-cloud/blastdb/available-blastdbs.html)

`-d` `--descriptions`
Limits number of descriptions returned. Default: 500.

`-a` `--alignments`
Limits number of alignments to return. Default: 500.

`-hs` `--hitlist_size`
Limits number of hits to return. Default: 50.

`-e` `--expect`
Defines the [expect value](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=FAQ#expect) cutoff. Default: 10.0.

`-o` `--out`  
Path to the csv the results will be saved in, e.g. path/to/directory/results.csv. Default: Standard out.  
Jupyter Lab / Google Colab: `save=True` will save the output in the current working directory.

#### Flags
`-ng` `--ncbi_gi`
Also returns NCBI GI identifiers.

`-lcf` `--low_comp_filt`
Turns on [low complexity filter](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=FAQ#LCR).

`-mbo` `--megablast_off`
Turns off MegaBLAST algorithm. Default: MegaBLAST on (blastn only).

`-q` `--quiet` 
Prevents progress information from being displayed.

`wrap_text`
Jupyter Lab / Google Colab only. `wrap_text=True` displays data frame with wrapped text for easy reading (default: False). 
