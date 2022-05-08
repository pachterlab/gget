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

# Instruction Manual
## gget ref
Fetch FTPs and metadata for reference genomes and annotations from [Ensembl](https://www.ensembl.org/) by species. Returns a dictionary/json containing the requested URLs with their respective metadata (or use flag `ftp` to only return the links).

### Options
`-l` `--list`  
List all available species.

`-s` `--species`  
Species for which the FTPs will be fetched in the format genus_species, e.g. homo_sapiens.

`-w` `--which`  
Defines which results to return.
Default: 'all' -> Returns all available results.
Possible entries are one or a combination (as a list of strings) of the following:
'gtf' - Returns the annotation (GTF).
'cdna' - Returns the trancriptome (cDNA).
'dna' - Returns the genome (DNA).
'cds - Returns the coding sequences corresponding to Ensembl genes. (Does not contain UTR or intronic sequence.)
'cdrna' - Returns transcript sequences corresponding to non-coding RNA genes (ncRNA).
'pep' - Returns the protein translations of Ensembl genes.

`-r` `--release`  
Defines the Ensembl release number from which the files are fetched, e.g. 104.

`-ftp` `--ftp`  
Return only the requested FTP links in a list (default: False).

`-d` `--download`  
Download the requested FTPs to the current directory (requires [curl](https://curl.se/docs/)).

`-o` `--out`  
Path to the json file the results will be saved in, e.g. path/to/directory/results.json. Default: Standard out.
For Jupyter Lab / Google Colab: `save=True` will save the output to the current working directory.

### Examples
#### Use in combination with [kallisto | bustools](https://www.kallistobus.tools/kb_usage/kb_ref/) to build a reference

```bash
kb ref -i INDEX -g T2G -f1 FASTA $(gget ref --ftp -w dna,gtf -s homo_sapiens)
```
&rarr; kb ref builds a reference using the latest DNA and GTF files of species **Homo sapiens** passed to it by gget ref

#### Show all available species

```python
# Jupyter Lab / Google Colab:
!gget ref --list

# Terminal:
$ gget ref --list
```
&rarr; Returns a list with all available species from the latest Ensembl release.

#### Fetch GTF, DNA, and cDNA FTP links for a specific species

```python
# Jupyter Lab / Google Colab:
ref("homo_sapiens")

# Terminal:
$ gget ref -s homo_sapiens
```
&rarr; Returns a json with the latest links to human GTF and FASTA files, their respective release dates and time, and the Ensembl release from which the links were fetched, in the format:
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
  
Note: To download the files linked to by the FTPs into the current directory, add flag `-d`.

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

___
## gget search 
Query [Ensembl](https://www.ensembl.org/) for genes or transcripts from a defined species using free form search words.  
> :warning: **gget search currently only supports genes listed in the Ensembl core API, which includes limited external references.** Manually searching the [Ensembl website](https://uswest.ensembl.org/) might yield more results.

### Options
`-sw` `--searchwords`  
One or more free form searchwords for the query, e.g. gaba, nmda. Searchwords are not case-sensitive.

`-s` `--species`  
Species or database to be searched.  
Species can be passed in the format 'genus_species', e.g. 'homo_sapiens'. To pass a specific CORE database (e.g. a specific mouse strain), enter the name of the CORE database, e.g. 'mus_musculus_dba2j_core_105_1'. All availabale species databases can be found here: http://ftp.ensembl.org/pub/release-105/mysql/

`-t` `--d_type`  
Possible entries: 'gene' (default), 'transcript'
Returns either genes or transcripts, respectively, which match the searchwords.

`-ao` `--andor`  
Possible entries: 'or', 'and'
'or': ID descriptions must include at least one of the searchwords (default). 
'and': Only return IDs whose descriptions include all searchwords.

`-l` `--limit`  
Limits the number of search results to the top `[limit]` genes found (default: None).

`-o` `--out`  
Path to the file the results will be saved in, e.g. path/to/directory/results.csv (default: None &rarr; just prints results).  
For Jupyter Lab / Google Colab: `save=True` will save the output to the current working directory.

### Examples
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


#### Query Ensembl for transcripts from a specific species which include ALL searchwords
```python
# Jupyter Lab / Google Colab:
search(["gaba", "gamma-aminobutyric"], "nothobranchius_furzeri", d_type="transcript", andor="and")

# Terminal 
$ gget search -sw gaba,gamma-aminobutyric -s nothobranchius_furzeri -t transcript -ao and
```
&rarr; Returns all killifish transcripts that contain all of the searchwords in their Ensembl or external reference description.

#### Query Ensembl for genes from a specific species using a single searchword and while limiting the number of returned search results 

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

___
## gget info  
Look up gene or transcript Ensembl IDs. Returns their common name, description, homologs, synonyms, corresponding transcript/gene, transcript isoforms and more from the Ensembl database as well as external references.

### Options
`-id` `--ens_ids`  
One or more Ensembl IDs.

`-e` `--expand`  
Expand returned information (default: False). For genes: add isoform information. For transcripts: add translation and exon information.

`-H` `--homology`  
Returns homology information of ID (default: False).

`-x` `--xref`  
Returns information from external references (default: False).

`-o` `--out`  
Path to the file the results will be saved in, e.g. path/to/directory/results.json (default: None &rarr; just prints results).  
For Jupyter Lab / Google Colab: `save=True` will save the output to the current working directory.

### Examples
#### Look up a list of gene Ensembl IDs including information on all isoforms

```python
# Jupyter Lab / Google Colab:
info(["ENSG00000034713", "ENSG00000104853", "ENSG00000170296"], expand=True)

# Terminal 
$ gget info -id ENSG00000034713,ENSG00000104853,ENSG00000170296 -e
```
&rarr; Returns a json containing information about each ID, amongst others the common name, description, and corresponding transcript/gene, in the format:
```
{
            "Ensembl ID": {
                        "species": genus_species,
                        "object_type": e.g. Gene,
                        "biotype": Gene biotype, e.g. protein_coding,
                        "display_name": Common gene name,
                        "description": Ensemble description,
                        "assembly_name": Name of species assmebly,
                        "seq_region_name": Sequence region,
                        "start": Sequence start position,
                        "end": Sequence end position,
                        "strand": Strand
                        "canonical_transcript": Transcript ID,
                        # All transcript isoforms:
                        "Transcript": [{'display_name': Transcript name,
					'biotype': Transcript biotype,
					'id': Transcript ID}, ...]
                        },
}
```
Note: When looking up Ensembl IDs of transcripts instead of genes, the "Transcript" entry above will be replaced by "Translation" and "Exon" information. 

#### Look up a transcript Ensembl ID and include external reference descriptions

```python
# Jupyter Lab / Google Colab:
info("ENSDART00000135343", xref=True)

# Terminal 
$ gget info -id ENSDART00000135343 -x
```
&rarr; Returns a json containing the homology information, and external reference description of each ID in addition to the standard information mentioned above.

___
## gget seq  
Fetch DNA sequences from gene or transcript Ensembl IDs. 

### Options
`-id` `--ens_ids`  
One or more Ensembl IDs.

`-i` `--isoforms`  
If a gene Ensembl ID is passed, this returns sequences of all known transcript isoforms.

`-o` `--out`  
Path to the file the results will be saved in, e.g. path/to/directory/results.fa (default: None &rarr; just prints results).  
For Jupyter Lab / Google Colab: `save=True` will save the output FASTA to the current working directory.

### Examples
#### Fetch the sequences of several transcript Ensembl IDs
```python
# Jupyter Lab / Google Colab:
seq(["ENST00000441207","ENST00000587537"])

# Terminal 
$ gget seq -id ENST00000441207,ENST00000587537
```
&rarr; Returns a FASTA containing the sequence of each ID, in the format:
```
>Ensembl_ID chromosome:assembly:seq_region_name:seq_region_start:seq_region_end:strand
GGGAATGGAAATCTGTCCCTCGTGCTGGAAGCCAACCAGTGGTGATGACTCTGTGTGCCACTCCGCCTCCTACAGCGCGGATCCTCTG  
CGTGTGTCCTCGCAAGACAAGCTCGATGAAATGGCCGAGTCCAGTCAAGCAAACTTTGAGGGAA...
```

#### Fetch the sequences of a gene Ensembl ID and all its transcript isoforms
```python
# Jupyter Lab / Google Colab:
seq("ENSMUSG00000025040", isoforms=True)

# Terminal 
$ gget seq -id ENSMUSG00000025040 -i
```
&rarr; Returns a FASTA containing the sequence of the gene ID and the sequences of all of each transcripts.
