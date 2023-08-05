[![pypi version](https://img.shields.io/pypi/v/gget)](https://pypi.org/project/gget)
[![image](https://anaconda.org/bioconda/gget/badges/version.svg)](https://anaconda.org/bioconda/gget)
[![Downloads](https://static.pepy.tech/personalized-badge/gget?period=total&units=international_system&left_color=grey&right_color=brightgreen&left_text=Downloads)](https://pepy.tech/project/gget)
[![Conda](https://img.shields.io/conda/dn/bioconda/gget?logo=Anaconda)](https://anaconda.org/bioconda/gget)
[![license](https://img.shields.io/pypi/l/gget)](LICENSE)
![status](https://github.com/pachterlab/gget/workflows/CI/badge.svg)
![status](https://github.com/lauraluebbert/test_gget_alphafold/workflows/CI_alphafold/badge.svg)
[![Star on GitHub](https://img.shields.io/github/stars/pachterlab/gget.svg?style=social)](https://github.com/pachterlab/gget/)  

<img align="left" width="50%" height="50%" src="https://github.com/pachterlab/gget/blob/main/figures/gget_overview.png?raw=true" />
  
# Welcome!
  
`gget` is a free, open-source command-line tool and Python package that enables efficient querying of genomic databases.  
<br>
`gget` consists of a collection of separate but interoperable modules, each designed to facilitate one type of database querying in a single line of code.   
<br>  
<br>  
<br>  
<br>  

[<img src="https://github.com/pachterlab/gget/assets/56094636/fbeb2662-7e97-47eb-a0eb-e27b68a6f256" width="30%" height="30%" />](alphafold.md)
[<img src="https://github.com/pachterlab/gget/assets/56094636/feda5cec-89b4-4ebc-83c3-40bab363ca72" width="30%" height="30%" />](archs4.md)
[<img src="https://github.com/pachterlab/gget/assets/56094636/02336057-bda2-4555-9696-d65789533331" width="30%" height="30%" />](blast.md)  





| Module |     Description  |
|:--------------|:------------------|

| [`gget blat`](blat.md) | *Find the location of a nucleotide or amino acid sequence in a genome* using [BLAT](https://genome.ucsc.edu/cgi-bin/hgBlat). BLAT is designed <br> to quickly find sequences of >95% similarity.  |
| [`gget cellxgene`](cellxgene.md)   |  *Query single-cell transcriptomics data by tissue, gene, and more* from [CZ CELLxGENE Discover](https://cellxgene.cziscience.com/). <br> [CZ CELLxGENE Discover](https://cellxgene.cziscience.com/) comprises hundreds of standardized data collections and millions of cells <br>  characterizing the functionality of mouse and human tissues. `gget cellxgene` returns an AnnData <br>  object containing the requested count matrix and metadata.|
| [`gget enrichr`](enrichr.md)  | *Perform an enrichment analysis on a list of target and background genes* using [Enrichr](https://maayanlab.cloud/Enrichr/). Gene set <br> enrichment analysis (GSEA) identifies classes of genes that are over-represented in a large set of<br>  genes or proteins associated with different phenotypes (e.g. cell types, diseases, pathways, etc.). <br> [Enrichr](https://maayanlab.cloud/Enrichr/) facilitates enrichment analysis against a diverse collection of >200 GSEA databases. |
| [`gget gpt`](gpt.md) | *Generate text based on a text prompt* using [OpenAI](https://openai.com/). |
|[`gget info`](info.md)  | *Fetch extensive gene and transcript metadata associated with an Ensembl ID* including <br> transcript and isoform information. `gget info` combines data from [Ensembl](https://www.ensembl.org/), [UniProt](https://www.uniprot.org/), <br> and [NCBI](https://www.ncbi.nlm.nih.gov/).  |
|[`gget muscle`](muscle.md)  | *Align multiple nucleotide or amino acid sequences to each other* using the [Muscle5](https://www.drive5.com/muscle/) algorithm.  <br> [Muscle5](https://www.drive5.com/muscle/)  constructs an ensemble of high-accuracy alignment with diverse biases by perturbing <br> a hidden Markov model and permuting its guide tree.|
 |[`gget pdb`](pdb.md)  | *Get the structure and metadata of a protein* from the [RCSB Protein Data Bank](https://www.rcsb.org/). `gget pdb` returns <br> 3D structure data deposited in the [RCSB PDB](https://www.rcsb.org/) database in Protein Data Bank (PDB) format.
 |[`gget ref`](ref.md)  | *Fetch File Transfer Protocol (FTP) links and metadata for a species' reference genome*  from <br> [Ensembl](https://www.ensembl.org/). `gget ref` will fetch data on the latest genome by default, but also supports the <br> specification of an Ensembl release. It can also download the genome directly.
 |[`gget search`](search.md)   | *Search for genes and transcripts from [Ensembl](https://www.ensembl.org/) using free-form search terms.*  <br>  `gget search` supports fetching data from a specific Ensembl release or database <br> (e.g. different mouse strains).
  |[`gget seq`](seq.md)  | *Fetch the nucleotide or amino acid sequence(s) of a gene or transcript* from its Ensembl ID.<br>  The nucleotide or amino acid sequences are fetched from [Ensembl](https://www.ensembl.org/) or [UniProt](https://www.uniprot.org/), respectively. <br> `gget seq` can fetch the canonical transcript or all isoforms.

<br>


If you use `gget` in a publication, please [cite*](cite.md):    
```
Luebbert, L., & Pachter, L. (2023). Efficient querying of genomic reference databases with gget. Bioinformatics. https://doi.org/10.1093/bioinformatics/btac836
```
Read the article here: [https://doi.org/10.1093/bioinformatics/btac836](https://doi.org/10.1093/bioinformatics/btac836)

<br>
<br>

<img src="https://user-images.githubusercontent.com/56094636/222949999-0b89cba2-134f-4cbe-acbb-8f20b3f52684.jpg" alt="" width="250" height="160" />
