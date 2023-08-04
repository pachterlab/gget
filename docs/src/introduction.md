[![pypi version](https://img.shields.io/pypi/v/gget)](https://pypi.org/project/gget)
[![image](https://anaconda.org/bioconda/gget/badges/version.svg)](https://anaconda.org/bioconda/gget)
[![Downloads](https://static.pepy.tech/personalized-badge/gget?period=total&units=international_system&left_color=grey&right_color=brightgreen&left_text=Downloads)](https://pepy.tech/project/gget)
[![Conda](https://img.shields.io/conda/dn/bioconda/gget?logo=Anaconda)](https://anaconda.org/bioconda/gget)
[![license](https://img.shields.io/pypi/l/gget)](LICENSE)
![status](https://github.com/pachterlab/gget/workflows/CI/badge.svg)
![status](https://github.com/lauraluebbert/test_gget_alphafold/workflows/CI_alphafold/badge.svg)
[![Star on GitHub](https://img.shields.io/github/stars/pachterlab/gget.svg?style=social)](https://github.com/pachterlab/gget/)  

<p text-align=justify>`gget` is a free, open-source command-line tool and Python package that enables efficient querying of genomic databases. `gget`  consists of a collection of separate but interoperable modules, each designed to facilitate one type of database querying in a single line of code.  </p> 
 
<img src="https://github.com/pachterlab/gget/blob/main/figures/gget_overview.png?raw=true)"  display="block"
        margin-left="auto" margin-right= "auto" width= 50% height= 50%  flex-direction=column> 

    
| Module        |     Description  |
|---------------|------------------|
| [`gget alphafold`](alphafold.md) |    Predict the 3D structure of a protein from its amino acid sequence using a simplified version of [DeepMind](https://www.deepmind.com/)’s [AlphaFold2](https://github.com/deepmind/alphafold). |
|  [`gget archs4`](archs4.md)   |  Find the most correlated genes to a gene of interest or find the gene's tissue expression atlas using [ARCHS4](https://maayanlab.cloud/archs4/). Archs4 is a resource that provides access to gene and transcript counts uniformly processed from all human and mouse RNA-seq experiments from the Gene Expression Omnibus (GEO) and the Sequence Read Archive (SRA). | 
|  [`gget blast`](blast.md)   | BLAST a nucleotide or amino acid sequence to any [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) database.  BLAST finds regions of similarity between biological sequences. The program compares nucleotide or protein sequences to sequence databases and calculates the statistical significance.|
| [`gget blat`](blat.md) | Find the genomic location of a nucleotide or amino acid sequence using [BLAT](https://genome.ucsc.edu/cgi-bin/hgBlat). BLAT on DNA is designed to quickly find sequences of 95% and greater similarity of length 25 bases or more.  |
| [`gget cellxgene`](cellxgene.md)   |  Query data from [CZ CELLxGENE Discover](https://cellxgene.cziscience.com/) using the [CZ CELLxGENE Discover Census](https://github.com/chanzuckerberg/cellxgene-census). Cellxgene is an interactive, performant explorer for single cell transcriptomics data.|
| [`gget enrichr`](enrichr.md)  | Perform an enrichment analysis on a list of genes using [Enrichr](https://maayanlab.cloud/Enrichr/). Enrichr is an interactive and collaborative HTML5 gene list enrichment analysis tool.|
| [`gget gpt`](gpt.md) |Generates text based on a given prompt using the [OpenAI API](https://openai.com/)'s 'openai.ChatCompletion.create' endpoint. |
|[`gget info`](info.md)  |Fetch extensive gene and transcript metadata from [Ensembl](https://www.ensembl.org/), [UniProt](https://www.uniprot.org/), and [NCBI](https://www.ncbi.nlm.nih.gov/) using Ensembl IDs.  |
|[`gget muscle`](muscle.md)  | Align multiple nucleotide or amino acid sequences to each other using [Muscle5](https://www.drive5.com/muscle/). Muscle5 is a novel algorithm which constructs an ensemble of high-accuracy alignment with diverse biases by perturbing a hidden Markov model and permuting its guide tree.|
 |[`gget pdb`](pdb.md)  | Get the structure and metadata of a protein from the [RCSB Protein Data Bank](https://www.rcsb.org/).  
 |[`gget ref`](ref.md)  | Fetch File Transfer Protocols (FTPs) and metadata for reference genomes and annotations from [Ensembl](https://www.ensembl.org/) by species.  
 |[`gget search`](search.md)   | Fetch genes and transcripts from [Ensembl](https://www.ensembl.org/) using free-form search terms.  
  [`gget setup`](setup.md)  |Helper module to install/download third-party dependencies for a specified gget module.  
  |[`gget seq`](seq.md)  |Fetch nucleotide or amino acid sequences of genes or transcripts from [Ensembl](https://www.ensembl.org/) or [UniProt](https://www.uniprot.org/), respectively.  

If you use `gget` in a publication, please [cite*](cite.md):    
```
Luebbert, L., & Pachter, L. (2023). Efficient querying of genomic reference databases with gget. Bioinformatics. https://doi.org/10.1093/bioinformatics/btac836
```
Read the article here: [https://doi.org/10.1093/bioinformatics/btac836](https://doi.org/10.1093/bioinformatics/btac836)

<br>
<br>
<br>
<br>
<br>

<img src="https://user-images.githubusercontent.com/56094636/222949999-0b89cba2-134f-4cbe-acbb-8f20b3f52684.jpg" alt="" width="250" height="160" />

