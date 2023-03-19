`gget` currently consists of the following modules:  
- [`gget alphafold`](alphafold.md)  
Predict the 3D structure of a protein from its amino acid sequence using a simplified version of [DeepMind](https://www.deepmind.com/)’s [AlphaFold2](https://github.com/deepmind/alphafold).  
  
- [`gget archs4`](archs4.md)  
Find the most correlated genes to a gene of interest or find the gene's tissue expression atlas using [ARCHS4](https://maayanlab.cloud/archs4/).  
  
- [`gget blast`](blast.md)  
BLAST a nucleotide or amino acid sequence to any [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) database.  
  
- [`gget blat`](blat.md)  
Find the genomic location of a nucleotide or amino acid sequence using [BLAT](https://genome.ucsc.edu/cgi-bin/hgBlat).  
  
- [`gget enrichr`](enrichr.md)  
Perform an enrichment analysis on a list of genes using [Enrichr](https://maayanlab.cloud/Enrichr/).  
  
- [`gget gpt`](gpt.md)  
Generates text based on a given prompt using the [OpenAI API](https://openai.com/)'s 'openai.ChatCompletion.create' endpoint.  
  
- [`gget info`](info.md)  
Fetch extensive gene and transcript metadata from [Ensembl](https://www.ensembl.org/), [UniProt](https://www.uniprot.org/), and [NCBI](https://www.ncbi.nlm.nih.gov/) using Ensembl IDs.  
  
- [`gget muscle`](muscle.md)  
Align multiple nucleotide or amino acid sequences to each other using [Muscle5](https://www.drive5.com/muscle/).  
  
- [`gget pdb`](pdb.md)  
Get the structure and metadata of a protein from the [RCSB Protein Data Bank](https://www.rcsb.org/).  
  
- [`gget ref`](ref.md)  
Fetch File Transfer Protocols (FTPs) and metadata for reference genomes and annotations from [Ensembl](https://www.ensembl.org/) by species.  
  
- [`gget search`](search.md)   
Fetch genes and transcripts from [Ensembl](https://www.ensembl.org/) using free-form search terms.  
  
- [`gget setup`](setup.md)  
Helper module to install/download third-party dependencies for a specified gget module.  
  
- [`gget seq`](seq.md)  
Fetch nucleotide or amino acid sequences of genes or transcripts from [Ensembl](https://www.ensembl.org/) or [UniProt](https://www.uniprot.org/), respectively.  