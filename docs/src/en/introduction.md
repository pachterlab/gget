[![pypi version](https://img.shields.io/pypi/v/gget)](https://pypi.org/project/gget)
[![image](https://anaconda.org/bioconda/gget/badges/version.svg)](https://anaconda.org/bioconda/gget)
[![Downloads](https://static.pepy.tech/personalized-badge/gget?period=total&units=international_system&left_color=grey&right_color=brightgreen&left_text=Downloads)](https://pepy.tech/project/gget)
[![Conda](https://img.shields.io/conda/dn/bioconda/gget?logo=Anaconda)](https://anaconda.org/bioconda/gget)
[![license](https://img.shields.io/pypi/l/gget)](LICENSE)
![status](https://github.com/pachterlab/gget/actions/workflows/ci.yml/badge.svg)
![status](https://github.com/lauraluebbert/test_gget_alphafold/actions/workflows/CI_alphafold.yml/badge.svg)
[![Star on GitHub](https://img.shields.io/github/stars/pachterlab/gget.svg?style=social)](https://github.com/pachterlab/gget/)  

[<img align="right" width="50%" height="50%" src="https://github.com/pachterlab/gget/blob/main/docs/assets/website_v2_gget_overview.png?raw=true" />](https://raw.githubusercontent.com/pachterlab/gget/main/figures/gget_overview.png)

# Welcome!
  
`gget` is a free, open-source command-line tool and Python package that enables efficient querying of genomic databases.  
<br>
`gget` consists of a collection of separate but interoperable modules, each designed to facilitate one type of database querying in a single line of code.   
<br>
NOTE: The databases queried by `gget` are continuously being updated which sometimes changes their structure. `gget` modules are tested automatically on a biweekly basis and updated to match new database structures when necessary. If you encounter a problem, please upgrade to the latest `gget` version using `pip install --upgrade gget`. If the problem persists, please [report the issue](https://github.com/pachterlab/gget/issues/new/choose).  
<br>
[<kbd> <br> Request a new feature <br> </kbd>](https://github.com/pachterlab/gget/issues/new/choose)
<br>
<br>
> `gget info` and `gget seq` are currently unable to fetch information for WormBase and FlyBase IDs (all other IDs are functioning normally). This issue arose due to a bug in Ensembl release 112. We appreciate Ensembl's efforts in addressing this issue and expect a fix soon. Thank you for your patience.

# gget modules

These are the `gget` core modules. Click on any module to access detailed documentation.

<table style="width:100%; table-layout:fixed;">
  <tr>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/en/alphafold.md"><strong>gget alphafold</strong></a><br>Predict 3D protein structure from an amino acid sequence.</td>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/en/archs4.md"><strong>gget archs4</strong></a><br>What is the expression of my gene in tissue X?</td>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/en/bgee.md"><strong>gget bgee</strong></a><br>Find all orthologs of a gene.</td>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/en/blast.md"><strong>gget blast</strong></a><br>BLAST a nucleotide or amino acid sequence.</td>
  </tr>
  <tr>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/en/blat.md"><strong>gget blat</strong></a><br>Find the genomic location of a nucleotide or amino acid sequence.</td>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/en/cbio.md"><strong>gget cbio</strong></a><br>Explore a gene's expression in the specified cancers.</td>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/en/cellxgene.md"><strong>gget cellxgene</strong></a><br>Get ready-to-use single-cell RNA seq count matrices from certain tissues/ diseases/ etc.</td>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/en/cosmic.md"><strong>gget cosmic</strong></a><br>Search for genes, mutations, and other factors associated with certain cancers.</td>
  </tr>
  <tr>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/en/diamond.md"><strong>gget diamond</strong></a><br>Align amino acid sequences to a reference.</td>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/en/elm.md"><strong>gget elm</strong></a><br>Find protein interaction domains and functions in an amino acid sequence.</td>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/en/enrichr.md"><strong>gget enrichr</strong></a><br>Check if a list of genes is associated with a specific celltype/ pathway/ disease/ etc.</td>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/en/info.md"><strong>gget info</strong></a><br>Fetch all of the information associated with an Ensembl ID.</td>
  </tr>
  <tr>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/en/muscle.md"><strong>gget muscle</strong></a><br>Align multiple nucleotide or amino acid sequences to each other.</td>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/en/mutate.md"><strong>gget mutate</strong></a><br>Mutate nucleotide sequences based on specified mutations.</td>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/en/opentargets.md"><strong>gget opentargets</strong></a><br>Explore which diseases and drugs a gene is associated with.</td>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/en/pdb.md"><strong>gget pdb</strong></a><br>Fetch data from the Protein Data Bank (PDB) based on a PDB ID.</td>
  </tr>
  <tr>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/en/ref.md"><strong>gget ref</strong></a><br>Get reference genomes from Ensembl.</td>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/en/search.md"><strong>gget search</strong></a><br>Find Ensembl IDs associated with the specified search word.</td>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/en/seq.md"><strong>gget seq</strong></a><br>Fetch the nucleotide or amino acid sequence of a gene.</td>
  </tr>
</table>

<br>  

If you use `gget` in a publication, please [cite*](/gget/en/cite.md):    
```
Luebbert, L., & Pachter, L. (2023). Efficient querying of genomic reference databases with gget. Bioinformatics. https://doi.org/10.1093/bioinformatics/btac836
```
Read the article here: [https://doi.org/10.1093/bioinformatics/btac836](https://doi.org/10.1093/bioinformatics/btac836)

<br>  
<iframe width="560" height="315" src="https://www.youtube.com/embed/cVR0k6Mt97o?si=BJwRyaymmxF9w65f" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen></iframe>  

![logo-bmbf](https://github.com/user-attachments/assets/68449ba0-09d1-4f1e-b747-5ae2fec08a21) 
![logo-okfn](https://github.com/user-attachments/assets/452ae8d8-69f0-4d0d-848c-ddfb40357eb2)
