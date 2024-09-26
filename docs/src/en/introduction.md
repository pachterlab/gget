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

This grid shows the `gget` core modules, each with a short description. Click on any module to access detailed documentation.

<table style="width:100%; table-layout:fixed;">
  <tr>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="docs/module1.md"><strong>gget alphafold</strong></a><br>Predict 3D protein structure from an amino acid sequence.</td>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="docs/module2.md"><strong>gget archs4</strong></a><br>What is the expression of my gene in tissue X?.</td>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="docs/module3.md"><strong>gget bgee</strong></a><br>Find gene orthologs.</td>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="docs/module4.md"><strong>gget blast</strong></a><br>BLAST a nucleotide or amino acid sequence.</td>
  </tr>
  <tr>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="docs/module5.md"><strong>gget blat</strong></a><br>Find the genome and genomic location of a nucleotide or amino acid sequence.</td>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="docs/module6.md"><strong>gget cbio</strong></a><br>Explore a gene's expression in certain cancers.</td>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="docs/module7.md"><strong>gget cellxgene</strong></a><br>Get ready-to-use single-cell RNA seq count matrices for certain tissues/disease/etc.</td>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="docs/module8.md"><strong>gget cosmic</strong></a><br>Search for genes, mutations, and other factors associated with certain cancers.</td>
  </tr>
  <tr>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="docs/module9.md"><strong>gget diamond</strong></a><br>Align amino acid sequences to a reference.</td>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="docs/module10.md"><strong>gget elm</strong></a><br>Find protein interaction domains and functions in an amino acid sequence.</td>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="docs/module11.md"><strong>gget enrichr</strong></a><br>Check if a list of genes is associated with a specific celltype/pathway/disease/etc</td>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="docs/module12.md"><strong>gget info</strong></a><br>Fetch all of the information associated with an Ensembl ID.</td>
  </tr>
  <tr>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="docs/module13.md"><strong>gget muscle</strong></a><br>Align multiple nucleotide or amino acid sequences to each other.</td>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="docs/module14.md"><strong>gget mutate</strong></a><br>Mutate nucleotide sequences based on specified mutations.</td>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="docs/module15.md"><strong>gget opentargets</strong></a><br>Explore which diseases and drugs a gene is associated with.</td>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="docs/module16.md"><strong>gget pdb</strong></a><br>Fetch data from the Protein Data Bank (PDB) based on a PDB ID.</td>
  </tr>
  <tr>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="docs/module17.md"><strong>gget ref</strong></a><br>Get reference genomes from Ensembl.</td>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="docs/module18.md"><strong>gget search</strong></a><br>Find Ensembl IDs associated with the specified search word.</td>
    <td style="width:25%; padding:20px; text-align:center; vertical-align:top;"><a href="docs/module19.md"><strong>gget seq</strong></a><br>Fetch the nucleotide or amino acid sequence(s) of a gene.</td>
  </tr>
</table>


### [More tutorials](https://github.com/pachterlab/gget_examples)

<br>  

If you use `gget` in a publication, please [cite*](/gget/en/cite.md):    
```
Luebbert, L., & Pachter, L. (2023). Efficient querying of genomic reference databases with gget. Bioinformatics. https://doi.org/10.1093/bioinformatics/btac836
```
Read the article here: [https://doi.org/10.1093/bioinformatics/btac836](https://doi.org/10.1093/bioinformatics/btac836)

<br>  
<iframe width="560" height="315" src="https://www.youtube.com/embed/cVR0k6Mt97o?si=BJwRyaymmxF9w65f" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen></iframe>  



