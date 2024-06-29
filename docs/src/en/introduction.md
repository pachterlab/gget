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
The databases queried by `gget` are continuously being updated which sometimes changes their structure. `gget` modules are tested automatically on a biweekly basis and updated to match new database structures when necessary. If you encounter a problem, please upgrade to the latest `gget` version using `pip install --upgrade gget`. If the problem persists, please [report the issue](https://github.com/pachterlab/gget/issues/new/choose).  
<br>
[<kbd> <br> Request a new feature <br> </kbd>](https://github.com/pachterlab/gget/issues/new/choose)
<br>
<br>
> `gget info` and `gget seq` are currently unable to fetch information for WormBase and FlyBase IDs (all other IDs are functioning normally). This issue arose due to a bug in Ensembl release 112. We appreciate Ensembl's efforts in addressing this issue and expect a fix soon. Thank you for your patience.
<br>

[<img src="https://github.com/pachterlab/gget/blob/main/docs/assets/website_v2_gget_alphafold.png?raw=true" width="32%" height="32%" />](/gget/en/alphafold.md)
[<img src="https://github.com/pachterlab/gget/blob/main/docs/assets/website_v2_gget_archs4.png?raw=true" width="32%" height="32%" />](/gget/en/archs4.md)
[<img src="https://github.com/pachterlab/gget/blob/main/docs/assets/website_v2_gget_blast.png?raw=true" width="32%" height="32%" />](/gget/en/blast.md)  

[<img src="https://github.com/pachterlab/gget/blob/main/docs/assets/website_v2_gget_blat.png?raw=true" width="32%" height="32%" />](/gget/en/blat.md)
[<img src="https://github.com/pachterlab/gget/blob/main/docs/assets/website_v2_gget_cellxgene.png?raw=true" width="32%" height="32%" />](/gget/en/cellxgene.md)
[<img src="https://github.com/pachterlab/gget/blob/main/docs/assets/website_v2_gget_cosmic.png?raw=true" width="32%" height="32%" />](/gget/en/cosmic.md)  

[<img src="https://github.com/pachterlab/gget/blob/main/docs/assets/website_v2_gget_diamond.png?raw=true" width="32%" height="32%" />](/gget/en/diamond.md)
[<img src="https://github.com/pachterlab/gget/blob/main/docs/assets/website_v2_gget_elm.png?raw=true" width="32%" height="32%" />](/gget/en/elm.md)
[<img src="https://github.com/pachterlab/gget/blob/main/docs/assets/website_v2_gget_enrichr.png?raw=true" width="32%" height="32%" />](/gget/en/enrichr.md)  

[<img src="https://github.com/pachterlab/gget/blob/main/docs/assets/website_v2_gget_info.png?raw=true" width="32%" height="32%" />](/gget/en/info.md)
[<img src="https://github.com/pachterlab/gget/blob/main/docs/assets/website_v2_gget_muscle.png?raw=true" width="32%" height="32%" />](/gget/en/muscle.md)
[<img src="https://github.com/pachterlab/gget/blob/main/docs/assets/website_v2_gget_pdb.png?raw=true" width="32%" height="32%" />](/gget/en/pdb.md)  

[<img src="https://github.com/pachterlab/gget/blob/main/docs/assets/website_v2_gget_ref.png?raw=true" width="32%" height="32%" />](/gget/en/ref.md)
[<img src="https://github.com/pachterlab/gget/blob/main/docs/assets/website_v2_gget_search.png?raw=true" width="32%" height="32%" />](/gget/en/search.md)
[<img src="https://github.com/pachterlab/gget/blob/main/docs/assets/website_v2_gget_seq.png?raw=true" width="32%" height="32%" />](/gget/en/seq.md) 

### [More tutorials](https://github.com/pachterlab/gget_examples)

<br>  

If you use `gget` in a publication, please [cite*](/gget/en/cite.md):    
```
Luebbert, L., & Pachter, L. (2023). Efficient querying of genomic reference databases with gget. Bioinformatics. https://doi.org/10.1093/bioinformatics/btac836
```
Read the article here: [https://doi.org/10.1093/bioinformatics/btac836](https://doi.org/10.1093/bioinformatics/btac836)

<br>  
<iframe width="560" height="315" src="https://www.youtube.com/embed/cVR0k6Mt97o?si=BJwRyaymmxF9w65f" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen></iframe>  

<br>  
<iframe src="https://pachterlab.github.io/gget/docs/assets/get_modules_overview.html" width="800" height="600"></iframe>



