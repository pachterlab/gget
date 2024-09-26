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
    <td style="width:25%; padding:20px; text-align:center;"><a href="docs/module1.md"><strong>Module 1</strong></a><br>Short description of Module 1.</td>
    <td style="width:25%; padding:20px; text-align:center;"><a href="docs/module2.md"><strong>Module 2</strong></a><br>Short description of Module 2.</td>
    <td style="width:25%; padding:20px; text-align:center;"><a href="docs/module3.md"><strong>Module 3</strong></a><br>Short description of Module 3.</td>
    <td style="width:25%; padding:20px; text-align:center;"><a href="docs/module4.md"><strong>Module 4</strong></a><br>Short description of Module 4.</td>
  </tr>
  <tr>
    <td style="width:25%; padding:20px; text-align:center;"><a href="docs/module5.md"><strong>Module 5</strong></a><br>Short description of Module 5.</td>
    <td style="width:25%; padding:20px; text-align:center;"><a href="docs/module6.md"><strong>Module 6</strong></a><br>Short description of Module 6.</td>
    <td style="width:25%; padding:20px; text-align:center;"><a href="docs/module7.md"><strong>Module 7</strong></a><br>Short description of Module 7.</td>
    <td style="width:25%; padding:20px; text-align:center;"><a href="docs/module8.md"><strong>Module 8</strong></a><br>Short description of Module 8.</td>
  </tr>
  <tr>
    <td style="width:25%; padding:20px; text-align:center;"><a href="docs/module9.md"><strong>Module 9</strong></a><br>Short description of Module 9.</td>
    <td style="width:25%; padding:20px; text-align:center;"><a href="docs/module10.md"><strong>Module 10</strong></a><br>Short description of Module 10.</td>
    <td style="width:25%; padding:20px; text-align:center;"><a href="docs/module11.md"><strong>Module 11</strong></a><br>Short description of Module 11.</td>
    <td style="width:25%; padding:20px; text-align:center;"><a href="docs/module12.md"><strong>Module 12</strong></a><br>Short description of Module 12.</td>
  </tr>
  <tr>
    <td style="width:25%; padding:20px; text-align:center;"><a href="docs/module13.md"><strong>Module 13</strong></a><br>Short description of Module 13.</td>
    <td style="width:25%; padding:20px; text-align:center;"><a href="docs/module14.md"><strong>Module 14</strong></a><br>Short description of Module 14.</td>
    <td style="width:25%; padding:20px; text-align:center;"><a href="docs/module15.md"><strong>Module 15</strong></a><br>Short description of Module 15.</td>
    <td style="width:25%; padding:20px; text-align:center;"><a href="docs/module16.md"><strong>Module 16</strong></a><br>Short description of Module 16.</td>
  </tr>
  <tr>
    <td style="width:25%; padding:20px; text-align:center;"><a href="docs/module17.md"><strong>Module 17</strong></a><br>Short description of Module 17.</td>
    <td style="width:25%; padding:20px; text-align:center;"><a href="docs/module18.md"><strong>Module 18</strong></a><br>Short description of Module 18.</td>
    <td style="width:25%; padding:20px; text-align:center;"><a href="docs/module19.md"><strong>Module 19</strong></a><br>Short description of Module 19.</td>
    <td style="width:25%; padding:20px; text-align:center;"><a href="docs/module20.md"><strong>Module 20</strong></a><br>Short description of Module 20.</td>
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



