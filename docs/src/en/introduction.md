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

<!-- Module Names and Descriptions -->
{% set modules = [
    {"name": "Module 1", "description": "Short description of Module 1", "link": "docs/module1.md"},
    {"name": "Module 2", "description": "Short description of Module 2", "link": "docs/module2.md"},
    {"name": "Module 3", "description": "Short description of Module 3", "link": "docs/module3.md"},
    {"name": "Module 4", "description": "Short description of Module 4", "link": "docs/module4.md"},
    {"name": "Module 5", "description": "Short description of Module 5", "link": "docs/module5.md"},
    {"name": "Module 6", "description": "Short description of Module 6", "link": "docs/module6.md"},
    {"name": "Module 7", "description": "Short description of Module 7", "link": "docs/module7.md"},
    {"name": "Module 8", "description": "Short description of Module 8", "link": "docs/module8.md"},
    {"name": "Module 9", "description": "Short description of Module 9", "link": "docs/module9.md"},
    {"name": "Module 10", "description": "Short description of Module 10", "link": "docs/module10.md"},
    {"name": "Module 11", "description": "Short description of Module 11", "link": "docs/module11.md"},
    {"name": "Module 12", "description": "Short description of Module 12", "link": "docs/module12.md"},
    {"name": "Module 13", "description": "Short description of Module 13", "link": "docs/module13.md"},
    {"name": "Module 14", "description": "Short description of Module 14", "link": "docs/module14.md"},
    {"name": "Module 15", "description": "Short description of Module 15", "link": "docs/module15.md"},
    {"name": "Module 16", "description": "Short description of Module 16", "link": "docs/module16.md"},
    {"name": "Module 17", "description": "Short description of Module 17", "link": "docs/module17.md"},
    {"name": "Module 18", "description": "Short description of Module 18", "link": "docs/module18.md"},
    {"name": "Module 19", "description": "Short description of Module 19", "link": "docs/module19.md"},
    {"name": "Module 20", "description": "Short description of Module 20", "link": "docs/module20.md"}
] %}

<!-- Table Layout for Modules -->
<table style="width:100%; table-layout:fixed;">
  <tr>
    <td style="width:25%; padding:20px; text-align:center;">
      <a href="{{ modules[0].link }}"><strong>{{ modules[0].name }}</strong></a><br>{{ modules[0].description }}
    </td>
    <td style="width:25%; padding:20px; text-align:center;">
      <a href="{{ modules[1].link }}"><strong>{{ modules[1].name }}</strong></a><br>{{ modules[1].description }}
    </td>
    <td style="width:25%; padding:20px; text-align:center;">
      <a href="{{ modules[2].link }}"><strong>{{ modules[2].name }}</strong></a><br>{{ modules[2].description }}
    </td>
    <td style="width:25%; padding:20px; text-align:center;">
      <a href="{{ modules[3].link }}"><strong>{{ modules[3].name }}</strong></a><br>{{ modules[3].description }}
    </td>
  </tr>
  <tr>
    <td style="width:25%; padding:20px; text-align:center;">
      <a href="{{ modules[4].link }}"><strong>{{ modules[4].name }}</strong></a><br>{{ modules[4].description }}
    </td>
    <td style="width:25%; padding:20px; text-align:center;">
      <a href="{{ modules[5].link }}"><strong>{{ modules[5].name }}</strong></a><br>{{ modules[5].description }}
    </td>
    <td style="width:25%; padding:20px; text-align:center;">
      <a href="{{ modules[6].link }}"><strong>{{ modules[6].name }}</strong></a><br>{{ modules[6].description }}
    </td>
    <td style="width:25%; padding:20px; text-align:center;">
      <a href="{{ modules[7].link }}"><strong>{{ modules[7].name }}</strong></a><br>{{ modules[7].description }}
    </td>
  </tr>
  <tr>
    <td style="width:25%; padding:20px; text-align:center;">
      <a href="{{ modules[8].link }}"><strong>{{ modules[8].name }}</strong></a><br>{{ modules[8].description }}
    </td>
    <td style="width:25%; padding:20px; text-align:center;">
      <a href="{{ modules[9].link }}"><strong>{{ modules[9].name }}</strong></a><br>{{ modules[9].description }}
    </td>
    <td style="width:25%; padding:20px; text-align:center;">
      <a href="{{ modules[10].link }}"><strong>{{ modules[10].name }}</strong></a><br>{{ modules[10].description }}
    </td>
    <td style="width:25%; padding:20px; text-align:center;">
      <a href="{{ modules[11].link }}"><strong>{{ modules[11].name }}</strong></a><br>{{ modules[11].description }}
    </td>
  </tr>
  <tr>
    <td style="width:25%; padding:20px; text-align:center;">
      <a href="{{ modules[12].link }}"><strong>{{ modules[12].name }}</strong></a><br>{{ modules[12].description }}
    </td>
    <td style="width:25%; padding:20px; text-align:center;">
      <a href="{{ modules[13].link }}"><strong>{{ modules[13].name }}</strong></a><br>{{ modules[13].description }}
    </td>
    <td style="width:25%; padding:20px; text-align:center;">
      <a href="{{ modules[14].link }}"><strong>{{ modules[14].name }}</strong></a><br>{{ modules[14].description }}
    </td>
    <td style="width:25%; padding:20px; text-align:center;">
      <a href="{{ modules[15].link }}"><strong>{{ modules[15].name }}</strong></a><br>{{ modules[15].description }}
    </td>
  </tr>
  <tr>
    <td style="width:25%; padding:20px; text-align:center;">
      <a href="{{ modules[16].link }}"><strong>{{ modules[16].name }}</strong></a><br>{{ modules[16].description }}
    </td>
    <td style="width:25%; padding:20px; text-align:center;">
      <a href="{{ modules[17].link }}"><strong>{{ modules[17].name }}</strong></a><br>{{ modules[17].description }}
    </td>
    <td style="width:25%; padding:20px; text-align:center;">
      <a href="{{ modules[18].link }}"><strong>{{ modules[18].name }}</strong></a><br>{{ modules[18].description }}
    </td>
    <td style="width:25%; padding:20px; text-align:center;">
      <a href="{{ modules[19].link }}"><strong>{{ modules[19].name }}</strong></a><br>{{ modules[19].description }}
    </td>
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



