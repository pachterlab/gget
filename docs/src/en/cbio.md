> Python arguments are equivalent to long-option arguments (`--arg`), unless otherwise specified. Flags are True/False arguments in Python. The manual for any gget tool can be called from the command-line using the `-h` `--help` flag.  
## gget cbio 📖
Plot cancer genomics heatmaps using data from [cBioPortal](https://www.cbioportal.org/) using Ensembl IDs or gene names.

**Positional argument**  
`subcommand`  
Either `search` or `plot`

### `search` subcommand (Python: `gget.cbio_search`)
Find cBioPortal study IDs by keyword.   
Return format: JSON (command-line) or string list (Python).     
**Note: This does not return studies with mixed cancer types.**

**Positional argument**     
`keywords`  
Space-separated list of keywords to search for, e.g. <code>breast&nbsp;lung</code>.

### `plot` subcommand (Python: `gget.cbio`)
Plot cancer genomics heatmaps using data from cBioPortal.
Return format: PNG (command-line and Python)

**Required arguments**  
`-s` `--study_ids`    
Space-separated list of cBioPortal study IDs, e.g. <code>msk_impact_2017&nbsp;egc_msk_2023</code>.

`-g` `--genes`  
Space-separated list of gene names or Ensembl IDs, e.g. <code>NOTCH3&nbsp;ENSG00000108375</code>.

**Optional arguments**  
`-st` `--stratification`    
Column to stratify the data by. Default: `tissue`.  
Options:
- `tissue`
- `cancer_type`
- `cancer_type_detailed`
- `study_id`
- `sample`

`-vt` `--variation_type`    
Type of variation to plot. Default: `mutation_occurrences`.     
Options:
- `mutation_occurrences`
- `cna_nonbinary` Note: `stratification` must be `sample` for this option.
- `sv_occurrences`
- `cna_occurrences`
- `Consequence` Note: `stratification` must be `sample` for this option.

`-f` `--filter`     
Filter the data by a specific value in a specific column, e.g. `study_id:msk_impact_2017`   
Python: `filter_=(column, value)`

`-dd` `--data_dir`  
Directory to store data files. Default: `./gget_cbio_cache`.

`-fd` `--figure_dir`    
Directory to output figures. Default: `./gget_cbio_figures`.

`-dpi` `--dpi`  
DPI of the output figure. Default: 100.

**Flags**   

`-q` `--quiet`   
Command-line only. Prevents progress information from being displayed.  
Python: Use `verbose=False` to prevent progress information from being displayed. 

`-y` `--no_confirm`     
Command-line only. Skip download confirmation prompts.  
Python: Use `confirm_download=True` to enable download confirmation prompts.

`-sh` `--show`  
Show the plot in a window/in a Jupyter notebook.
  
  
### Examples

**Find all cBioPortal studies with cancer types matching specific keywords:**   
```bash
gget cbio search esophag ovary ovarian
```
```python
# Python
import gget
gget.cbio_search('esophag', 'ovary', 'ovarian')
```
&rarr; Returns a list of studies with cancer types matching the keywords `esophag`, `ovary`, or `ovarian`.

```
['egc_tmucih_2015', 'egc_msk_2017', ..., 'msk_spectrum_tme_2022']
```

<br/><br/>

**Plot a heatmap of mutation occurrences for specific genes in a specific study:**   
```bash
gget cbio plot -s msk_impact_2017 -g AKT1 ALK FLT4 MAP3K1 MLL2 MLL3 NOTCH3 NOTCH4 PDCD1 RNF43 -st tissue -vt mutation_occurrences -dpi 200 -y
```
```python
import gget
gget.cbio(
    ['msk_impact_2017'],
    ['AKT1', 'ALK', 'FLT4', 'MAP3K1', 'MLL2', 'MLL3', 'NOTCH3', 'NOTCH4', 'PDCD1', 'RNF43'],
    stratification='tissue',
    variation_type='mutation_occurrences',
    dpi=200
)
```

&rarr; Saves a heatmap of mutation occurrences for the specified genes in the specified study to `./gget_cbio_figures/Heatmap_tissue.png`.

![Heatmap](https://raw.githubusercontent.com/techno-sam/gget/f6b4465eecae0f07c71558f8f6e7b60d36c8d41a/docs/assets/gget_cbio_figure_1.png)

<br/><br/>

**Plot a heatmap of mutation types for specific genes in a specific study:**   
```bash
gget cbio plot -s msk_impact_2017 -g AKT1 ALK FLT4 MAP3K1 MLL2 MLL3 NOTCH3 NOTCH4 PDCD1 RNF43 -st sample -vt Consequence -dpi 200 -y
```
```python
import gget
gget.cbio(
    ['msk_impact_2017'],
    ['AKT1', 'ALK', 'FLT4', 'MAP3K1', 'MLL2', 'MLL3', 'NOTCH3', 'NOTCH4', 'PDCD1', 'RNF43'],
    stratification='sample',
    variation_type='Consequence',
    dpi=200,
)
```

&rarr; Saves a heatmap of mutation types for the specified genes in the specified study to `./gget_cbio_figures/Heatmap_sample.png`.

![Heatmap](https://raw.githubusercontent.com/techno-sam/gget/f6b4465eecae0f07c71558f8f6e7b60d36c8d41a/docs/assets/gget_cbio_figure_2.png)
    
#### [More examples](https://github.com/pachterlab/gget_examples)
