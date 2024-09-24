> Python arguments are equivalent to long-option arguments (`--arg`), unless otherwise specified. Flags are True/False arguments in Python. The manual for any gget tool can be called from the command-line using the `-h` `--help` flag.  
## gget bgee 🐝
Fetch orthology and gene expression data from [Bgee](https://www.bgee.org/) using Ensembl IDs.   
Return format: JSON/CSV (command-line) or data frame (Python).

> If you are specifically interested in human gene expression data, consider using [gget opentargets](./opentargets.md) instead.
> **gget bgee** has less data, but supports more species.

**Positional argument**  
`ens_id`  
Ensembl gene ID, e.g. ENSG00000169194 or ENSSSCG00000014725.

**Required arguments**  
`-t` `--type`  
Type of data to fetch. Options: `orthologs`, `expression`.  

**Optional arguments**  
`-o` `--out`    
Path to the JSON file the results will be saved in, e.g. path/to/directory/results.json. Default: Standard out.

**Flags**   
`-csv` `--csv`  
Command-line only. Returns the output in CSV format, instead of JSON format.  
Python: Use `json=True` to return output in JSON format.

`-q` `--quiet`   
Command-line only. Prevents progress information from being displayed.  
Python: Use `verbose=False` to prevent progress information from being displayed.
  
  
### Examples

**Get orthologs for a gene**

```bash
gget bgee ENSSSCG00000014725 -t orthologs
```

```python
import gget
gget.bgee("ENSSSCG00000014725", type="orthologs")
```

&rarr; Returns orthologs for the gene with Ensembl ID ENSSSCG00000014725.

| gene_id            | gene_name    | species_id | genus   | species    |
|--------------------|--------------|------------|---------|------------|
| 734881             | hbb1         | 8355       | Xenopus | laevis     |
| ENSFCAG00000038029 | LOC101098159 | 9685       | Felis   | catus      |
| ENSBTAG00000047356 | LOC107131172 | 9913       | Bos     | taurus     |
| ENSOARG00000019163 | LOC101105437 | 9940       | Ovis    | aries      |
| ENSXETG00000025667 | hbg1         | 8364       | Xenopus | tropicalis |
| ...                | ...          | ...        | ...     | ...        |

<br/><br/>

**Get gene expression data for a gene**

```bash
gget bgee ENSSSCG00000014725 -t expression
```
```python
import gget
gget.bgee("ENSSSCG00000014725", type="expression")
```

&rarr; Returns gene expression data for the gene with Ensembl ID ENSSSCG00000014725.

| anat_entity_id | anat_entity_name            | score | score_confidence | expression_state |
|----------------|-----------------------------|-------|------------------|------------------|
| UBERON:0000178 | blood                       | 99.98 | high             | expressed        |
| UBERON:0002106 | spleen                      | 99.96 | high             | expressed        |
| UBERON:0002190 | subcutaneous adipose tissue | 99.70 | high             | expressed        |
| UBERON:0005316 | endocardial endothelium     | 99.61 | high             | expressed        |
| UBERON:0002107 | liver                       | 99.27 | high             | expressed        |
| ...            | ...                         | ...   | ...              | ...              |

    
#### [More examples](https://github.com/pachterlab/gget_examples)
