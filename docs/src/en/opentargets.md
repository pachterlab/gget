> Python arguments are equivalent to long-option arguments (`--arg`), unless otherwise specified. Flags are True/False arguments in Python. The manual for any gget tool can be called from the command-line using the `-h` `--help` flag.  
## gget opentargets 📖
Fetch associated diseases or drugs from [OpenTargets](https://platform.opentargets.org/) using Ensembl IDs.   
Return format: JSON/CSV (command-line) or data frame (Python).

**Positional argument**  
`ens_id`  
Ensembl gene ID, e.g ENSG00000169194.

**Optional arguments**  
`-r` `--resource`   
Defines the type of information to return in the output. Default: 'diseases'.   
Possible resources are:     
'diseases' - Returns associated diseases.
'drugs' - Returns associated drugs.
'tractability' - Returns tractability data.
'pharmacogenetics' - Returns pharmacogenetic responses.
'expression' - Returns gene expression data (by tissues, organs, and anatomical systems).
'depmap' - Returns DepMap gene-disease effect data.

`-l` `--limit`  
Limit the number of results, e.g 10. Default: No limit.     
Note: Not compatible with the 'tractability' and 'depmap' resources.

`-o` `--out`    
Path to the JSON file the results will be saved in, e.g. path/to/directory/results.json. Default: Standard out.  
Python: `save=True` will save the output in the current working directory.

**Flags**   
`-csv` `--csv`  
Command-line only. Returns the output in CSV format, instead of JSON format.

`-q` `--quiet`   
Command-line only. Prevents progress information from being displayed.  
Python: Use `verbose=False` to prevent progress information from being displayed. 

`wrap_text`  
Python only. `wrap_text=True` displays data frame with wrapped text for easy reading (default: False).  
  
  
### Examples

**Get associated diseases for a specific gene:**   
```bash
gget opentargets ENSG00000169194 -l 1
```
```python
# Python
import gget
gget.opentargets('ENSG00000169194', limit=1)
```
&rarr; Returns the top disease associated with the gene ENSG00000169194.

| id            | name          | description                                                           | score            |
|---------------|---------------|-----------------------------------------------------------------------|------------------|
| EFO_0000274   | atopic eczema | A chronic inflammatory genetically determined disease of the skin ... | 0.66364347241831 |

<br/><br/>

**Get associated drugs for a specific gene:**   
```bash
gget opentargets ENSG00000169194 -r drugs -l 2
```
```python
import gget
gget.opentargets('ENSG00000169194', resource='drugs', limit=2)
```

&rarr; Returns the top 2 drugs associated with the gene ENSG00000169194.

| id            | name         | type     | action_mechanism          | description                                                  | synonyms                                           | trade_names           | disease_id  | disease_name             | trial_phase | trial_status | trial_ids     | approved |
|---------------|--------------|----------|---------------------------|--------------------------------------------------------------|----------------------------------------------------|-----------------------|-------------|--------------------------|-------------|--------------|---------------|----------|
| CHEMBL1743081 | TRALOKINUMAB | Antibody | Interleukin-13 inhibitor  | Antibody drug with a maximum clinical trial phase of IV ...  | ['CAT-354', 'Tralokinumab']                        | ['Adbry', 'Adtralza'] | EFO_0000274 | atopic eczema            | 4           |              | []            | True     |
| CHEMBL4297864 | CENDAKIMAB   | Antibody | Interleukin-13 inhibitor  | Antibody drug with a maximum clinical trial phase of III ... | [ABT-308, Abt-308, CC-93538, Cendakimab, RPC-4046] | []                    | EFO_0004232 | eosinophilic esophagitis | 3           | Recruiting   | [NCT04991935] | False    |

*Note: Returned `trial_ids` are [ClinicalTrials.gov](https://clinicaltrials.gov) identifiers*

<br/><br/>

**Get tractability data for a specific gene:**   
```bash
gget opentargets ENSG00000169194 -r tractability
```
```python
import gget
gget.opentargets('ENSG00000169194', resource='tractability')
```

&rarr; Returns tractability data for the gene ENSG00000169194.

| label                 | modality       |
|-----------------------|----------------|
| High-Quality Pocket   | Small molecule |
| Approved Drug         | Antibody       |
| GO CC high conf       | Antibody       |
| UniProt loc med conf  | Antibody       |
| UniProt SigP or TMHMM | Antibody       |

<br/><br/>

**Get pharmacogenetic responses for a specific gene:**
```bash
gget opentargets ENSG00000169194 -r pharmacogenetics -l 1
```
```python
import gget
gget.opentargets('ENSG00000169194', resource='pharmacogenetics', limit=1)
```

&rarr; Returns pharmacogenetic responses for the gene ENSG00000169194.

| rs_id     | genotype_id       | genotype | variant_consequence_id | variant_consequence_label | drugs                                                                                                                                   | phenotype                                                               | genotype_annotation                                                                                          | response_category | direct_target | evidence_level | source   | literature |
|-----------|-------------------|----------|------------------------|---------------------------|-----------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------|-------------------|---------------|----------------|----------|------------|
| rs1295686 | 5_132660151_T_T,T | TT       | SO:0002073             | no_sequence_alteration    | &nbsp;&nbsp;&nbsp;&nbsp;id&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;name<br/>0&nbsp;&nbsp;None&nbsp;&nbsp;hepatitis&nbsp;vaccines | increased risk for non&#8209;immune response to the hepatitis B vaccine | Patients with the TT genotype may be at increased risk for non-immune response to the hepatitis B vaccine... | efficacy          | False         | 3              | pharmgkb | [21111021] |

*Note: Returned `literature` ids are [Europe PMC](https://europepmc.org/article/med/) identifiers*

<br/><br/>

**Get tissues where a gene is most expressed:**
```bash
gget opentargets ENSG00000169194 -r expression -l 2
```
```python
import gget
gget.opentargets('ENSG00000169194', resource='expression', limit=2)
```

&rarr; Returns the top 2 tissues where the gene ENSG00000169194 is most expressed.

| tissue_id      | tissue_name                           | rna_zscore | rna_value | rna_unit | rna_level | anatomical_systems                                                   | organs                                                 |
|----------------|---------------------------------------|------------|-----------|----------|-----------|----------------------------------------------------------------------|--------------------------------------------------------|
| UBERON_0000473 | testis                                | 5          | 1026      |          | 3         | [reproductive&nbsp;system]                                           | [reproductive&nbsp;organ, reproductive&nbsp;structure] |
| CL_0000542     | EBV&#8209;transformed&nbsp;lymphocyte | 1          | 54        |          | 2         | [hemolymphoid&nbsp;system, immune&nbsp;system, lymphoid&nbsp;system] | [immune organ]                                         |

<br/><br/>

**Get DepMap gene-disease effect data for a specific gene:**
```bash
gget opentargets ENSG00000169194 -r depmap
```
```python
import gget
gget.opentargets('ENSG00000169194', resource='depmap')
```

&rarr; Returns DepMap gene-disease effect data for the gene ENSG00000169194.

| depmap_id        | expression | effect   | tissue_id      | tissue_name | cell_line_name | disease_cell_line_id | disease_name    | mutation |
|------------------|------------|----------|----------------|-------------|----------------|----------------------|-----------------|----------|
| ACH&#8209;001532 | 0.176323   | 0.054950 | UBERON_0002113 | kidney      | JMU-RTK-2      | None                 | Rhabdoid Cancer | None     |

    
#### [More examples](https://github.com/pachterlab/gget_examples)
