> Python arguments are equivalent to long-option arguments (`--arg`), unless otherwise specified. Flags are True/False arguments in Python. The manual for any gget tool can be called from the command-line using the `-h` `--help` flag.  
## gget elm 🎭
Locally predict Eukaryotic Linear Motifs from an amino acid sequence or UniProt ID using data from the [ELM database](http://elm.eu.org/).    
Return format: JSON (command-line) or data frame/CSV (Python). This module returns two data frames (or JSON formatted files) (see examples).     

Before using `gget elm` for the first time, run `gget setup elm` / `gget.setup("elm")` once (also see [`gget setup`](setup.md)).   

**Positional argument**  
`sequence`  
Amino acid sequence or Uniprot ID (str).  
When providing a Uniprot ID, use flag `--uniprot` (Python: `uniprot==True`).  

**Optional arguments**  
`-s` `--sensitivity`  
Sensitivity of DIAMOND alignment (str). Default: "very-sensitive".   
One of the following: fast, mid-sensitive, sensitive, more-sensitive, very-sensitive, or ultra-sensitive.  

`-t` `--threads`  
Number of threads used in DIAMOND alignment (int). Default: 1.  

`-bin` `--diamond_binary`  
Path to DIAMOND binary (str). Default: None -> Uses DIAMOND binary installed with `gget`.  

`-o` `--out`   
Path to the folder to save results in (str), e.g. "path/to/directory". Default: Standard out; temporary files are deleted.   

**Flags**  
`-u` `--uniprot`  
Set to True if `sequence` is a Uniprot ID instead of an amino acid sequence. Default: False.  

`-csv` `--csv`  
Command-line only. Returns results in CSV format.  
Python: Use `json=True` to return output in JSON format.

`-q` `--quiet`   
Command-line only. Prevents progress information from being displayed.  
Python: Use `verbose=False` to prevent progress information from being displayed.   

### Examples
Find ELMs in an amino acid sequence:  
```bash
gget setup elm          # Downloads/updates local ELM database
gget elm -o gget_elm_results LIAQSIGQASFV
```
```python
# Python
gget.setup(“elm”)      # Downloads/updates local ELM database
ortholog_df, regex_df = gget.elm("LIAQSIGQASFV")
```
  
Find ELMs giving a UniProt ID as input:  
```bash
gget setup elm          # Downloads/updates local ELM database
gget elm -o gget_elm_results --uniprot Q02410
```
```python
# Python
gget.setup(“elm”)      # Downloads/updates local ELM database
ortholog_df, regex_df = gget.elm("Q02410", uniprot=True)
```
&rarr; Returns two data frames (or JSON formatted dictionaries for command line) containing extensive information about linear motifs associated with orthologous proteins and motifs found in the input sequence directly based on their regex expressions:  

ortholog_df:  
  
|Ortholog_UniProt_ID|ProteinName|class_accession|ELMIdentifier  |FunctionalSiteName                   |Description                                                                                                                              |Organism    |…  |
|:-----------------:|:---------:|:-------------:|:-------------:|:-----------------------------------:|:---------------------------------------------------------------------------------------------------------------------------------------:|:----------:|:-:|
|Q02410             |APBA1_HUMAN|ELME000357     |LIG_CaMK_CASK_1|CASK CaMK domain binding ligand motif|Motif that mediates binding to the calmodulin-dependent protein kinase (CaMK) domain of the peripheral plasma membrane protein CASK/Lin2.|Homo sapiens|…  |
|Q02410             |APBA1_HUMAN|ELME000091     |LIG_PDZ_Class_2|PDZ domain ligands                   |The C-terminal class 2 PDZ-binding motif is classically represented by a pattern such as                                                 |Homo sapiens|…  |

regex_df:  
  
|Instance_accession|ELMIdentifier     |FunctionalSiteName             |ELMType|Description                                                                                                                                            |Instances (Matched Sequence)|Organism                      |…  |
|:----------------:|:----------------:|:-----------------------------:|:-----:|:-----------------------------------------------------------------------------------------------------------------------------------------------------:|:--------------------------:|:----------------------------:|:-:|
|ELME000321        |CLV_C14_Caspase3-7|Caspase cleavage motif         |CLV    |Caspase-3 and Caspase-7 cleavage site.                                                                                                                 |ERSDG                       |Mus musculus                  |…  |
|ELME000102        |CLV_NRD_NRD_1     |NRD cleavage site              |CLV    |N-Arg dibasic convertase (NRD/Nardilysin) cleavage site.                                                                                               |RRA                         |Rattus norvegicus             |…  |
|ELME000100        |CLV_PCSK_PC1ET2_1 |PCSK cleavage site             |CLV    |NEC1/NEC2 cleavage site.                                                                                                                               |KRD                         |Mus musculus                  |…  |
|ELME000146        |CLV_PCSK_SKI1_1   |PCSK cleavage site             |CLV    |Subtilisin/kexin isozyme-1 (SKI1) cleavage site.                                                                                                       |RLLTA                       |Homo sapiens                  |…  |
|ELME000231        |DEG_APCC_DBOX_1   |APCC-binding Destruction motifs|DEG    |An RxxL-based motif that binds to the Cdh1 and Cdc20 components of APC/C thereby targeting the protein for destruction in a cell cycle dependent manner|SRVKLNIVR                   |Saccharomyces cerevisiae S288c|…  |
|…                 |…                 |…                              |…      |…                                                                                                                                                      |…                           |…                             |…  |

(Motifs that occur in many different species might look repeated, but all rows should be unique.)

#### [More examples](https://github.com/pachterlab/gget_examples)
