# Examples

## gget FetchTP usage examples

### Fetch *Homo sapiens* GTF, DNA, and cDNA FTPs from latest Ensemble release
Jupyter Lab:
```python
fetchtp("homo_sapiens", save=True)
```
Ouput in file: [JL_fetchtp_human.json](https://github.com/lauraluebbert/gget/blob/main/examples/JL_fetchtp_human.json)

Terminal:
```
gget fetchtp -sp homo_sapiens
```
Ouput in file: [term_fetchtp_human.json](https://github.com/lauraluebbert/gget/blob/main/examples/term_fetchtp_human.json)

### Fetch from a previous Ensembl release, e.g. release 103
Jupyter Lab:
```python
fetchtp("homo_sapiens", release=103, save=True)
```
Ouput in file: [JL_fetchtp_human_r103.json](https://github.com/lauraluebbert/gget/blob/main/examples/JL_fetchtp_human_r103.json)

Terminal:
```
gget fetchtp -sp homo_sapiens -r 103
```

Ouput in file: [term_fetchtp_human_r103.json](https://github.com/lauraluebbert/gget/blob/main/examples/term_fetchtp_human_r103.json)


## gget search usage examples

### Query Ensembl for *Homo sapiens* genes whose description contains the searchwords "gaba", "nmda", or "ampa"
Jupyter Lab:
```python
gget(["gaba", "nmda", "ampa"], "human", save=True)
```
Output in file: [JL_gget_human_gaba_nmda_ampa.csv](https://github.com/lauraluebbert/gget/blob/main/examples/JL_gget_human_gaba_nmda_ampa.csv)

Terminal:
```
gget search -sw gaba nmda ampa -sp human
```
Output in file: [term_gget_human_gaba_nmda_ampa.csv](https://github.com/lauraluebbert/gget/blob/main/examples/term_gget_human_gaba_nmda_ampa.csv)

Note: gget search allows the following species abbreviations:
"homo_sapiens" - "human"
"mus_musculus" - "mouse"
"taeniopygia_guttata" - "zebra finch"
"caenorhabditis_elegans" - "roundworm"
All other species have to be called using their specific database, as shown in the example below. All 236 databases can be found [here](http://ftp.ensembl.org/pub/release-105/mysql/).

### Query Ensembl for *Mus musculus* genes whose description contains the searchword "mitochondrial", but limit the results to 10 genes
Jupyter Lab:
 ```python
gget("mitochondrial", "mus_musculus", limit=10, save=True)
```
Output in file: [JL_gget_mouse_mito.csv](https://github.com/lauraluebbert/gget/blob/main/examples/JL_gget_mouse_mito.csv)

Terminal:
```
gget search -sw mitochondrial -sp mus_musculus -l 10
```
Output in file: [term_gget_mouse_mito.csv](https://github.com/lauraluebbert/gget/blob/main/examples/term_gget_mouse_mito.csv)

### Query Ensembl for genes from the killifish (*Nothobranchius furzeri*) database whose description contains the searchword "brain"
Jupyter Lab:
```python
gget("brain", "nothobranchius_furzeri_core_105_2", save=True)
```
Output in file: [JL_gget_killifish_brain.csv](https://github.com/lauraluebbert/gget/blob/main/examples/JL_gget_killifish_brain.csv)

Terminal:
```
gget search -sw brain -sp nothobranchius_furzeri_core_105_2
```
Output in file: [term_gget_killifish_brain.csv](https://github.com/lauraluebbert/gget/blob/main/examples/term_gget_killifish_brain.csv)

