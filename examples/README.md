# Examples

## gget ref examples
### Fetch *Homo sapiens* GTF, DNA, and cDNA FTPs from latest Ensemble release
```python
# Jupyter Lab / Google Colab:
ref("human", save=True)

# Terminal
$ gget ref -s human -o ref_human.json
```
Ouput in file: [ref_human.json](https://github.com/lauraluebbert/gget/blob/main/examples/ref_human.json)

Note:  
gget ref and search allow the following species abbreviations:  
"homo_sapiens" &rarr;  "human"  
"mus_musculus" &rarr;  "mouse"  

### Fetch from a previous Ensembl release, e.g. release 103
```python
# Jupyter Lab / Google Colab:
ref("homo_sapiens", release=103, save=True)

# Terminal
$ gget ref -s homo_sapiens -r 103 -o ref_human_r103.json
```
Ouput in file: [ref_human_r103.json](https://github.com/lauraluebbert/gget/blob/main/examples/ref_human_r103.json)


## gget search examples
### Query Ensembl for *Homo sapiens* genes whose description contains one of the searchwords "gaba", "nmda", or "ampa"
```python
# Jupyter Lab / Google Colab:
search(["gaba", "nmda", "ampa"], "human", save=True)

# Terminal
$ gget search -sw gaba nmda ampa -s human -o search_human_gaba_nmda_ampa.csv
```
Output in file: [search_human_gaba_nmda_ampa.csv](https://github.com/lauraluebbert/gget/blob/main/examples/search_human_gaba_nmda_ampa.csv)

### Query Ensembl for killifish transcripts whose description contains the searchword "brain"
 ```python
 # Jupyter Lab / Google Colab:
search("brain", "nothobranchius_furzeri", d_type="transcript", save=True)

# Terminal
$ gget search -sw brain -s nothobranchius_furzeri -t transcript -o search_killifish_brain.csv
```
Output in file: [search_killifish_brain.csv](https://github.com/lauraluebbert/gget/blob/main/examples/search_killifish_brain.csv)

### Query Ensembl for goldfish genes whose description contains the searchwords "brain" AND "glycophorin"
 ```python
 # Jupyter Lab / Google Colab:
search(["blood","glycophorin"], "carassius_auratus", andor="and", save=True)

# Terminal
$ gget search -sw blood,glycophorin -s carassius_auratus -ao and -o search_goldfish_blood.csv
```
Output in file: [search_goldfish_blood.csv](https://github.com/lauraluebbert/gget/blob/main/examples/search_goldfish_blood.csv)

## gget info examples
### Look up a list of gene Ensembl IDs including information on all isoforms
 ```python
# Jupyter Lab / Google Colab:
info(["ENSG00000034713", "ENSG00000104853", "ENSG00000170296"], expand=True, save=True)

# Terminal 
$ gget info -id ENSG00000034713,ENSG00000104853,ENSG00000170296 -e -o info_results.json
```
Output in file: [info_results.json](https://github.com/lauraluebbert/gget/blob/main/examples/info_results.json)


## gget seq examples
### Fetch the sequences of several transcript Ensembl IDs
 ```python
# Jupyter Lab / Google Colab:
seq(["ENST00000441207","ENST00000587537"], save=True)

# Terminal 
$ gget seq -id ENST00000441207,ENST00000587537 -o seq_results.fa
```
Output in file: [seq_results.fa](https://github.com/lauraluebbert/gget/blob/main/examples/seq_results.fa)

### Fetch the sequences of a gene Ensembl ID and all its transcript isoforms
 ```python
# Jupyter Lab / Google Colab:
seq("ENSMUSG00000025040", isoforms=True, save=True)

# Terminal 
$ gget seq -id ENSMUSG00000025040 -i -o seq_iso_results.fa
```
Output in file: [seq_iso_results.fa](https://github.com/lauraluebbert/gget/blob/main/examples/seq_iso_results.fa)



