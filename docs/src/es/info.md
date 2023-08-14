> Parámetros de Python són iguales a los parámetros largos (`--parámetro`) de Terminal, si no especificado de otra manera. Banderas son parámetros de verdadero o falso (True/False) en Python. El manuál para cualquier modulo de gget se puede llamar desde la Terminal con la bandera `-h` `--help`.  
## gget info 💡
Obtenga información detallada sobre genes y transcripciones de [Ensembl](https://www.ensembl.org/), [UniProt](https://www.uniprot.org/) y [NCBI](https://www. ncbi.nlm.nih.gov/) utilizando sus IDs del tipo Ensembl.  
Regresa: Resultados en formato JSON (Terminal) o Dataframe/CSV (Python).  

**Parámetro posicional**  
`ens_ids`   
Uno o más ID del tipo Ensembl.  

**Parámetros optionales**  
`-o` `--out`   
Ruta al archivo en el que se guardarán los resultados, p. ruta/al/directorio/resultados.csv (o .json). Por defecto: salida estándar (STDOUT).  
Para Python, usa `save=True` para guardar los resultados en el directorio de trabajo actual.  

**Banderas**  
`-n` `--ncbi`  
DESACTIVA los resultados de [NCBI](https://www.ncbi.nlm.nih.gov/).  
Para Python: `ncbi=False` evita la incluida de datos de NCBI (por defecto: True).    

`-u` `--uniprot`  
DESACTIVA los resultados de [UniProt](https://www.uniprot.org/).  
Para Python: `uniprot=False` evita la incluida de datos de UniProt (por defecto: True).   

`-pdb` `--pdb`  
INCLUYE [PDB](https://www.ebi.ac.uk/pdbe/) IDs en los resultados (podría aumentar el tiempo de ejecución).  
Para Python: `pdb=True` incluye IDs de PDB en los resultados (por defecto: False). 

`-csv` `--csv`  
Solo para la Terminal. Regresa los resultados en formato CSV.    
Para Python, usa `json=True` para regresar los resultados en formato JSON.  

`-q` `--quiet`   
Solo para la Terminal. Impide la informacion de progreso de ser exhibida durante la corrida.  
Para Python, usa `verbose=False` para imipidir la informacion de progreso de ser exhibida durante la corrida.  

`wrap_text`  
Solo para Python. `wrap_text=True` muestra los resultados con texto envuelto para facilitar la lectura (por defecto: False).  


### Por ejemplo
```bash
gget info ENSG00000034713 ENSG00000104853 ENSG00000170296
```
```python
# Python
gget.info(["ENSG00000034713", "ENSG00000104853", "ENSG00000170296"])
```
&rarr; Regresa información detallada sobre cada ID de Ensembl ingresada:

|      | uniprot_id     | ncbi_gene_id     | primary_gene_name | synonyms | protein_names | ensembl_description | uniprot_description | ncbi_description | biotype | canonical_transcript | ... |
| -------------- |-------------------------| ------------------------| -------------- | ----------|-----|----|----|----|----|----|----|
| ENSG00000034713| P60520 | 11345 | GABARAPL2 | [ATG8, ATG8C, FLC3A, GABARAPL2, GATE-16, GATE16, GEF-2, GEF2] | Gamma-aminobutyric acid receptor-associated protein like 2 (GABA(A) receptor-associated protein-like 2)... | GABA type A receptor associated protein like 2 [Source:HGNC Symbol;Acc:HGNC:13291] | FUNCTION: Ubiquitin-like modifier involved in intra- Golgi traffic (By similarity). Modulates intra-Golgi transport through coupling between NSF activity and ... | Enables ubiquitin protein ligase binding activity. Involved in negative regulation of proteasomal protein catabolic process and protein... | protein_coding | ENST00000037243.7 |... |
| . . .            | . . .                     | . . .                     | . . .            | . . .       | . . . | . . . | . . . | . . . | . . . | . . . | ... |
  
#### [More examples](https://github.com/pachterlab/gget_examples)
