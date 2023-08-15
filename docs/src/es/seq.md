> Parámetros de Python són iguales a los parámetros largos (`--parámetro`) de Terminal, si no especificado de otra manera. Banderas son parámetros de verdadero o falso (True/False) en Python. El manuál para cualquier modulo de gget se puede llamar desde la Terminal con la bandera `-h` `--help`.  
## gget seq 🧬
Obtenga la(s) secuencia(s) nucleótidos o aminoácidos de un gen (y todas sus isoformas) con su ID de Ensembl.  
Regresa: Archivo de tipo FASTA.  

**Parámetro posicional**  
`ens_ids`   
One or more Ensembl IDs.

**Parámetros optionales**  
`-o` `--out`   
Ruta al archivo en el que se guardarán los resultados, p. ruta/al/directorio/resultados.fa. Por defecto: salida estándar (STDOUT).  
Para Python, usa `save=True` para guardar los resultados en el directorio de trabajo actual.  

**Banderas**  
`-t` `--translate`  
Regresa secuencias de aminoácidos (en lugar de nucleótidos).  
Las secuencias de nucleótidos se obtienen de [Ensembl](https://www.ensembl.org/).  
Las secuencias de aminoácidos se obtienen de [UniProt](https://www.uniprot.org/).  

`-iso` `--isoforms`   
Regresa las secuencias de todas las transcripciones conocidas.  
(Solo para IDs de genes).  

`-q` `--quiet`   
Solo para la Terminal. Impide la informacion de progreso de ser exhibida durante la corrida.  
Para Python, usa `verbose=False` para imipidir la informacion de progreso de ser exhibida durante la corrida.  


### Por ejemplo  
```bash
gget seq ENSG00000034713 ENSG00000104853 ENSG00000170296
```
```python
# Python
gget.seq(["ENSG00000034713", "ENSG00000104853", "ENSG00000170296"])
```
&rarr; Regresa las secuencias de nucleótidos de ENSG00000034713, ENSG00000104853, y ENSG00000170296 en formato FASTA.  

<br/><br/>

```bash
gget seq -t -iso ENSG00000034713
```
```python
# Python
gget.seq("ENSG00000034713", translate=True, isoforms=True)
```
&rarr; Regresa las secuencias de aminoácidos de todas las transcripciones conocidas de ENSG00000034713 en formato FASTA.  

#### [Más ejemplos](https://github.com/pachterlab/gget_examples)
