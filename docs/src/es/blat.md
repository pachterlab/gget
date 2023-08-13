> Parámetros de Python són iguales a los parámetros largos (`--parámetro`) de Terminal, si no especificado de otra manera. Banderas son parámetros de verdadero o falso (True/False) en Python. El manuál para cualquier modulo de gget se puede llamar desde el Terminal con la bandera `-h` `--help`.  
## gget blat 🎯
Encuentra la ubicación genómica de una secuencia de nucleótidos o aminoácidos usando [BLAT](https://genome.ucsc.edu/cgi-bin/hgBlat).   
Regresa: Resultados en formato JSON (Terminal) o Dataframe/CSV (Python).  

**Parámetro posicional**  
`sequence`   
Secuencia de nucleótidos o aminoácidos, o una ruta a un archivo tipo FASTA o .txt.  

**Parámetros optionales**  
`-st` `--seqtype`    
'DNA', 'protein', 'translated%20RNA', o 'translated%20DNA'.   
Por defecto: 'DNA' para secuencias de nucleótidos; 'protein' para secuencias de aminoácidos.  

`-a` `--assembly`    
Ensamblaje del genoma. 'human' (hg38) (esto se usa por defecto), 'mouse' (mm39) (ratón), 'zebrafinch' (taeGut2) (
pinzón cebra),   
o cualquiera de los ensamblajes de especies disponibles [aquí](https://genome.ucsc.edu/cgi-bin/hgBlat) (use el nombre corto del ensamblado, p. 'hg38').  

`-o` `--out`   
Ruta al archivo en el que se guardarán los resultados, p. ruta/al/directorio/resultados.csv (o .json). Por defecto: salida estándar (STDOUT).  
Para Python, usa `save=True` para guardar los resultados en el directorio de trabajo actual.  
  
**Banderas**  
`-csv` `--csv`  
Solo para la Terminal. Regresa los resultados en formato CSV.    
Para Python, usa `json=True` para regresar los resultados en formato JSON.  

`-q` `--quiet`   
Solo para la Terminal. Impide la informacion de progreso de ser exhibida durante la corrida.  
Para Python, usa `verbose=False` para imipidir la informacion de progreso de ser exhibida durante la corrida.  


### Por ejemplo
```bash
gget blat -a taeGut2 MKWMFKEDHSLEHRCVESAKIRAKYPDRVPVIVEKVSGSQIVDIDKRKYLVPSDITVAQFMWIIRKRIQLPSEKAIFLFVDKTVPQSR
```
```python
# Python
gget.blat("MKWMFKEDHSLEHRCVESAKIRAKYPDRVPVIVEKVSGSQIVDIDKRKYLVPSDITVAQFMWIIRKRIQLPSEKAIFLFVDKTVPQSR", assembly="taeGut2")
```
&rarr; Regresa los resultados de BLAT para el ensamblaje taeGut2 (pinzón cebra). En este ejemplo, `gget blat` automáticamente detecta esta secuencia como una secuencia de aminoácidos y, por lo tanto, establece el tipo de secuencia (`--seqtype`) como *proteína*. 

| genome     | query_size     | aligned_start     | aligned_end        | matches | mismatches | %_aligned | ... |
| -------------- |-------------------------| ------------------------| -------------- | ----------|-----|---|---|
| taeGut2| 88 | 	12 | 88 | 77 | 0 | 87.5 | ... |

#### [Màs ejemplos](https://github.com/pachterlab/gget_examples)
