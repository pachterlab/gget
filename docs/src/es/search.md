> Parámetros de Python són iguales a los parámetros largos (`--parámetro`) de Terminal, si no especificado de otra manera. Banderas son parámetros de verdadero o falso (True/False) en Python. El manuál para cualquier modulo de gget se puede llamar desde la Terminal con la bandera `-h` `--help`.  
## gget search 🔎
Obtenga genes y transcripciones de [Ensembl](https://www.ensembl.org/) usando términos de búsqueda de forma libre.     
Los resultados se comparan según las secciones "nombre del gen" y "descripción" en la base de datos de Ensembl. `gget` versión >= 0.27.9 también incluye resultados que coinciden con la sección "sinónimo" de Ensembl.    
Regresa: Resultados en formato JSON (Terminal) o Dataframe/CSV (Python).  

**Parámetro posicional**  
`searchwords`   
Una o más palabras de búsqueda de forma libre, p. ej. gaba nmda. (Nota: la búsqueda no distingue entre mayúsculas y minúsculas).  

**Otros parámetros requeridos**   
`-s` `--species`  
Especies o base de datos a buscar.   
Una especie se puede pasar en el formato 'género_especie', p. ej. 'homo_sapiens' o 'arabidopsis_thaliana'.  
Para pasar una base de datos específica, pase el nombre de la base de datos CORE, p. ej. 'mus_musculus_dba2j_core_105_1'.  
  
Todas las bases de datos disponibles para cada versión de Ensembl se pueden encontrar aquí:  
Vertebrados: [http://ftp.ensembl.org/pub/current/mysql/](http://ftp.ensembl.org/pub/current/mysql/)  
Invertebrados: [http://ftp.ensemblgenomes.org/pub/current/](http://ftp.ensemblgenomes.org/pub/current/) + selecciona reino animal + selecciona mysql/  
  
Accesos directos: 'human', 'mouse'  

**Parámetros optionales**  
`-r` `--release`   
Define el número de versión de Ensembl desde el que se obtienen los archivos, p. ej. 104. Por defecto: None -> se usa la última versión de Ensembl.  
Nota: *No se aplica a las especies invertebrados* (en su lugar, puede pasar una base de datos de una especies específica (incluyen un número de versión) al argumento `species`).    
Este argumento se sobrescribe si se pasa una base de datos específica (que incluye un número de publicación) al argumento `species`.  

`-t` `--id_type`  
'gene' (esto se use por defecto) o 'transcript'   
Regesa genes o transcripciones, respectivamente.  

`-ao` `--andor`  
'or' (esto se use por defecto) o 'and'  
'or' ('o'): Regresa todos los genes que INCLUYEN AL MENOS UNA de las palabras de búsqueda en su nombre/descripción.  
'and' ('y'): Regresa solo los genes que INCLUYEN TODAS las palabras de búsqueda en su nombre/descripción.  

`-l` `--limit`   
Limita el número de resultados de búsqueda, p. ej. 10. Por defecto: None.  

`-o` `--out`   
Ruta al archivo en el que se guardarán los resultados, p. ej. ruta/al/directorio/resultados.csv (o .json). Por defecto: salida estándar (STDOUT).  
Para Python, usa `save=True` para guardar los resultados en el directorio de trabajo actual.  

**Banderas**  
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
gget search -s human gaba gamma-aminobutyric
```
```python
# Python
gget.search(["gaba", "gamma-aminobutyric"], "homo_sapiens")
```
&rarr; Regresa todos los genes que contienen al menos una de las palabras de búsqueda en su nombre o descripción de Ensembl/referencia externa:  

| ensembl_id     | gene_name     | ensembl_description     | ext_ref_description        | biotype | url |
| -------------- |-------------------------| ------------------------| -------------- | ----------|-----|
| ENSG00000034713| GABARAPL2 | 	GABA type A receptor associated protein like 2 [Source:HGNC Symbol;Acc:HGNC:13291] | GABA type A receptor associated protein like 2 | protein_coding | https://uswest.ensembl.org/homo_sapiens/Gene/Summary?g=ENSG00000034713 |
| . . .            | . . .                     | . . .                     | . . .            | . . .       | . . . |
    
#### [Más ejemplos](https://github.com/pachterlab/gget_examples)
