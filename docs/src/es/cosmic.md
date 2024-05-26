> Parámetros de Python són iguales a los parámetros largos (`--parámetro`) de Terminal, si no especificado de otra manera. Las banderas son parámetros de verdadero o falso (True/False) en Python. El manuál para cualquier modulo de gget se puede llamar desde la Terminal con la bandera `-h` `--help`.   
## gget cosmic 🪐
Busque genes, mutaciones, etc. asociados con cánceres utilizando la base de datos [COSMIC](https://cancer.sanger.ac.uk/cosmic) (Catálogo de mutaciones somáticas en cáncer).  
Produce: Resultados en formato JSON (Terminal) o Dataframe/CSV (Python).  
`gget cosmic` fue escrito por [@AubakirovArman](https://github.com/AubakirovArman).

Se aplican tarifas de licencia para el uso comercial de COSMIC. Puede leer más sobre la concesión de licencias de datos COSMIC [aquí](https://cancer.sanger.ac.uk/cosmic/license).

**Parámetro posicional**  
`searchterm`   
Término de búsqueda. Puede ser una mutación, un nombre de gen (o ID de Ensembl), tipo de cáncer, sitio del tumor, ID de estudio, ID de PubMed o ID de muestra, tal como se define con el argumento `entity`. Ejemplo: 'EGFR'  

**Parámetros optionales**  
`-e` `--entity`  
'mutations' (mutación), 'genes' (nombre de gen / ID de Ensembl), 'cancer' (tipo de cáncer), 'tumour site' (sitio del tumor), 'studies' (ID de estudio), 'pubmed' (ID de PubMed), o 'samples' (ID de muestra). Por defecto: 'mutations'.  
Define el tipo de término de búsqueda (`searchterm`).  

`-l` `--limit`  
Limita el número de resultados producidos. Por defecto: 100.  

`-o` `--out`   
Ruta al archivo en el que se guardarán los resultados, p. ej. ruta/al/directorio/resultados.csv (o .json). Por defecto: salida estándar (STDOUT).  
Para Python, usa `save=True` para guardar los resultados en el directorio de trabajo actual.  

**Banderas**  
`-csv` `--csv`
Solo para Terminal. Produce los resultados en formato CSV.  
Para Python, usa json=True para producir los resultados en formato JSON.

`-q` `--quiet`   
Solo para Terminal. Impide la información de progreso de ser exhibida durante la ejecución del programa.  
Para Python, usa `verbose=False` para imipidir la informacion de progreso de ser exhibida durante la ejecución del programa.  

  
### Por ejemplo    
```bash
gget cosmic -e genes EGFR
```
```python
# Python
gget.cosmic("EGFR", entity="genes")
```
&rarr; Produce los resultados COSMIC para el gen 'EGFR':  

| Gene     | Alternate IDs     | Tested samples     | Simple Mutations        | Fusions | Coding Mutations | ... |
| -------------- |-------------------------| ------------------------| -------------- | ----------|-----|---|
| EGFR| EGFR,ENST00000275493.6,... | 210280 | 31900 | 0 | 31900 | ... |
| . . . | . . . | . . . | . . . | . . . | . . . | . . . | ... | 
