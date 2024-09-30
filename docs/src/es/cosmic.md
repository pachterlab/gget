> Parámetros de Python són iguales a los parámetros largos (`--parámetro`) de Terminal, si no especificado de otra manera. Las banderas son parámetros de verdadero o falso (True/False) en Python. El manuál para cualquier modulo de gget se puede llamar desde la Terminal con la bandera `-h` `--help`.   
## gget cosmic 🪐
Busque genes, mutaciones, etc. asociados con cánceres utilizando la base de datos [COSMIC](https://cancer.sanger.ac.uk/cosmic) (Catálogo de mutaciones somáticas en cáncer).  
Produce: JSON (línea de comandos) o marco de datos/CSV (Python) cuando `download_cosmic=False`. Cuando `download_cosmic=True`, descarga la base de datos solicitada en la carpeta especificada.  

Este módulo fue escrito en parte por [@AubakirovArman](https://github.com/AubakirovArman) (consulta de información) y [@josephrich98](https://github.com/josephrich98) (descarga de base de datos).  

NOTA: Se aplican tarifas de licencia para el uso comercial de COSMIC. Puede leer más sobre la concesión de licencias de datos COSMIC [aquí](https://cancer.sanger.ac.uk/cosmic/license).  

**Parámetro posicional (para consultar información)**  
`searchterm`   
Término de búsqueda, que puede ser una mutación, un nombre de gen (o ID de Ensembl), una muestra, etc.  
Ejemplos para los argumentos de searchterm y entidad:   

| searchterm   | entidad    | |
|--------------|------------|-|
| EGFR         | mutaciones | -> Encuentra mutaciones en el gen EGFR asociadas con el cáncer |
| v600e        | mutaciones | -> Encuentra genes para los cuales una mutación v600e está asociada con el cáncer |
| COSV57014428 | mutaciones | -> Encuentra mutaciones asociadas con esta ID de mutaciones COSMIC |
| EGFR         | genes      | -> Obtiene el número de muestras, mutaciones simples/codificantes y fusiones observadas en COSMIC para EGFR |
| prostate     | cáncer     | -> Obtiene el número de muestras probadas y mutaciones para el cáncer de próstata |
| prostate     | sitio_tumoral | -> Obtiene el número de muestras probadas, genes, mutaciones, fusiones, etc. con 'próstata' como sitio de tejido primario |
| ICGC         | estudios   | -> Obtiene el código de proyecto y descripciones de todos los estudios de ICGC (Consortio Internacional del Genoma del Cáncer) |
| EGFR         | pubmed     | -> Encuentra publicaciones de PubMed sobre EGFR y cáncer |
| ICGC         | muestras   | -> Obtiene metadatos sobre todas las muestras de ICGC (Consortio Internacional del Genoma del Cáncer) |
| COSS2907494  | muestras   | -> Obtiene metadatos sobre esta ID de muestra COSMIC (tipo de cáncer, tejido, # genes analizados, # mutaciones, etc.) |

NOTA: (Solo Python) Establezca en `None` cuando se descarguen bases de datos COSMIC con `download_cosmic=True`.  

**Parámetros opcionales (para consultar información)**  
`-e` `--entity`  
'mutations' (predeterminado), 'genes', 'cáncer', 'sitio_tumoral', 'estudios', 'pubmed' o 'muestras'.
Define el tipo de resultados a devolver.

`-l` `--limit`  
Limita el número de resultados a devolver. Predeterminado: 100.

**Banderas (para descargar bases de datos COSMIC)**  
`-d` `--download_cosmic`  
Conmuta al modo de descarga de base de datos.

`-gm` `--gget_mutate`  
DESACTIVA la creación de una versión modificada de la base de datos para usar con gget mutate.
Python: `gget_mutate` es Verdadero por defecto. Establezca `gget_mutate=False` para deshabilitar.

**Parámetros opcionales (para descargar bases de datos COSMIC)**  
`-mc` `--mutation_class`  
'cancer' (predeterminado), 'cell_line', 'census', 'resistance', 'genome_screen', 'targeted_screen', o 'cancer_example'  
Tipo de base de datos COSMIC para descargar.

`-cv` `--cosmic_version`  
Versión de la base de datos COSMIC. Predeterminado: Ninguno -> Se establece en la última versión por defecto.

`-gv` `--grch_version`  
Versión del genoma de referencia humano GRCh en el que se basó la base de datos COSMIC (37 o 38). Predeterminado: 37

**Parámetros opcionales (generales)**  
`-o` `--out`   
Ruta al archivo (o carpeta cuando se descargan bases de datos con el flag `download_cosmic`) donde se guardarán los resultados, p. ej. 'ruta/a/resultados.json'.  
Predeterminado: None  
-> Cuando download_cosmic=False: Los resultados se devolverán a la salida estándar  
-> Cuando download_cosmic=True: La base de datos se descargará en el directorio de trabajo actual  

**Banderas (generales)**  
`-csv` `--csv`  
Solo para Terminal. Produce los resultados en formato CSV.  
Para Python, usa json=True para producir los resultados en formato JSON.

`-q` `--quiet`   
Solo para Terminal. Impide la información de progreso de ser exhibida durante la ejecución del programa.  
Para Python, usa `verbose=False` para imipidir la informacion de progreso de ser exhibida durante la ejecución del programa.  

  
### Por ejemplo    
#### Consultar información
```bash
gget cosmic -e genes EGFR
```
```python
# Python
gget.cosmic("EGFR", entity="genes")
```
&rarr; Devuelve mutaciones en el gen EGFR asociadas con el cáncer en el formato:

| Gene     | Syntax     | Alternate IDs                  | Canonical  |
| -------- |------------| -------------------------------| ---------- |
| EGFR     | c.*2446A>G | EGFR c.*2446A>G, EGFR p.?, ... | y          |
| EGFR     | c.(2185_2283)ins(18) | EGFR c.(2185_2283)ins(18), EGFR p.?, ... | y          |
| . . .    | . . .      | . . .                          | . . .      | 

### Descargar bases de datos COSMIC
```bash
gget cosmic --download_cosmic
```
```python
# Python
gget.cosmic(searchterm=None, download_cosmic=True)  
```
&rarr; Descargue la base de datos sobre cáncer de COSMIC de la última versión de COSMIC.  

# Citar    
Si utiliza `gget alphafold` en una publicación, favor de citar los siguientes artículos:

- Luebbert, L., & Pachter, L. (2023). Efficient querying of genomic reference databases with gget. Bioinformatics. [https://doi.org/10.1093/bioinformatics/btac836](https://doi.org/10.1093/bioinformatics/btac836)

- Tate JG, Bamford S, Jubb HC, Sondka Z, Beare DM, Bindal N, Boutselakis H, Cole CG, Creatore C, Dawson E, Fish P, Harsha B, Hathaway C, Jupe SC, Kok CY, Noble K, Ponting L, Ramshaw CC, Rye CE, Speedy HE, Stefancsik R, Thompson SL, Wang S, Ward S, Campbell PJ, Forbes SA. COSMIC: the Catalogue Of Somatic Mutations In Cancer. Nucleic Acids Res. 2019 Jan 8;47(D1):D941-D947. doi: [10.1093/nar/gky1015](https://doi.org/10.1093/nar/gky1015). PMID: 30371878; PMCID: PMC6323903.
