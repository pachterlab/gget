[<kbd> Ver el codigo fuente de la pagina en GitHub </kbd>](https://github.com/pachterlab/gget/blob/main/docs/src/es/cosmic.md)

> Parámetros de Python són iguales a los parámetros largos (`--parámetro`) de Terminal, si no especificado de otra manera. Las banderas son parámetros de verdadero o falso (True/False) en Python. El manuál para cualquier modulo de gget se puede llamar desde la Terminal con la bandera `-h` `--help`.   
# gget cosmic 🪐  
Busca genes, mutaciones y otros factores asociados con el cáncer utilizando la base de datos [COSMIC](https://cancer.sanger.ac.uk/cosmic) (Catalogue Of Somatic Mutations In Cancer).  
Formato de retorno: JSON (línea de comandos) o data frame/CSV (Python) cuando `download_cosmic=False`. Cuando `download_cosmic=True`, se descarga la base de datos solicitada en la carpeta especificada.    

Este módulo fue escrito originalmente en parte por [@AubakirovArman](https://github.com/AubakirovArman) (consultas de información) y [@josephrich98](https://github.com/josephrich98) (descarga de bases de datos).  

NOTA: Se aplican tarifas de licencia para el uso comercial de COSMIC. Puedes leer más sobre la licencia de los datos de COSMIC [aquí](https://cancer.sanger.ac.uk/cosmic/license).  

NOTA: Al utilizar este módulo por primera vez, primero descarga una base de datos de COSMIC para obtener `cosmic_tsv_path` (ver ejemplos abajo).  

**Argumento posicional (para consultar información)**  
`searchterm`  
Término de búsqueda, que puede ser una mutación, un nombre de gen (o ID de Ensembl), una muestra, etc.  
Ejemplos: 'EGFR', 'ENST00000275493', 'c.650A>T', 'p.Q217L', 'COSV51765119', 'BT2012100223LNCTB' (ID de muestra)  
NOTA: (solo en Python) Establecer en `None` al descargar bases de datos de COSMIC con `download_cosmic=True`.  

**Argumento obligatorio (para consultar información)**  
`-ctp` `--cosmic_tsv_path`  
Ruta al archivo tsv de la base de datos de COSMIC, por ejemplo: 'path/to/CancerMutationCensus_AllData_v101_GRCh37.tsv'.  
Este archivo se descarga al usar los argumentos descritos debajo para descargar bases de datos.  
NOTA: Este argumento es obligatorio cuando `download_cosmic=False`.  

**Argumentos opcionales (para consultar información)**  
`-l` `--limit`  
Límite en la cantidad de resultados a devolver. Valor por defecto: 100.  

`-csv` `--csv`  
Solo para línea de comandos. Devuelve los resultados en formato CSV.  
Python: usa `json=True` para obtener la salida en formato JSON.  

**Banderas (para descargar bases de datos de COSMIC)**  
`-d` `--download_cosmic`  
Activa el modo de descarga de base de datos.  

`-gm` `--gget_mutate`  
Crea una versión modificada de la base de datos COSMIC para usar con [`gget mutate`](mutate.md).  

**Argumentos opcionales (para descargar bases de datos de COSMIC)**  
`-cp` `--cosmic_project`  
'cancer' (por defecto), 'cancer_example', 'census', 'resistance', 'cell_line', 'genome_screen', o 'targeted_screen'  
Tipo de base de datos COSMIC a descargar:  

| cosmic_project  | Descripción                                                            | Notas                                                                               | Tamaño  |
|-----------------|------------------------------------------------------------------------|-------------------------------------------------------------------------------------|---------|
| cancer          | Cancer Mutation Census (CMC) (conjunto más comúnmente usado de COSMIC) | Solo disponible para GRCh37. Esquema más completo (requiere más tiempo para buscar). | 2 GB    |
| cancer_example  | Subconjunto de CMC de ejemplo para pruebas y demostración              | Descargable sin cuenta COSMIC. Conjunto de datos mínimo.                           | 2.5 MB  |
| census          | Censo de mutaciones somáticas en genes conocidos relacionado al cáncer | Conjunto curado más pequeño de genes impulsores del cáncer.                        | 630 MB  |
| resistance      | Mutaciones asociadas con resistencia a fármacos                       | Útil para investigación en farmacogenómica.                                        | 1.6 MB  |
| cell_line       | Datos de mutaciones del proyecto de líneas celulares                  | A menudo incluye metadatos de muestras.                                            | 2.7 GB  |
| genome_screen   | Mutaciones de estudios de cribado genómico                            | Incluye datos menos curados, útiles para estudios a gran escala.                   |         |
| targeted_screen | Mutaciones de paneles de cribado dirigido                             | Datos centrados, útiles en contextos clínicos.                                     |         |

`-cv` `--cosmic_version`  
Versión de la base de datos COSMIC. Valor por defecto: None → se usa la versión más reciente.  

`-gv` `--grch_version`  
Versión del genoma de referencia humano GRCh en el que se basa la base de datos COSMIC (37 o 38). Por defecto: 37  

`--email`  
Correo electrónico para iniciar sesión en COSMIC. Útil para evitar la entrada manual al ejecutar gget COSMIC. Por defecto: None  

`--password`  
Contraseña para iniciar sesión en COSMIC. Útil para evitar la entrada manual al ejecutar gget COSMIC, pero se almacenará en texto plano en el script. Por defecto: None  

**Argumentos adicionales para la bandera `--gget_mutate`**  
`--keep_genome_info`  
Indica si se debe conservar la información genómica en la base modificada para usar con `gget mutate`. Por defecto: False  

`--remove_duplicates`  
Indica si se deben eliminar filas duplicadas de la base modificada para usar con `gget mutate`. Por defecto: False  

`--seq_id_column`  
(str) Nombre de la columna de ID de secuencia en el archivo CSV creado por `gget_mutate`. Por defecto: "seq_ID"  

`--mutation_column`  
(str) Nombre de la columna de mutaciones en el archivo CSV creado por `gget_mutate`. Por defecto: "mutation"  

`--mut_id_column`  
(str) Nombre de la columna de ID de mutación en el archivo CSV creado por `gget_mutate`. Por defecto: "mutation_id"  

**Argumentos opcionales (generales)**  
`-o` `--out`  
Ruta del archivo (o carpeta cuando se descargan bases de datos con la bandera `download_cosmic`) donde se guardarán los resultados, por ejemplo: 'path/to/results.json'.  
Valores por defecto:  
→ Cuando `download_cosmic=False`: los resultados se devuelven por la salida estándar  
→ Cuando `download_cosmic=True`: la base de datos se descarga en el directorio de trabajo actual  

**Banderas (generales)**  
`-q` `--quiet`  
Solo en línea de comandos. Evita que se muestren mensajes de progreso.  
Python: usa `verbose=False` para evitar mensajes de progreso.  

---

### Ejemplos
#### Descargar la base de datos "cancer" de COSMIC y consultar información
```bash
# The download_cosmic command will ask for your COSMIC email and password and only needs to be run once
gget cosmic --download_cosmic --cosmic_project cancer

gget cosmic EGFR --cosmic_tsv_path 'CancerMutationCensus_AllData_Tsv_v101_GRCh37/CancerMutationCensus_AllData_v101_GRCh37.tsv'
```
```python
# Python
# The download_cosmic command will ask for your COSMIC email and password and only needs to be run once
gget.cosmic(searchterm=None, download_cosmic=True, cosmic_project="cancer")

gget.cosmic("EGFR", cosmic_tsv_path="CancerMutationCensus_AllData_Tsv_v101_GRCh37/CancerMutationCensus_AllData_v101_GRCh37.tsv")
```

&rarr; El primer comando descarga la base de datos solicitada de la última versión de COSMIC en el directorio de trabajo actual. El segundo comando busca en la base de datos las mutaciones asociadas al gen 'EGFR' y devuelve los resultados en el siguiente formato:  

| GENE_NAME | ACCESSION_NUMBER | ONC_TSG | Mutation_CDS | Mutation_AA |  ... |
| ---- | ---- | ---- | ---- | ---- | ---- |
| EGFR | ENST00000275493.2 | oncogene | c.650A>T | p.Q217L | ... |
| EGFR | ENST00000275493.2 | oncogene | c.966C>T | p.G322= | ... |
| ... | ... | ... | ... | ... | ... |

#### Descargar la base de datos "census" de COSMIC y consultar información
```bash
# The download_cosmic command will ask for your COSMIC email and password and only needs to be run once
gget cosmic --download_cosmic --cosmic_project census

gget cosmic EGFR --cosmic_tsv_path 'Cosmic_MutantCensus_Tsv_v101_GRCh37/Cosmic_MutantCensus_v101_GRCh37.tsv'
```
```python
# Python
# The download_cosmic command will ask for your COSMIC email and password and only needs to be run once
gget.cosmic(searchterm=None, download_cosmic=True, cosmic_project="cancer")

gget.cosmic("EGFR", cosmic_tsv_path="Cosmic_MutantCensus_Tsv_v101_GRCh37/Cosmic_MutantCensus_v101_GRCh37.tsv")
```

&rarr; El primer comando descarga la base de datos solicitada de la última versión de COSMIC en el directorio de trabajo actual. El segundo comando busca en la base de datos las mutaciones asociadas al gen 'EGFR' y devuelve los resultados en el siguiente formato:  

| GENE_SYMBOL | COSMIC_GENE_ID | MUTATION_DESCRIPTION | MUTATION_CDS | Mutation_AA | MUTATION_SOMATIC_STATUS | ... |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| EGFR | COSG35617 | inframe_deletion | c.2235_2249del | 	p.E746_A750del | Reported in another cancer sample as somatic | ... |
| EGFR | COSG35617 | missense_variant | c.2573T>G | p.L858R | Reported in another cancer sample as somatic | ... |
| ... | ... | ... | ... | ... | ... | ... |


# Citar    
Si utiliza `gget cosmic` en una publicación, favor de citar los siguientes artículos:

- Luebbert, L., & Pachter, L. (2023). Efficient querying of genomic reference databases with gget. Bioinformatics. [https://doi.org/10.1093/bioinformatics/btac836](https://doi.org/10.1093/bioinformatics/btac836)

- Tate JG, Bamford S, Jubb HC, Sondka Z, Beare DM, Bindal N, Boutselakis H, Cole CG, Creatore C, Dawson E, Fish P, Harsha B, Hathaway C, Jupe SC, Kok CY, Noble K, Ponting L, Ramshaw CC, Rye CE, Speedy HE, Stefancsik R, Thompson SL, Wang S, Ward S, Campbell PJ, Forbes SA. COSMIC: the Catalogue Of Somatic Mutations In Cancer. Nucleic Acids Res. 2019 Jan 8;47(D1):D941-D947. doi: [10.1093/nar/gky1015](https://doi.org/10.1093/nar/gky1015). PMID: 30371878; PMCID: PMC6323903.




