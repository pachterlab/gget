> Parámetros de Python són iguales a los parámetros largos (`--parámetro`) de Terminal, si no es especificado de otra manera. Las banderas son designadas como cierto o falso (True/False) en Python. El manuál para cualquier modulo de gget se puede obtener desde Terminal con la bandera `-h` `--help`.  
# gget virus 🦠

Descarga secuencias de virus y metadatos asociados desde la [base de datos NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/). `gget virus` aplica filtros tanto del lado del servidor como locales para descargar de forma eficiente conjuntos de datos personalizados.  
Regresa: archivos FASTA, CSV y JSONL guardados en una carpeta de salida.  

**Argumento posicional**

`virus`  
Nombre taxonómico del virus (por ejemplo, "Zika virus"), ID taxonómico (por ejemplo, 2697049) o número de acceso (por ejemplo, `NC_045512.2`).

Línea de comandos: `gget virus "Zika virus" ...`  
Python: `gget.virus("Zika virus", ...)`

**Argumentos opcionales**

### Argumentos de salida

`-o` `--out`  
Ruta a la carpeta donde se guardarán los resultados. Valor por defecto: directorio de trabajo actual.  
Python: `outfolder="ruta/a/carpeta"`

### Filtros de hospedador

`--host`  
Filtra por nombre del organismo hospedador o por ID de taxonomía NCBI (por ejemplo, `human`, `Aedes aegypti`, `1335626`).

### Filtros de secuencia y genes

`--nuc_completeness`  
Filtra por integridad de la secuencia nucleotídica. Valores posibles: `complete` o `partial`.

`--min_seq_length`  
Filtra por longitud mínima de la secuencia.

`--max_seq_length`  
Filtra por longitud máxima de la secuencia.

`--min_gene_count`  
Filtra por número mínimo de genes.

`--max_gene_count`  
Filtra por número máximo de genes.

`--min_protein_count`  
Filtra por número mínimo de proteínas.

`--max_protein_count`  
Filtra por número máximo de proteínas.

`--min_mature_peptide_count`  
Filtra por número mínimo de péptidos maduros.

`--max_mature_peptide_count`  
Filtra por número máximo de péptidos maduros.

`--max_ambiguous_chars`  
Filtra por el número máximo de caracteres nucleotídicos ambiguos (N).

`--has_proteins`  
Filtra por secuencias que contengan proteínas o genes específicos (por ejemplo, `spike`, `ORF1ab`). Puede ser un único nombre de proteína o una lista de nombres de proteínas.  
Python: `has_proteins="spike"` o `has_proteins=["spike", "ORF1ab"]`

### Filtros de fechas

`--min_collection_date`  
Filtra por la fecha mínima de toma de muestra (YYYY-MM-DD).

`--max_collection_date`  
Filtra por la fecha máxima de toma de muestra (YYYY-MM-DD).

`--min_release_date`  
Filtra por la fecha mínima de publicación de la secuencia (YYYY-MM-DD).

`--max_release_date`  
Filtra por la fecha máxima de publicación de la secuencia (YYYY-MM-DD).

### Filtros de localización y remitente

`--geographic_location`  
Filtra por la localización geográfica de la toma de muestra (por ejemplo, `USA`, `Asia`).

`--submitter_country`  
Filtra por el país del remitente de la secuencia.

`--source_database`  
Filtra por la base de datos de origen. Valores posibles: `genbank` o `refseq`.

### Filtros específicos de SARS-CoV-2

`--lineage`  
Filtra por linaje de SARS-CoV-2 (por ejemplo, `B.1.1.7`, `P.1`).

**Flags**

`-a` `--is_accession`  
Indica que el argumento posicional `virus` es un número de acceso.

`--refseq_only`  
Restringe la búsqueda únicamente a genomas RefSeq (secuencias curadas de mayor calidad).

`--is_sars_cov2`  
Usa los paquetes de datos en caché optimizados de NCBI para una consulta de SARS-CoV-2. Esto permite descargas más rápidas y fiables. El sistema puede detectar automáticamente consultas de SARS-CoV-2 por nombre taxonómico, pero para consultas basadas en número de acceso se debe activar explícitamente esta opción.

`--is_alphainfluenza`  
Usa los paquetes de datos en caché optimizados de NCBI para consultas de Alphainfluenza (virus de la gripe A). Esto permite descargas más rápidas y fiables para conjuntos grandes de datos de gripe A. El sistema puede detectar automáticamente consultas de Alphainfluenza por nombre taxonómico, pero para consultas basadas en número de acceso se debe activar explícitamente esta opción.

`--genbank_metadata`  
Recupera y guarda metadatos detallados adicionales de GenBank, incluyendo fechas de toma de muestra, detalles del hospedador y referencias de publicación, en un archivo separado `{virus}_genbank_metadata.csv` (además de volcados completos en XML/CSV).

`--genbank_batch_size`  
Tamaño de lote para las peticiones de metadatos de GenBank. Valor por defecto: 200. Lotes más grandes pueden ser más rápidos pero también más propensos a *timeouts*.  
Python: `genbank_batch_size=200`

`--annotated`  
Filtra para devolver solo secuencias que tengan anotación de genes/proteínas.  
Línea de comandos: `--annotated true` o `--annotated false`.  
Python: `annotated=True` o `annotated=False`.

`--lab_passaged`  
Filtra a favor o en contra de muestras pasadas por laboratorio.  
Línea de comandos: `--lab_passaged true` para obtener solo muestras pasadas por laboratorio, o `--lab_passaged false` para excluirlas.  
Python: `lab_passaged=True` o `lab_passaged=False`.

`--proteins_complete`  
Incluye solo secuencias en las que todas las proteínas anotadas estén completas.

`--keep_temp`  
Conserva todos los archivos intermedios/temporales generados durante el procesamiento. Por defecto solo se mantienen los archivos finales.

`--download_all_accessions`  
⚠️ **ADVERTENCIA**: descarga **TODOS** los accesos de virus de NCBI (taxonomía completa Viruses, ID taxonómico 10239). Es un conjunto de datos extremadamente grande que puede tardar muchas horas en descargarse y requerir mucho espacio en disco. Úsalo con precaución y asegúrate de tener suficiente almacenamiento y ancho de banda. Cuando esta opción está activa, se ignora el argumento `virus`.

`-q` `--quiet`  
Solo línea de comandos. Evita que se muestre información de progreso.  
Python: usa `verbose=False` para ocultar la información de progreso.

### Ejemplo básico

```bash
gget virus "Zika virus" --nuc_completeness complete --host human --out zika_data
```

```python
# Python
import gget

gget.virus(
	"Zika virus",
	nuc_completeness="complete",
	host="human",
	outfolder="zika_data"
)
```

&rarr; Descarga genomas completos del virus del Zika de hospedadores humanos. Los resultados se guardan en la carpeta `zika_data` como `Zika_virus_sequences.fasta`, `Zika_virus_metadata.csv`, `Zika_virus_metadata.jsonl` y `command_summary.txt`.

El archivo CSV de metadatos tendrá un aspecto similar a:

| accession | Organism Name | GenBank/RefSeq | Release date | Length | Nuc Completeness | Geographic Location | Host | ... |
|---|---|---|---|---|---|---|---|---|
| KX198135.1 | Zika virus | GenBank | 2016-05-18 | 10807 | complete | Americas:Haiti | Homo sapiens | ... |
| . . . | . . . | . . . | . . . | . . . | . . . | . . . | . . . | ... |

El archivo de resumen de comando (`command_summary.txt`) contendrá, por ejemplo:

```
================================================================================
GGET VIRUS COMMAND SUMMARY
================================================================================

Execution Date: 2025-12-15 13:33:39
Output Folder: zika_data

--------------------------------------------------------------------------------
COMMAND LINE
--------------------------------------------------------------------------------
gget virus "Zika virus" --nuc_completeness complete --host human --out zika_data

--------------------------------------------------------------------------------
EXECUTION STATUS
--------------------------------------------------------------------------------
✓ Command completed successfully

--------------------------------------------------------------------------------
SEQUENCE STATISTICS
--------------------------------------------------------------------------------
Total records from API: 234
After metadata filtering: 234
Final sequences (after all filters): 234

--------------------------------------------------------------------------------
DETAILED STATISTICS
--------------------------------------------------------------------------------
Unique hosts: 1
	- Homo sapiens

Unique geographic locations: 15
	- Americas:Brazil
	- Americas:Colombia
	- ... (showing top 20)

Sequence length range: 10272 - 11155 bp
Average sequence length: 10742 bp

Completeness breakdown:
	- complete: 234

Source database breakdown:
	- GenBank: 233
	- RefSeq: 1

Unique submitter countries: 12
	- USA
	- Brazil
	- ... (showing top 20)

--------------------------------------------------------------------------------
OUTPUT FILES
--------------------------------------------------------------------------------
FASTA Sequences: Zika_virus_sequences.fasta (2.45 MB)
JSONL Metadata: Zika_virus_metadata.jsonl (0.53 MB)
CSV Metadata: Zika_virus_metadata.csv (0.42 MB)

================================================================================
END OF SUMMARY
================================================================================
```

**Nota**: si alguna operación falla durante la ejecución (timeouts de la API, errores al descargar secuencias, fallos al recuperar metadatos de GenBank), el resumen incluirá una sección "FAILED OPERATIONS - RETRY COMMANDS" con comandos y URLs exactos que se pueden ejecutar manualmente para reintentar las operaciones fallidas.

<br><br>
**Descargar un genoma de referencia específico de SARS-CoV-2 usando su número de acceso:**

```bash
gget virus NC_045512.2 --is_accession --is_sars_cov2
```

```python
# Python
import gget

gget.virus("NC_045512.2", is_accession=True, is_sars_cov2=True)
```

&rarr; Utiliza el método de descarga optimizado para SARS-CoV-2 para obtener el genoma de referencia y sus metadatos.

<br><br>
**Descargar secuencias de SARS-CoV-2 con optimización en caché Y metadatos de GenBank:**

```bash
gget virus "SARS-CoV-2" --host human --nuc_completeness complete --min_seq_length 29000 --genbank_metadata
```

```python
# Python
import gget

gget.virus(
	"SARS-CoV-2", 
	host="human", 
	nuc_completeness="complete",
	min_seq_length=29000,
	genbank_metadata=True,
	is_sars_cov2=True,
	outfolder="covid_data"
)
```

&rarr; Usa la descarga en caché para mayor velocidad (mediante los paquetes de datos de SARS-CoV-2 de NCBI cuando están disponibles), aplica el filtro de longitud de secuencia tras la descarga y recupera metadatos detallados de GenBank para todas las secuencias filtradas.

<br><br>
**Descargar secuencias de virus de la gripe A con caché optimizada y filtrado posterior a la descarga:**

```bash
gget virus "Influenza A virus" --host human --nuc_completeness complete --max_seq_length 15000 --genbank_metadata --is_alphainfluenza
```

```python
# Python
import gget

gget.virus(
	"Influenza A virus", 
	host="human", 
	nuc_completeness="complete",
	max_seq_length=15000,
	genbank_metadata=True,
	is_alphainfluenza=True,
	outfolder="influenza_a_data"
)
```

&rarr; Utiliza los paquetes de datos en caché de NCBI para Alphainfluenza para descargar genomas completos de gripe A en humanos mucho más rápido que con el método estándar de la API, y después aplica el filtro de longitud de secuencia y recupera metadatos de GenBank.

#### [Más ejemplos](https://github.com/pachterlab/gget_examples)

# Citar    
Si utiliza `gget virus` en una publicación, favor de citar los siguientes artículos:

- Luebbert, L., & Pachter, L. (2023). Efficient querying of genomic reference databases with gget. Bioinformatics. [https://doi.org/10.1093/bioinformatics/btac836](https://doi.org/10.1093/bioinformatics/btac836)

- O’Leary, N.A., Cox, E., Holmes, J.B. et al (2024). Exploring and retrieving sequence and metadata for species across the tree of life with NCBI Datasets. Sci Data 11, 732. [https://doi.org/10.1038/s41597-024-03571-y](https://doi.org/10.1038/s41597-024-03571-y)

