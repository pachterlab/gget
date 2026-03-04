[<kbd> Ver el codigo fuente de la pagina en GitHub </kbd>](https://github.com/pachterlab/gget/blob/main/docs/src/es/virus.md)

> Parámetros de Python són iguales a los parámetros largos (`--parámetro`) de Terminal, si no especificado de otra manera. Banderas son parámetros de verdadero o falso (True/False) en Python. El manuál para cualquier modulo de gget se puede llamar desde la Terminal con la bandera `-h` `--help`.  
# gget virus 🦠  

Descargue secuencias nucleotídicas virales, junto con metadatos ricos y vinculados, de toda la International Nucleotide Sequence Database Collaboration ([INSDC](https://www.insdc.org/)), incluyendo NCBI, [ENA](https://www.ebi.ac.uk/ena/browser/) y [DDBJ](https://www.ddbj.nig.ac.jp/index-e.html) (a través de [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/)), con la opción de enriquecer adicionalmente los resultados usando metadatos de NCBI GenBank (por ejemplo, anotaciones de genes y proteínas, secuencias de aminoácidos y más). `gget virus` aplica filtros secuenciales tanto del lado del servidor como locales para descargar de forma eficiente conjuntos de datos personalizados.

Formato de salida: archivos FASTA, CSV y JSONL guardados en una carpeta de salida.  

[Cuaderno de Google Colab sin código y compartible para descargar secuencias virales.](https://colab.research.google.com/github/pachterlab/gget_examples/blob/main/gget_virus/gget_virus_colab.ipynb)

Este módulo fue escrito por [Ferdous Nasri](https://github.com/ferbsx).

**Nota**: Para consultas de SARS-CoV-2 y Alphainfluenza (Influenza A), `gget virus` utiliza los paquetes de datos optimizados en caché de NCBI mediante la [NCBI datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/). El binario de la CLI de datasets se incluye con gget para las principales plataformas—no se requiere instalación adicional. Si ya tienes la CLI `datasets` instalada en tu sistema, gget usará automáticamente tu instalación existente.

**Argumento posicional**  
`virus`
Nombre del taxón viral (p. ej., 'Zika virus'), ID de taxón (p. ej., '2697049'), número de acceso de NCBI (p. ej., 'NC_045512.2'), lista de números de acceso separados por espacios (p. ej., 'NC_045512.2 MN908947.3 MT020781.1') o ruta a un archivo de texto que contenga números de acceso (uno por línea) (p. ej., 'path/to/text.txt').

Añada `--is_accession` al proporcionar un número de acceso de NCBI. Añada `--is_sars_cov2` o `--is_alphainfluenza` para la descarga optimizada de secuencias de SARS-CoV-2 o Alphainfluenza, respectivamente.

Para descargas en caché de SARS-CoV-2 y Alphainfluenza, se admite:

* Acceso único: `NC_045512.2`
* Lista separada por espacios: `NC_045512.2 MN908947.3 MT020781.1`
* Ruta a archivo de texto: `accessions.txt` (un número de acceso por línea)

Use la opción `--download_all_accessions` para aplicar filtros sin buscar un virus específico.

**Argumentos opcionales**   

_Filtros de hospedador_  

`--host`  
Filtra por nombre del organismo hospedador o ID de Taxonomía de NCBI (p. ej. 'human', 'Aedes aegypti', `1335626`).

_Filtros de Secuencia y Gen_  

`--nuc_completeness`  
Filtrar por integridad del nucleótido. Una de las siguientes opciones: 'complete' o 'partial'.  
Establezca 'complete' para devolver únicamente secuencias de nucleótidos marcadas como completas; establezca 'partial' para devolver únicamente secuencias marcadas como parciales.

`--min_seq_length`  
Filtra por longitud mínima de secuencia.

`--max_seq_length`  
Filtra por longitud máxima de secuencia.

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
Filtra por número máximo de caracteres nucleotídicos ambiguos (N).

`--has_proteins`  
Filtra por secuencias que contengan proteínas o genes específicos (p. ej. 'spike', 'ORF1ab'). Puede ser un solo nombre de proteína o una lista de nombres de proteínas.
Python: `has_proteins="spike"` o `has_proteins=["spike", "ORF1ab"]`

`--annotated`  
`'true'` o `'false'`. Filtra por secuencias que han sido anotadas con información de genes/proteínas.  
Línea de comandos: `--annotated true` para obtener únicamente secuencias anotadas con información de genes/proteínas, o `--annotated false` para excluirlas.  
Python: `annotated=True` o `annotated=False` (`annotated=None` para no aplicar ningún filtro).  

`--lab_passaged`  
`'true'` o `'false'`. Filtra a favor o en contra de muestras pasadas por laboratorio (*lab-passaged*).  
Línea de comandos: `--lab_passaged true` para obtener únicamente muestras pasadas por laboratorio, o `--lab_passaged false` para excluirlas.  
Python: `lab_passaged=True` o `lab_passaged=False` (`lab_passaged=None` para no aplicar ningún filtro).

`--vaccine_strain`  
Filtra a favor o en contra de secuencias de cepas de vacunas.  
Línea de comandos: `--vaccine_strain true` para obtener solo cepas de vacunas, o `--vaccine_strain false` para excluirlas.  
Python: `vaccine_strain=True` o `vaccine_strain=False` (`vaccine_strain=None` para no aplicar ningún filtro).

`--segment`  
Filtra por secuencias con segmento(s) específico(s) (p. ej. 'HA', 'NA'). Puede ser un solo nombre de segmento o una lista de nombres de segmentos.
Python: `segment="HA"` o `segment=["HA", "NA", "PB1"]`

_Filtros de fecha_  

`--min_collection_date`  
Filtra por fecha mínima de recolección de la muestra (YYYY-MM-DD).

`--max_collection_date`  
Filtra por fecha máxima de recolección de la muestra (YYYY-MM-DD).

`--min_release_date`  
Filtra por fecha mínima de liberación de la secuencia (YYYY-MM-DD).

`--max_release_date`  
Filtra por fecha máxima de liberación de la secuencia (YYYY-MM-DD).

_Filtros de ubicación y remitente_

`--geographic_location`  
Filtra por ubicación geográfica de la recolección de la muestra (p. ej. 'USA', 'Asia').

`--submitter_country`  
Filtra por el país del remitente de la secuencia. Puede ser un solo país o una lista separada por comas.

`--source_database`  
Filtra por base de datos de origen. Uno de: 'genbank' o 'refseq'.

_Filtros específicos de SARS-CoV-2_

`--lineage`  
Filtra por linaje de SARS-CoV-2 (p. ej. 'B.1.1.7', 'P.1'). Puede ser un solo linaje o una lista de linajes.
Python: `lineage="B.1.1.7"` o `lineage=["B.1.1.7", "P.1"]`

_Configuración del pipeline_

`--genbank_batch_size`  
Tamaño de lote para solicitudes a la API de metadatos de GenBank. Por defecto: 200. Lotes más grandes son más rápidos pero pueden ser más propensos a timeouts.  

`-o` `--out`  
Ruta a la carpeta donde se guardarán los resultados. Por defecto: directorio de trabajo actual.  
Python: `outfolder="path/to/folder"`

**Banderas**  
`-a` `--is_accession`  
Bandera para indicar que el argumento posicional `virus` es un número de acceso, una lista de accesiones separadas por espacios, o una ruta a un archivo de texto que contiene números de acceso (uno por línea).

`--download_all_accessions`  
Use esta bandera al aplicar filtros sin buscar un virus específico (deje el argumento `virus` vacío).  
⚠️ **ADVERTENCIA**: Si no especifica filtros adicionales, esta bandera descargará TODAS las secuencias virales disponibles de NCBI (toda la taxonomía de Virus, taxon ID 10239). Este es un conjunto de datos extremadamente grande que puede tardar muchas horas en descargarse y requerir un espacio considerable en disco. Úsela con precaución y asegúrese de contar con suficiente almacenamiento y ancho de banda. Cuando esta bandera está activada, el argumento `virus` se ignora.

`--is_sars_cov2`  
Usa los paquetes de datos optimizados en caché de NCBI para una consulta de SARS-CoV-2. Esto proporciona descargas más rápidas y confiables. El sistema puede detectar automáticamente consultas por nombre de taxón de SARS-CoV-2, pero para consultas basadas en accesiones debes establecer esta bandera explícitamente.

`--is_alphainfluenza`  
Usa los paquetes de datos optimizados en caché de NCBI para una consulta de Alphainfluenza (virus de la Influenza A). Esto proporciona descargas más rápidas y confiables para grandes conjuntos de datos de Influenza A. El sistema puede detectar automáticamente consultas por nombre de taxón de Alphainfluenza, pero para consultas basadas en accesiones debes establecer esta bandera explícitamente.

`-g` `--genbank_metadata`  
Obtiene y guarda metadatos adicionales detallados desde GenBank, incluyendo fechas de recolección, detalles del hospedador y referencias de publicaciones, en un archivo separado `{virus}_genbank_metadata.csv` (además de volcados completos XML/CSV dumps).

`--proteins_complete`  
Bandera para incluir solo secuencias donde todas las proteínas anotadas estén completas.  

`-kt` `--keep_temp`  
Bandera para conservar todos los archivos intermedios/temporales generados durante el procesamiento. Por defecto, solo se conservan los archivos de salida finales.

`-q` `--quiet`  
Uso limitado para Terminal. Impide la información de progreso de ser exhibida durante la ejecución del programa.  
Para Python, usa `verbose=False`.  

### Ejemplo

```bash
gget virus "Zika virus" --nuc_completeness complete --host human --out zika_data
````

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

→ Descarga genomas completos de Zika virus de hospedadores humanos. Los resultados se guardan en la carpeta `zika_data` como `Zika_virus_sequences.fasta`, `Zika_virus_metadata.csv`, `Zika_virus_metadata.jsonl` y `command_summary.txt`.


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

→ Uses the optimized download method for SARS-CoV-2 to fetch the reference genome and its metadata.

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

→ Uses cached download for speed (via NCBI's SARS-CoV-2 data packages when available), applies the sequence length filter post-download, and fetches detailed GenBank metadata for all filtered sequences.

<br><br>
**Descargar secuencias del virus de la Influenza A con caché optimizada y filtrado posterior a la descarga:**

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

→ Uses NCBI's cached data packages for Alphainfluenza to download complete Influenza A genomes from human hosts much faster than the standard API method, then applies the sequence length filter and fetches GenBank metadata.

#### [Más ejemplos](https://github.com/pachterlab/gget_examples/tree/main/gget_virus)

# Citar

Si utilizas `gget virus` en una publicación, por favor cita los siguientes artículos:

  - Nasri, F. et al (2026). En preparación.

  - Luebbert, L., & Pachter, L. (2023). Efficient querying of genomic reference databases with gget. Bioinformatics. [https://doi.org/10.1093/bioinformatics/btac836](https://doi.org/10.1093/bioinformatics/btac836)

  - O'Leary, N.A., Cox, E., Holmes, J.B. et al (2024). Exploring and retrieving sequence and metadata for species across the tree of life with NCBI Datasets. Sci Data 11, 732. [https://doi.org/10.1038/s41597-024-03571-y](https://doi.org/10.1038/s41597-024-03571-y)
