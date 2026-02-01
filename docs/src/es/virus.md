> Parámetros de Python són iguales a los parámetros largos (`--parámetro`) de Terminal, si no especificado de otra manera. Banderas son parámetros de verdadero o falso (True/False) en Python. El manuál para cualquier modulo de gget se puede llamar desde la Terminal con la bandera `-h` `--help`.  
# gget virus 🦠  

Descargue secuencias nucleotídicas virales, junto con metadatos ricos y vinculados, de toda la International Nucleotide Sequence Database Collaboration ([INSDC](https://www.insdc.org/)), incluyendo NCBI, [ENA](https://www.ebi.ac.uk/ena/browser/) y [DDBJ](https://www.ddbj.nig.ac.jp/index-e.html) (a través de [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/)), con la opción de enriquecer adicionalmente los resultados usando metadatos de NCBI GenBank (por ejemplo, anotaciones de genes y proteínas, secuencias de aminoácidos y más). `gget virus` aplica filtros secuenciales tanto del lado del servidor como locales para descargar de forma eficiente conjuntos de datos personalizados.

Formato de salida: archivos FASTA, CSV y JSONL guardados en una carpeta de salida.  

Este módulo fue escrito por [Ferdous Nasri](https://github.com/ferbsx).

**Nota**: Para consultas de SARS-CoV-2 y Alphainfluenza (Influenza A), `gget virus` utiliza los paquetes de datos optimizados en caché de NCBI mediante la [NCBI datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/). El binario de la CLI de datasets se incluye con gget para las principales plataformas—no se requiere instalación adicional. Si ya tienes la CLI `datasets` instalada en tu sistema, gget usará automáticamente tu instalación existente.

**Argumento posicional**  
`virus`  
Nombre del taxón del virus (p. ej. 'Zika virus'), ID taxonómico (p. ej. 2697049) o número de acceso (p. ej. 'NC\_045512.2').  

**Argumentos opcionales**   
`-o` `--out`  
Ruta a la carpeta donde se guardarán los resultados. Por defecto: directorio de trabajo actual.
Python: `outfolder="path/to/folder"`

_Filtros de hospedador_  

`--host`  
Filtra por nombre del organismo hospedador o ID de Taxonomía de NCBI (p. ej. 'human', 'Aedes aegypti', `1335626`).

_Filtros de Secuencia y Gen_  

`--nuc_completeness`  
Filtra por completitud nucleotídica. Uno de: 'complete' o 'partial'.

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
Filtra por el país del remitente de la secuencia.

`--source_database`  
Filtra por base de datos de origen. Uno de: 'genbank' o 'refseq'.

_Filtros específicos de SARS-CoV-2_

`--lineage`  
Filtra por linaje de SARS-CoV-2 (p. ej. 'B.1.1.7', 'P.1').

**Banderas**  
`-a` `--is_accession`  
Bandera para indicar que el argumento posicional `virus` es un número de acceso.

`--refseq_only`  
Bandera para limitar la búsqueda solo a genomas RefSeq (secuencias de mayor calidad, curadas).

`--is_sars_cov2`  
Usa los paquetes de datos optimizados en caché de NCBI para una consulta de SARS-CoV-2. Esto proporciona descargas más rápidas y confiables. El sistema puede detectar automáticamente consultas por nombre de taxón de SARS-CoV-2, pero para consultas basadas en accesiones debes establecer esta bandera explícitamente.

`--is_alphainfluenza`  
Usa los paquetes de datos optimizados en caché de NCBI para una consulta de Alphainfluenza (virus de la Influenza A). Esto proporciona descargas más rápidas y confiables para grandes conjuntos de datos de Influenza A. El sistema puede detectar automáticamente consultas por nombre de taxón de Alphainfluenza, pero para consultas basadas en accesiones debes establecer esta bandera explícitamente.

`-g` `--genbank_metadata`  
Obtiene y guarda metadatos adicionales detallados desde GenBank, incluyendo fechas de recolección, detalles del hospedador y referencias de publicaciones, en un archivo separado `{virus}_genbank_metadata.csv` (además de volcados completos XML/CSV dumps).

`--genbank_batch_size`  
Tamaño de lote para solicitudes a la API de metadatos de GenBank. Por defecto: 200. Lotes más grandes son más rápidos pero pueden ser más propensos a timeouts.  

`--annotated`  
Filtra por secuencias que han sido anotadas con información de genes/proteínas.  
Línea de comandos: `--annotated true` o `--annotated false`.   
Python: `annotated=True` o `annotated=False`.

`--lab_passaged`  
Filtra a favor o en contra de muestras pasadas en laboratorio.   
Línea de comandos: `--lab_passaged true` para obtener solo muestras pasadas en laboratorio, o `--lab_passaged false` para excluirlas.  
Python: `lab_passaged=True` o `lab_passaged=False`.

`--proteins_complete`  
Bandera para incluir solo secuencias donde todas las proteínas anotadas estén completas.  

`-kt` `--keep_temp`  
Bandera para conservar todos los archivos intermedios/temporales generados durante el procesamiento. Por defecto, solo se conservan los archivos de salida finales.

`--download_all_accessions`  
⚠️ **ADVERTENCIA**: Descarga TODAS las accesiones de virus desde NCBI (toda la taxonomía de Virus, taxon ID 10239). Este es un conjunto de datos extremadamente grande que puede tardar muchas horas en descargarse y requerir un espacio considerable en disco. Úsalo con precaución y asegúrate de tener almacenamiento y ancho de banda adecuados. Cuando esta bandera está activa, el argumento `virus` se ignora.

`-q` `--quiet`  
Solo línea de comandos. Evita que se muestre información de progreso.

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

El archivo CSV de metadatos se verá así:

| accession  | Organism Name | GenBank/RefSeq | Release date | Length | Nuc Completeness | Geographic Location | Host         | ... |
| ---------- | ------------- | -------------- | ------------ | ------ | ---------------- | ------------------- | ------------ | --- |
| KX198135.1 | Zika virus    | GenBank        | 2016-05-18   | 10807  | complete         | Americas:Haiti      | Homo sapiens | ... |
| . . .      | . . .         | . . .          | . . .        | . . .  | . . .            | . . .               | . . .        | ... |

El archivo de resumen del comando (`command_summary.txt`) contendrá, por ejemplo:

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

**Nota**: Si alguna operación falla durante la ejecución (timeouts de API, fallos de descarga de secuencias, fallos de metadatos de GenBank), el resumen incluirá una sección "FAILED OPERATIONS - RETRY COMMANDS" con comandos y URLs exactas que pueden ejecutarse manualmente para reintentar las operaciones fallidas. Por ejemplo:

```
--------------------------------------------------------------------------------
FAILED OPERATIONS - RETRY COMMANDS
--------------------------------------------------------------------------------
Some operations failed during execution. You can retry them manually:

[Failed Sequence Download Batches]
Total failed batches: 2

Batch 15: 200 sequences
Error: HTTPError: 500 Server Error
Accessions: NC_045512.2, MN908947.3, MT020781.1 ... and 197 more
Retry URL: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=NC_045512.2,MN908947.3,...&rettype=fasta&retmode=text

[Failed GenBank Metadata Batches]
Total failed batches: 1
See detailed log file: genbank_failed_batches.log

Accessions: NC_045512.2, MN908947.3, MT020781.1 ... and 2 more
Retry URL: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=NC_045512.2,MN908947.3,...&rettype=gb&retmode=xml

--------------------------------------------------------------------------------
```

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

# Citar

Si utilizas `gget virus` en una publicación, por favor cita los siguientes artículos:

  - Nasri, F. et al (2026). En preparación.

  - Luebbert, L., & Pachter, L. (2023). Efficient querying of genomic reference databases with gget. Bioinformatics. [https://doi.org/10.1093/bioinformatics/btac836](https://doi.org/10.1093/bioinformatics/btac836)

  - O’Leary, N.A., Cox, E., Holmes, J.B. et al (2024). Exploring and retrieving sequence and metadata for species across the tree of life with NCBI Datasets. Sci Data 11, 732. [https://doi.org/10.1038/s41597-024-03571-y](https://doi.org/10.1038/s41597-024-03571-y)

---
---

# Detalles Adicionales: Flujo de trabajo de recuperación de virus

## Visión general

La función `gget.virus()` implementa un flujo de trabajo optimizado de 10 pasos para recuperar secuencias virales y metadatos asociados desde NCBI. El sistema está diseñado para minimizar la sobrecarga de descarga filtrando primero los metadatos y luego descargando solo las secuencias que pasan los filtros iniciales, con recuperación opcional de metadatos detallados de GenBank. Para consultas de SARS-CoV-2 y Alphainfluenza, el flujo de trabajo puede usar paquetes de datos optimizados en caché mientras sigue aplicando todos los filtros y obteniendo metadatos de GenBank.

## Arquitectura

```
┌─────────────────────────────┐
│           Usuarios          │
│                             │
│  • Consulta de virus        │
│    (Taxón/Acc)              │
│  • Criterios de filtrado    │
│    (Hospedador, fechas,     │
│     longitud...)            │
│  • Banderas de salida       │
│    (`--genbank_metadata`)   │
└──────────────┬──────────────┘
               │
               ▼
┌─────────────────────────────┐
│  Verificación de descarga   │
│  en caché                   │
│  (SARS-CoV-2/Alphainfluenza)│
│                             │
│  • Autodetección o banderas │
│  • Descarga de paquetes     │
│    en caché                 │
│  • Aplicar filtros básicos  │
│    (host, complete, lineage)│
│  • Guardar para el pipeline │
└──────────────┬──────────────┘
               │
               ▼
┌─────────────────────────────┐
│   API y prefiltrado         │
│   (o usar metadatos en      │
│    caché)                   │
│                             │
│  • Llama a NCBI Datasets API│
│    O usa metadatos en caché │
│  • Aplica filtros del lado  │
│    del servidor (host,      │
│    refseq)                  │
└──────────────┬──────────────┘
               │
               ▼
┌─────────────────────────────┐
│ Filtrado local de metadatos │
│ y manejo de secuencias      │
│                             │
│  • Aplica TODOS los filtros │
│    locales restantes        │
│    (fechas, recuentos de    │
│     genes, etc.)            │
│  • Genera la lista final de │
│    números de acceso        │
│  • Usa secuencias en caché  │
│    O descarga vía           │
│    E-utilities              │
└──────────────┬──────────────┘
               │
   ┌───────────┴──────────────────────────────────────────┐
   │                                                      │
   ▼                                                      ▼
┌─────────────────────────────┐      ┌───────────────────────────────────┐
│   Procesamiento final       │      │   Metadatos de GenBank (Opcional) │
│                             │      │                                   │
│  • Aplica filtros a nivel   │      │ • Se obtienen incluso para        │
│    de secuencia (p. ej.,    │      │   descargas en caché cuando       │
│    max N's)                 │      │   se solicita                     │
│  • Formatea metadatos       │      │ • Usa la lista final de           │
│    estándar                 │      │   números de acceso               │
└──────────────┬──────────────┘      │ • Se obtienen vía E-utilities API │
               │                     └──────────────────┬────────────────┘
               │                                        │
               └──────────────────┬─────────────────────┘
                                  │
                                  ▼
                    ┌───────────────────────────────┐
                    │ Guardar archivos de salida    │
                    │ finales                       │
                    │                               │
                    │  • _sequences.fasta           │
                    │  • _metadata.csv & .jsonl     │
                    │  • _genbank_metadata.csv      │
                    │    (si se solicita)           │
                    └──────────────┬────────────────┘
                                   │
                                   ▼
                    ┌───────────────────────────────┐
                    │   Resumen y limpieza          │
                    │                               │
                    │  • command_summary.txt        │
                    │  • Mostrar resultados al      │
                    │    usuario                    │
                    │  • Limpiar archivos temporales│
                    └───────────────────────────────┘
```

## Pasos del flujo de trabajo

### Paso 1: Validación de entrada y configuración

* **Función**: función principal `virus()`
* **Propósito**: Validar todos los parámetros del usuario y configurar el logging
* **Operaciones clave**:

  * Validar el formato del taxón/accesión del virus
  * Verificar rangos y formatos de los parámetros de filtrado
  * Configurar la estructura del directorio de salida
  * Configurar el logging según el nivel de verbosidad
  * Verificar oportunidades de optimización para SARS-CoV-2 o Alphainfluenza

### Paso 2: Descarga optimizada en caché (SARS-CoV-2 y Alphainfluenza)

* **Funciones**: `download_sars_cov2_optimized()`, `download_alphainfluenza_optimized()`
* **Propósito**: Usar los paquetes de datos en caché precomputados por NCBI para descargas más rápidas
* **NCBI datasets CLI**: gget incluye el binario de la CLI de NCBI datasets para las principales plataformas (macOS, Linux, Windows). Si ya tienes la CLI `datasets` instalada en tu sistema, gget usará automáticamente la instalación del sistema.
* **Operaciones clave**:

  * Autodetectar o usar banderas explícitas para consultas de SARS-CoV-2/Alphainfluenza
  * Descargar paquetes comprimidos en caché mediante NCBI datasets CLI
  * Aplicar filtros básicos soportados por descargas en caché (host, complete_only, annotated, lineage)
  * Extraer secuencias y metadatos básicos
  * **Guardar datos para continuar el pipeline** (no retorna temprano)
  * Retroceso jerárquico a la API estándar si falla la descarga en caché
* **Filtros aplicados**:

  * ✅ `host` - Aplicado durante la descarga
  * ✅ `complete_only` - Aplicado durante la descarga
  * ✅ `annotated` - Aplicado durante la descarga
  * ✅ `lineage` (solo COVID) - Aplicado durante la descarga
  * ⏭️ Todos los demás filtros se aplican en pasos posteriores

### Paso 3: Recuperación de metadatos

* **Función**: `fetch_virus_metadata()`
* **Propósito**: Recuperar metadatos desde NCBI Datasets API con filtrado del lado del servidor, o usar metadatos en caché
* **Operaciones clave**:

  * **Si se usa descarga en caché**: Omitir llamada a la API, usar metadatos en caché
  * **De lo contrario**: Llamar a NCBI Datasets API con filtros del lado del servidor
  * Aplicar filtros del lado del servidor (host, ubicación geográfica, fecha de liberación, completitud)
  * Manejar paginación de la API con connection pooling
  * Implementar exponential backoff con jitter para reintentos
  * Parsear respuestas JSON con streaming para conjuntos de datos grandes
  * Almacenar metadatos en un formato estructurado con validación

### Paso 4: Filtrado solo de metadatos

* **Función**: `filter_metadata_only()`
* **Propósito**: Aplicar TODOS los filtros locales que no requieren datos de secuencia
* **Operaciones clave**:

  * Filtrar por rangos de fechas con análisis inteligente de fechas
  * Filtrar por completitud del genoma e indicadores de calidad
  * Aplicar filtros de rango numérico (recuentos de genes/proteínas, longitud de secuencia)
  * Manejar metadatos faltantes o malformados de forma robusta
  * Generar lista optimizada de accesiones para procesamiento dirigido
  * **Nota**: Los filtros no aplicados durante la descarga en caché se aplican aquí

### Paso 5: Manejo de secuencias

* **Función**: `download_sequences_by_accessions()`
* **Propósito**: Usar secuencias en caché o descargar secuencias FASTA para accesiones filtradas
* **Operaciones clave**:

  * **Si se usa descarga en caché**: Filtrar secuencias en caché por la lista de accesiones del Paso 4
  * **De lo contrario**: Descargar mediante la API de E-utilities con optimización por lotes
  * Implementar tamaños de lote configurables (por defecto: 200)
  * Hacer streaming de respuestas grandes para gestionar memoria
  * Manejar reintentos de descarga con exponential backoff
  * Devolver la ruta al archivo FASTA para procesamiento

### Paso 6: Filtrado dependiente de la secuencia

* **Función**: `filter_sequences()`
* **Propósito**: Aplicar filtros finales que requieren análisis de secuencia
* **Operaciones clave**:

  * Parsear secuencias FASTA y calcular métricas de secuencia
  * Filtrar por recuento de caracteres ambiguos (`max_ambiguous_chars`)
  * Filtrar por presencia de proteína/gen (`has_proteins`)
  * Filtrar por indicadores de completitud de proteínas (`proteins_complete`)
  * Devolver secuencias filtradas y metadatos actualizados

### Paso 7: Guardar los archivos de salida finales

* **Funciones**: `save_metadata_to_csv()`, `FastaIO.write()`
* **Propósito**: Guardar secuencias filtradas y metadatos en archivos de salida
* **Operaciones clave**:

  * Escribir secuencias filtradas en un archivo FASTA
  * Guardar metadatos en formatos CSV y JSONL
  * Registrar tamaños de archivos de salida para el resumen
  * Validar que los archivos se hayan creado correctamente

### Paso 8: Recuperación de metadatos de GenBank (Opcional)

* **Función**: `fetch_genbank_metadata()`
* **Propósito**: Obtener registros detallados de GenBank para el conjunto final de secuencias
* **Operaciones clave**:

  * **Disponible tanto para descargas en caché como sin caché**
  * Recuperar registros completos de GenBank
  * Extraer 23+ campos de metadatos por registro
  * Procesar en tamaños de lote configurables
  * Implementar rate limiting y reintentos
  * Parsear y validar XML de GenBank
  * Combinar con metadatos existentes

### Paso 9: Resumen final y generación del resumen del comando

* **Función**: `save_command_summary()`
* **Propósito**: Crear un resumen detallado de la ejecución y mostrar resultados
* **Operaciones clave**:

  * Registrar la línea de comandos y parámetros
  * Seguir estadísticas de filtrado en cada etapa
  * Listar archivos de salida con tamaños
  * Documentar operaciones fallidas con comandos de reintento
  * Mostrar un resumen de resultados completo al usuario

### Paso 10: Limpieza

* **Propósito**: Limpiar archivos temporales y finalizar la ejecución
* **Operaciones clave**:

  * Eliminar el directorio temporal de procesamiento (a menos que `keep_temp=True`)
  * Eliminar archivos de metadatos intermedios
  * Conservar el CSV de metadatos de GenBank cuando se recupera con éxito
  * Registrar el estado de finalización

## Dependencias de funciones

```
virus()
├── check_min_max()                          [Paso 1: Validación de entrada]
│   └── Valida pares de parámetros min/max
├── is_sars_cov2_query()                     [Paso 2: Detección de SARS-CoV-2]
│   └── Autodetecta consultas de SARS-CoV-2
├── download_sars_cov2_optimized()           [Paso 2: Descarga en caché]
│   ├── _get_datasets_path()
│   ├── Llamadas a NCBI datasets CLI
│   └── Descarga de paquetes en caché
├── is_alphainfluenza_query()                [Paso 2b: Detección de Alphainfluenza]
│   └── Autodetecta consultas de Alphainfluenza
├── download_alphainfluenza_optimized()      [Paso 2b: Descarga en caché]
│   ├── _get_datasets_path()
│   ├── Llamadas a NCBI datasets CLI
│   └── Descarga de paquetes en caché
├── unzip_file()                             [Paso 2/2b: Extraer datos en caché]
│   └── Utilidades de extracción ZIP
├── fetch_virus_metadata()                   [Paso 3: Recuperación de metadatos API]
│   ├── Cliente de NCBI Datasets API
│   ├── Manejo de paginación
│   ├── Lógica de reintento con backoff
│   └── _get_modified_virus_name() para reintento
├── fetch_virus_metadata_chunked()           [Paso 3: Fallback para conjuntos grandes]
│   └── Estrategia de descarga por bloques de fecha
├── load_metadata_from_api_reports()         [Paso 3: Conversión de metadatos]
│   └── Convierte el formato de la API al formato interno
├── filter_metadata_only()                   [Paso 4: Filtrado de metadatos]
│   ├── parse_date() para comparaciones de fecha
│   ├── Validación numérica
│   └── Manejo de datos faltantes
├── download_sequences_by_accessions()       [Paso 5: Descarga de secuencias]
│   ├── Cliente de la API de E-utilities
│   ├── Procesamiento por lotes (por defecto: 200)
│   └── Manejo de streaming
├── filter_sequences()                       [Paso 6: Filtrado de secuencias]
│   ├── Parser FastaIO
│   └── Validación de secuencias
├── save_metadata_to_csv()                   [Paso 7: Guardar salidas]
│   └── Formateo y escritura CSV
├── fetch_genbank_metadata()                 [Paso 8: Datos opcionales de GenBank]
│   ├── _fetch_genbank_batch()
│   ├── _clean_xml_declarations()
│   ├── Utilidades de parseo XML
│   └── Rate limiting
├── save_genbank_metadata_to_csv()           [Paso 8: Guardar datos de GenBank]
│   └── Combina con metadatos del virus
└── save_command_summary()                   [Paso 9: Resumen de ejecución]
    └── Seguimiento de operaciones fallidas
```

## Características de optimización

### 1. **Filtrado del lado del servidor**

* Aplica filtros a nivel de la API de NCBI para reducir la transferencia de datos
* Filtros soportados: host, ubicación geográfica, fecha de liberación, completitud del genoma
* Validación automática de compatibilidad y valores de filtros

### 2. **Filtrado en múltiples etapas**

* **Etapa 1**: Filtros solo de metadatos (rápido, sin descarga de secuencias)
* **Etapa 2**: Filtros dependientes de la secuencia (conjunto prefiltrado)
* **Etapa 3**: Integración y filtrado de metadatos de GenBank
* **Etapa 4**: Validación final y controles de calidad

### 3. **Descargas optimizadas**

* Tamaños de lote configurables para diferentes tipos de datos
* Connection pooling para mejorar el rendimiento
* Manejo de streaming para descargas grandes
* Mecanismos de rate limiting y reintento

### 4. **Descargas optimizadas en caché**

* Manejo especial para consultas de SARS-CoV-2 y Alphainfluenza usando paquetes de datos en caché de NCBI
* Detección automática o banderas explícitas (`--is_sars_cov2`, `--is_alphainfluenza`)
* Estrategias de fallback jerárquicas a la API estándar si falla la descarga en caché
* Descargas significativamente más rápidas para conjuntos de datos grandes
* **Continuación del pipeline**: Las descargas en caché ahora continúan por todos los pasos del flujo de trabajo
* **Filtrado posterior a la descarga**: Los filtros no aplicados durante la descarga en caché se aplican después
* **Metadatos de GenBank**: Disponible para descargas en caché cuando se usa la bandera `--genbank_metadata`
* **Categorías de filtros**:

  * Aplicados durante la descarga: `host`, `complete_only`, `annotated`, `lineage` (COVID)
  * Aplicados después: todos los demás filtros (longitud de secuencia, recuentos de genes, fechas, etc.)

### 5. **Estructuras de datos eficientes**

* Diccionarios basados en accesión para búsquedas O(1)
* Parsers en streaming para JSON y XML
* Manejo de FASTA eficiente en memoria
* Combinación de metadatos optimizada

## Archivos de salida

### 1. **Secuencias FASTA** (`{virus}_sequences.fasta`)

* Contiene secuencias de nucleótidos para resultados filtrados
* Formato FASTA estándar con encabezados detallados
* Se conserva la orientación original de NCBI
* Anotaciones opcionales de proteínas/segmentos en los encabezados

### 2. **Metadatos CSV** (`{virus}_metadata.csv`)

* Formato tabular para análisis en hojas de cálculo
* Estructura de columnas estandarizada
* Información geográfica y taxonómica
* Detalles de recolección y envío
* Métricas de calidad y anotaciones

### 3. **Metadatos de GenBank** (`{virus}_genbank_metadata.csv`) [Opcional]

* 23+ columnas detalladas de metadatos
* Referencias de publicaciones
* Anotaciones de características
* Referencias cruzadas a otras bases de datos
* Detalles de cepa e aislamiento

### 4. **Metadatos JSONL** (`{virus}_metadata.jsonl`)

* Formato JSON Lines para metadatos de virus después del filtrado solo de metadatos
* Formato amigable para streaming y acceso programático
* Un objeto JSON por secuencia con los mismos campos que el CSV de metadatos
* Los campos específicos de GenBank se almacenan por separado en `{virus}_genbank_metadata.csv` cuando se utiliza `--genbank_metadata`

### 5. **Resumen del comando** (`command_summary.txt`)

* Resumen generado automáticamente de la ejecución del comando
* Registra la línea de comandos exacta que se ejecutó
* Estado de ejecución (éxito/fallo con mensajes de error)
* Estadísticas de filtrado en cada etapa:

  * Registros totales desde la API
  * Registros después del filtrado de metadatos
  * Secuencias finales después de todos los filtros
* Estadísticas detalladas:

  * Hospedadores únicos con recuentos (hasta 20 principales)
  * Ubicaciones geográficas únicas con recuentos (hasta 20 principales)
  * Rango de longitudes de secuencia y promedio
  * Desglose de completitud (complete vs partial)
  * Desglose de base de datos de origen (GenBank vs RefSeq)
  * Países de remitentes únicos con recuentos (hasta 20 principales)
* Lista de todos los archivos de salida generados con tamaños
* **Seguimiento de operaciones fallidas** (cuando aplica):

  * **Fallos por timeout de API**: URL exacta que expiró con sugerencias alternativas
  * **Lotes de descarga de secuencias fallidos**: números de lote, listas de accesiones y URLs de reintento
  * **Lotes de metadatos de GenBank fallidos**: listas de accesiones con URLs de reintento individuales
  * Todas las operaciones fallidas incluyen comandos/URLs exactos que pueden ejecutarse manualmente para reintentar

## Características de rendimiento

### Escalabilidad

* **Conjuntos pequeños** (< 1,000 secuencias): procesamiento casi instantáneo
* **Conjuntos medianos** (1,000 - 10,000 secuencias): minutos para completarse
* **Conjuntos grandes** (> 10,000 secuencias): paginación y filtrado optimizados

### Uso de memoria

* El procesamiento en streaming minimiza el uso de memoria
* Metadatos cacheados en memoria para operaciones de filtrado
* Archivos FASTA grandes procesados en bloques

### Eficiencia de red

* Mínimas llamadas a la API gracias al filtrado del lado del servidor
* Descargas dirigidas reducen el uso de ancho de banda
* Reintento automático con exponential backoff

## Manejo de errores

### Fallos de API

* Estrategia de reintento inteligente con exponential backoff y jitter
* Detección de errores del lado del servidor con guía específica:

  * Manejo de timeouts para conjuntos grandes
  * Sugerencias para optimizar filtros geográficos
  * Ajustes de tamaño de lote para metadatos de GenBank
* Connection pooling y gestión de sesiones
* Registro detallado de errores con pasos de troubleshooting

### Validación de datos

* Validación integral de parámetros de entrada:

  * Verificación de tipos para todos los parámetros
  * Validación de rangos para valores numéricos
  * Validación de formato y rango de fechas
  * Normalización de parámetros booleanos
* Verificación de integridad de secuencias:

  * Validación de formato FASTA
  * Detección de caracteres ambiguos
  * Chequeos de completitud de proteínas/genes
* Validación de consistencia de metadatos:

  * Verificación de presencia de campos requeridos
  * Validación de tipo de datos
  * Validación de referencias cruzadas
  * Validación de registros de GenBank

### Mecanismos de recuperación

* Limpieza automática de archivos temporales
* Preservación parcial de resultados:

  * Guardado de metadatos intermedios
  * Guardado progresivo del estado de filtrado
  * Caché de metadatos de GenBank
* Estrategias de fallback jerárquicas:

  * Paquetes optimizados de SARS-CoV-2
  * Fallback a datos en caché
  * Fallback a recuperación basada en API
* Reporte detallado de errores:

  * Análisis de causa raíz
  * Sugerencias de comandos alternativos
  * Recomendaciones para relajar filtros
  * Consejos de optimización de rendimiento

## Ejemplos de uso

### Ejemplos de línea de comandos

```bash
# Get help and see all available parameters
$ gget virus --help

$ gget virus "Nipah virus"

# Download Zika virus sequences with basic filtering (API + metadata filtering)
$ gget virus "Zika virus" --host human --min_seq_length 10000 --max_seq_length 11000

# Download with metadata and sequence filtering
$ gget virus "Ebolavirus" --max_seq_length 20000 --genbank_metadata -o ./ebola_data

# Download SARS-CoV-2 with cached optimization
$ gget virus "SARS-CoV-2" --host dog --nuc_completeness complete

# Download Influenza A with post-download sequence filtering (warning: big data size)
$ gget virus "Influenza A virus" --host human --max_ambiguous_chars 50 --has_proteins spike

# Using accession ID to get data
$ gget virus -a "MK947457" --host deer --min_collection_date "2020-01-01"
```

### Ejemplos en Python

```python
  import gget
  import pandas as pd
  from Bio import SeqIO

  # Basic download with GenBank metadata
  gget.virus(
    "Zika virus",
    host="human",
    genbank_metadata=True,
    outfolder="zika_data"
  )

  # Access different data types from output files
  sequences = list(SeqIO.parse("zika_data/Zika_virus_sequences.fasta", "fasta"))
  virus_metadata = pd.read_csv("zika_data/Zika_virus_metadata.csv")
  genbank_metadata = pd.read_csv("zika_data/Zika_virus_genbank_metadata.csv")

  # Print GenBank metadata summary
  for _, row in genbank_metadata.head().iterrows():
    print(f"Sequence: {row['accession']}")
    print(f"  Length: {row['sequence_length']} bp")
    print(f"  Host: {row.get('host', 'Unknown')}")
    print(f"  Location: {row.get('geographic_location', 'Unknown')}")
    print(f"  Collection date: {row.get('collection_date', 'Unknown')}")

  # Advanced filtering with GenBank data
  gget.virus(
    "SARS-CoV-2", 
    host="human",
    min_seq_length=29000,
    max_seq_length=30000,
    min_collection_date="2020-03-01",
    max_collection_date="2020-03-31",
    geographic_location="North America",
    genbank_metadata=True,
    genbank_batch_size=200,
    outfolder="covid_march2020"
  )

  # Process and analyze results

  # Read virus metadata
  virus_df = pd.read_csv("covid_march2020/SARS-CoV-2_metadata.csv")
  print(f"Total sequences: {len(virus_df)}")
  print(f"Unique hosts: {virus_df['Host'].nunique()}")
  print(f"Date range: {virus_df['Collection Date'].min()} to {virus_df['Collection Date'].max()}")

  # Read GenBank metadata for detailed analysis
  genbank_df = pd.read_csv("covid_march2020/SARS-CoV-2_genbank_metadata.csv")
  print(f"Sequences with GenBank data: {len(genbank_df)}")
  print("\nPublication summary:")
  print(genbank_df['reference_count'].describe())

  # Custom sequence analysis
  sequences = list(SeqIO.parse("covid_march2020/SARS-CoV-2_sequences.fasta", "fasta"))
  for record in sequences:
    gc_content = (str(record.seq).count('G') + str(record.seq).count('C')) / len(record.seq)
    print(f"{record.id}: GC content = {gc_content:.2%}")

  # Merge metadata sources
  merged_df = pd.merge(
    virus_df,
    genbank_df,
    on='accession',
    how='left',
    suffixes=('_virus', '_genbank')
  )

  # Save merged analysis
  merged_df.to_csv("covid_march2020/combined_analysis.csv", index=False)
```

### Ejemplos de estrategia de análisis

Los ejemplos anteriores demuestran diferentes enfoques de análisis:

1. **Integración básica de GenBank**: Obtener secuencias con metadatos de GenBank para un análisis completo
2. **Filtrado avanzado**: Combinar metadatos de virus y datos de GenBank con filtros personalizados
3. **Análisis personalizado**: Procesar secuencias y metadatos usando BioPython y Pandas
4. **Integración de datos**: Unir metadatos de virus y GenBank para un análisis detallado

### Acceso programático

```python
# Access filtered metadata and sequences
metadata_file = "covid_data/SARS-CoV-2_metadata.jsonl"
sequences_file = "covid_data/SARS-CoV-2_sequences.fasta"

# Process results with custom analysis
import json
with open(metadata_file) as f:
    for line in f:
        record = json.loads(line)
        # Custom analysis here
```
