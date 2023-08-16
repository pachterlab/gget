## ✨ ¡Lo más reciente!
**Versión ≥ 0.27.9** (7 de agosto de 2023):
- Nuevos argumentos para [`gget enrichr`](es/enrichr.md): use el argumento `background_list` para proporcionar una lista de genes 'background'
- [`gget search`](es/search.md) ahora también busca sinónimos [Ensembl](https://ensembl.org/) (además de nombres y descripciones de genes) para obtener resultados de búsqueda más completos (gracias a [Samuel Klein](https://github.com/KleinSamuel) por la [sugerencia](https://github.com/pachterlab/gget/issues/90))
  
**Versión ≥ 0.27.8** (12 de julio de 2023):
- Nuevo argumento para [`gget search`](es/search.md): especifique la versión de Ensembl desde la cual se obtiene la información con `-r` `--release`
- Se corrigió un [error](https://github.com/pachterlab/gget/issues/91) en [`gget pdb`](es/pdb.md) (este error se introdujo en la versión 0.27.5)

**Versión ≥ 0.27.7** (15 de mayo de 2023):
- Se movieron las dependencias para los módulos [`gget gpt`](es/gpt.md) y [`gget cellxgene`](es/cellxgene.md) de los requisitos instalados automáticamente a [`gget setup`](es/setup.md)
- Dependencias [`gget alphafold`](es/alphafold.md) actualizadas para compatibilidad con Python >= 3.10
- Se agregó el argumento `census_version` a [`gget cellxgene`](es/cellxgene.md)

**Versión ≥ 0.27.6** (1 de mayo de 2023) (TIRO debido a problemas con las dependencias -> reemplazada por la versión 0.27.7):  
- Gracias a el PR de [Tomás Di Domenico](https://github.com/tdido): [`gget search`](es/search.md) ahora también puede consultar los ID de plantas 🌱 Ensembl  
- Nuevo módulo: [`gget cellxgene`](es/cellxgene.md)

**Versión ≥ 0.27.5** (6 de abril de 2023):
- Se actualizó [`gget search`](es/search.md) para que funcione correctamente con la nueva versión de [Pandas](https://pypi.org/project/pandas/2.0.0/) 2.0.0 (lanzado el 3 de abril de 2023), además de versiones anteriores de Pandas
- Se actualizó [`gget info`](es/info.md) con nuevos banderas `uniprot` y `ncbi` que permiten desactivar los resultados de estas bases de datos de forma independiente para ahorrar tiempo de ejecución (nota: el indicador `ensembl_only` quedó obsoleto)
- Todos los módulos gget ahora tienen una bandera `-q / --quiet` (para Python: `verbose=False`) para desactivar la información de progreso

**Versión ≥ 0.27.4** (19 de marzo de 2023):
- Nuevo módulo: [`gget gpt`](es/gpt.md) 

**Versión ≥ 0.27.3** (11 de marzo de 2023):
- [`gget info`](es/info.md) excluye los ID de PDB de forma predeterminada para aumentar la velocidad (los resultados de PDB se pueden incluir usando la marca `--pdb` / `pdb=True`).

**Versión ≥ 0.27.2** (1 de enero de 2023):
- Se actualizó [`gget alphafold`](es/alphafold.md) a [DeepMind's AlphaFold v2.3.0](https://github.com/deepmind/alphafold/releases/tag/v2.3.0) (incluidos los nuevos argumentos `multimer_for_monomer ` y `multimer_recycles`)

**Versión ≥ 0.27.0** (10 de diciembre de 2022):
- Se actualizó [`gget alphafold`](es/alphafold.md) para que coincida con los cambios recientes de DeepMind
- Número de versión actualizado para que coincida con la edad de [el creador de gget](https://github.com/lauraluebbert) siguiendo una larga tradición de laboratorio de Pachter

**Versión ≥ 0.3.13** (11 de noviembre de 2022):
- Tiempo de ejecución reducido para [`gget enrichr`](es/enrichr.md) y [`gget archs4`](es/archs4.md) cuando se usa con ID de Ensembl

**Versión ≥ 0.3.12** (10 de noviembre de 2022):
- [`gget info`](es/info.md) ahora también devuelve datos de localización subcelular de UniProt
- El nuevo indicador [`gget info`](es/info.md) `ensembl_only` devuelve solo los resultados de Ensembl
- Tiempo de ejecución reducido para [`gget info`](es/info.md) y [`gget seq`](es/seq.md)

**Versión ≥ 0.3.11** (7 de septiembre de 2022):
- Nuevo módulo: [`gget pdb`](es/pdb.md)

**Versión ≥ 0.3.10** (2 de septiembre de 2022):
- [`gget alphafold`](es/alphafold.md) ahora también devuelve valores pLDDT para generar gráficos sin volver a ejecutar el programa (consulte también las [preguntas frecuentes de gget alphafold](https://github.com/pachterlab/gget/discusiones/39))

**Versión ≥ 0.3.9** (25 de agosto de 2022):
- Instrucciones de instalación de openmm actualizadas para [`gget alphafold`](es/alphafold.md)

**Versión ≥ 0.3.8** (12 de agosto de 2022):
- Se corrigieron los requisitos de versión de mysql-connector-python

**Versión ≥ 0.3.7** (9 de agosto de 2022):
- **NOTA:** El [sitio FTP de Ensembl](http://ftp.ensembl.org/pub/) cambió su estructura el 8 de agosto de 2022. Actualice a la versión `gget` ≥ 0.3.7 si usa [`obtener ref`](es/ref.md)

**Versión ≥ 0.3.5** (6 de agosto de 2022):
- Nuevo módulo: [`gget alphafold`](es/alphafold.md)

**Versión ≥ 0.2.6** (7 de julio de 2022):
- ¡[`gget ref`](es/ref.md) ahora admite genomas de plantas! 🌱

**Versión ≥ 0.2.5** (30 de junio de 2022):
- **NOTA:** [UniProt](https://www.uniprot.org/) cambió la estructura de su API el 28 de junio de 2022. Actualice a la versión `gget` ≥ 0.2.5 si usa alguno de los módulos que consultan datos de UniProt ([`gget info`](es/info.md) y [`gget seq`](es/seq.md)).

**Versión ≥ 0.2.3:** (26 de junio de 2022):
- JSON ahora es el formato de regreso predeterminado para la Terminal para los módulos que anteriormente devolvían el formato de data frame (CSV) (el formato se puede convertir a data frame/CSV usando la bandera `[-csv][--csv]`). El formato data frame/CSV sigue siendo el formato de regreso predeterminada para Python (Jupyter Lab/Google Colab) (y se puede convertir a JSON con `json=True`).
- Para todos los módulos, el primer parámetro requerido se convirtió en un parámetro posicional y ya no debe nombrarse en la línea de comandos, p. ej. `gget ref -s human` &rarr; `gget ref human`.
- `gget info`: `[--expand]` está en desuso. El módulo ahora siempre devolverá toda la información disponible.
- Ligeros cambios en la salida devuelta por `gget info`, incluida la devolución de los ID de Ensembl versionados.
- `gget info` y `gget seq` ahora son compatibles con las IDs de WormBase y FlyBase.
- Ahora también se pueden ingresar IDs de tipo Ensembl a `gget archs4` y `gget enrichr` con la bandera `[-e][--ensembl]` (`ensembl=True` para Python (Jupyter Lab / Google Colab)).
- El parámetro `seqtype` de `gget seq` fue reemplazado por la bandera `[-t][--translate]` (`translate=True/False` para Python (Jupyter Lab / Google Colab)) que devolverá secuencias de nucleótidos (`False`) o aminoácidos (`True`).
- El parámetro `seqtype` de `gget search` se renombró a `id_type` (aún tomando los mismos parámetros 'gene' o 'transcript').
