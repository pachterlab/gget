## ✨ ¡Lo más reciente!  
**Versión ≥ 0.28.0** (5 de noviembre de 2023):
- Documentación actualizada de [`gget muscle`](./muscle.md) para agregar un tutorial sobre cómo visualizar secuencias con diferentes longitudes de nombres de secuencia + ligero cambio en la visualización devuelta para que sea un poco más sólida ante diferentes nombres de secuencia  
- [`gget muscle`](./muscle.md) ahora también permite una lista de secuencias como entrada (como alternativa a proporcionar la ruta a un archivo FASTA)
- Permitir filtro de genes faltante para [`gget cellxgene`](cellxgene.md) (corrige [error](https://github.com/pachterlab/gget/issues/110))
- [`gget seq`](./seq.md): permite nombres de genes faltantes (correccione [https://github.com/pachterlab/gget/issues/107](https://github.com/pachterlab/gget /números/107))  
- Nuevos argumentos para [`gget enrichr`](enrichr.md): use el argumento `kegg_out` y `kegg_rank` para crear una imagen de la vía KEGG con los genes del análisis de enriquecimiento resaltados (gracias a [este PR](https ://github.com/pachterlab/gget/pull/106) por [Noriaki Sato](https://github.com/noriakis))  
- Nuevos módulos: [`gget elm`](elm.md) y [`gget Diamond`](diamond.md)
  
**Versión ≥ 0.27.9** (7 de agosto de 2023):
- Nuevos argumentos para [`gget enrichr`](enrichr.md): use el argumento `background_list` para proporcionar una lista de genes 'background'
- [`gget search`](search.md) ahora también busca sinónimos [Ensembl](https://ensembl.org/) (además de nombres y descripciones de genes) para obtener resultados de búsqueda más completos (gracias a [Samuel Klein](https://github.com/KleinSamuel) por la [sugerencia](https://github.com/pachterlab/gget/issu90))
  
**Versión ≥ 0.27.8** (12 de julio de 2023):
- Nuevo argumento para [`gget search`](search.md): especifique la versión de Ensembl desde la cual se obtiene la información con `-r` `--release`
- Se corrigió un [error](https://github.com/pachterlab/gget/issu91) en [`gget pdb`](pdb.md) (este error se introdujo en la versión 0.27.5)

**Versión ≥ 0.27.7** (15 de mayo de 2023):
- Se movieron las dependencias para los módulos [`gget gpt`](gpt.md) y [`gget cellxgene`](cellxgene.md) de los requisitos instalados automáticamente a [`gget setup`](setup.md)
- Dependencias [`gget alphafold`](alphafold.md) actualizadas para compatibilidad con Python >= 3.10
- Se agregó el argumento `census_version` a [`gget cellxgene`](cellxgene.md)

**Versión ≥ 0.27.6** (1 de mayo de 2023) (TIRO debido a problemas con las dependencias -> reemplazada por la versión 0.27.7):  
- Gracias a el PR de [Tomás Di Domenico](https://github.com/tdido): [`gget search`](search.md) ahora también puede consultar los ID de plantas 🌱 Ensembl  
- Nuevo módulo: [`gget cellxgene`](cellxgene.md)

**Versión ≥ 0.27.5** (6 de abril de 2023):
- Se actualizó [`gget search`](search.md) para que funcione correctamente con la nueva versión de [Pandas](https://pypi.org/project/pandas/2.0.0/) 2.0.0 (lanzado el 3 de abril de 2023), además de versiones anteriores de Pandas
- Se actualizó [`gget info`](info.md) con nuevos banderas `uniprot` y `ncbi` que permiten desactivar los resultados de estas bases de datos de forma independiente para ahorrar tiempo de ejecución (nota: el indicador `ensembl_only` quedó obsoleto)
- Todos los módulos gget ahora tienen una bandera `-q / --quiet` (para Python: `verbose=False`) para desactivar la información de progreso

**Versión ≥ 0.27.4** (19 de marzo de 2023):
- Nuevo módulo: [`gget gpt`](gpt.md) 

**Versión ≥ 0.27.3** (11 de marzo de 2023):
- [`gget info`](info.md) excluye los ID de PDB de forma predeterminada para aumentar la velocidad (los resultados de PDB se pueden incluir usando la marca `--pdb` / `pdb=True`).

**Versión ≥ 0.27.2** (1 de enero de 2023):
- Se actualizó [`gget alphafold`](alphafold.md) a [DeepMind's AlphaFold v2.3.0](https://github.com/deepmind/alphafold/releases/tag/v2.3.0) (incluidos los nuevos argumentos `multimer_for_monomer ` y `multimer_recycles`)

**Versión ≥ 0.27.0** (10 de diciembre de 2022):
- Se actualizó [`gget alphafold`](alphafold.md) para que coincida con los cambios recientes de DeepMind
- Número de versión actualizado para que coincida con la edad de [el creador de gget](https://github.com/lauraluebbert) siguiendo una larga tradición de laboratorio de Pachter

**Versión ≥ 0.3.13** (11 de noviembre de 2022):
- Tiempo de ejecución reducido para [`gget enrichr`](enrichr.md) y [`gget archs4`](archs4.md) cuando se usa con ID de Ensembl

**Versión ≥ 0.3.12** (10 de noviembre de 2022):
- [`gget info`](info.md) ahora también devuelve datos de localización subcelular de UniProt
- El nuevo indicador [`gget info`](info.md) `ensembl_only` devuelve solo los resultados de Ensembl
- Tiempo de ejecución reducido para [`gget info`](info.md) y [`gget seq`](seq.md)

**Versión ≥ 0.3.11** (7 de septiembre de 2022):
- Nuevo módulo: [`gget pdb`](pdb.md)

**Versión ≥ 0.3.10** (2 de septiembre de 2022):
- [`gget alphafold`](alphafold.md) ahora también devuelve valores pLDDT para generar gráficos sin volver a ejecutar el programa (consulte también las [preguntas frecuentes de gget alphafold](https://github.com/pachterlab/gget/discusion39))

**Versión ≥ 0.3.9** (25 de agosto de 2022):
- Instrucciones de instalación de openmm actualizadas para [`gget alphafold`](alphafold.md)

**Versión ≥ 0.3.8** (12 de agosto de 2022):
- Se corrigieron los requisitos de versión de mysql-connector-python

**Versión ≥ 0.3.7** (9 de agosto de 2022):
- **NOTA:** El [sitio FTP de Ensembl](http://ftp.ensembl.org/pub/) cambió su estructura el 8 de agosto de 2022. Actualice a la versión `gget` ≥ 0.3.7 si usa [`obtener ref`](ref.md)

**Versión ≥ 0.3.5** (6 de agosto de 2022):
- Nuevo módulo: [`gget alphafold`](alphafold.md)

**Versión ≥ 0.2.6** (7 de julio de 2022):
- ¡[`gget ref`](ref.md) ahora admite genomas de plantas! 🌱

**Versión ≥ 0.2.5** (30 de junio de 2022):
- **NOTA:** [UniProt](https://www.uniprot.org/) cambió la estructura de su API el 28 de junio de 2022. Actualice a la versión `gget` ≥ 0.2.5 si usa alguno de los módulos que consultan datos de UniProt ([`gget info`](info.md) y [`gget seq`](seq.md)).

**Versión ≥ 0.2.3:** (26 de junio de 2022):
- JSON ahora es el formato de regreso predeterminado para la Terminal para los módulos que anteriormente devolvían el formato de data frame (CSV) (el formato se puede convertir a data frame/CSV usando la bandera `[-csv][--csv]`). El formato data frame/CSV sigue siendo el formato de regreso predeterminada para Python (Jupyter Lab/Google Colab) (y se puede convertir a JSON con `json=True`).
- Para todos los módulos, el primer parámetro requerido se convirtió en un parámetro posicional y ya no debe nombrarse en la línea de comandos, p. ej. `gget ref -s human` &rarr; `gget ref human`.
- `gget info`: `[--expand]` está en desuso. El módulo ahora siempre devolverá toda la información disponible.
- Ligeros cambios en la salida devuelta por `gget info`, incluida la devolución de los ID de Ensembl versionados.
- `gget info` y `gget seq` ahora son compatibles con las IDs de WormBase y FlyBase.
- Ahora también se pueden ingresar IDs de tipo Ensembl a `gget archs4` y `gget enrichr` con la bandera `[-e][--ensembl]` (`ensembl=True` para Python (Jupyter Lab / Google Colab)).
- El parámetro `seqtype` de `gget seq` fue reemplazado por la bandera `[-t][--translate]` (`translate=True/False` para Python (Jupyter Lab / Google Colab)) que devolverá secuencias de nucleótidos (`False`) o aminoácidos (`True`).
- El parámetro `seqtype` de `gget search` se renombró a `id_type` (aún tomando los mismos parámetros 'gene' o 'transcript').
