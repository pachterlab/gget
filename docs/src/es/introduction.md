[![](https://img.shields.io/static/v1?label=Sponsor&message=%E2%9D%A4&logo=GitHub&color=%23fe8e86)](https://github.com/sponsors/lauraluebbert)
[![pypi version](https://img.shields.io/pypi/v/gget)](https://pypi.org/project/gget)
[![Downloads](https://static.pepy.tech/personalized-badge/gget?period=total&units=international_system&left_color=grey&right_color=brightgreen&left_text=Downloads)](https://pepy.tech/project/gget)
[![Conda](https://img.shields.io/conda/dn/bioconda/gget?logo=Anaconda)](https://anaconda.org/bioconda/gget)
[![license](https://img.shields.io/pypi/l/gget)](LICENSE)
![status](https://github.com/pachterlab/gget/actions/workflows/ci.yml/badge.svg)
[![Star on GitHub](https://img.shields.io/github/stars/pachterlab/gget.svg?style=social)](https://github.com/pachterlab/gget/)  

[<img align="right" width="50%" height="50%" src="https://github.com/pachterlab/gget/blob/main/docs/assets/website_v2_gget_overview.png?raw=true" />](https://raw.githubusercontent.com/pachterlab/gget/main/figures/gget_overview.png)

# ¡Bienvenidos!
  
`gget` es un programa gratuito de código fuente abierta de Terminal y Python que permite la consulta eficiente de bases de datos genómicas.  
<br>
`gget` consiste en un conjunto de módulos separados pero interoperables, cada uno diseñado para facilitar un tipo de consulta de base de datos en una sola línea de código.  
<br>
Las bases de datos consultadas por `gget` se actualizan continuamente, lo que a veces cambia su estructura. Los módulos `gget` se prueban automáticamente cada dos semanas y se actualizan para que coincidan con las nuevas estructuras de la base de datos cuando es necesario. Si encuentra algún problema, actualice a la última versión de `gget` usando `pip install --upgrade gget`. Si el problema persiste, [informa el problema](https://github.com/pachterlab/gget/issues/new/choose).  
<br>
[<kbd> <br> Solicitar una nueva función <br> </kbd>](https://github.com/pachterlab/gget/issues/new/choose)
<br>

# Módulos gget

Estos son los módulos principales de `gget`. Haga clic en cualquier módulo para acceder a la documentación detallada.

<table style="width:100%; table-layout:fixed;">
  <tr>
    <td style="width:33.33%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/es/alphafold.md"><strong>gget alphafold</strong></a><br>Predecir la estructura 3D de una proteína a partir de una secuencia de aminoácidos.</td>
    <td style="width:33.33%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/es/archs4.md"><strong>gget archs4</strong></a><br>¿Cuál es la expresión de mi gen en el tejido X?</td>
    <td style="width:33.33%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/es/bgee.md"><strong>gget bgee</strong></a><br>Encontrar todos los ortólogos de un gen.</td>
  </tr>
  <tr>
    <td style="width:33.33%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/es/blast.md"><strong>gget blast</strong></a><br>Realizar un BLAST de una secuencia de nucleótidos o aminoácidos.</td>
    <td style="width:33.33%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/es/blat.md"><strong>gget blat</strong></a><br>Encontrar la ubicación genómica de una secuencia de nucleótidos o aminoácidos.</td>
    <td style="width:33.33%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/es/cbio.md"><strong>gget cbio</strong></a><br>Explorar la expresión de un gen en los cánceres especificados.</td>
  </tr>
  <tr>
    <td style="width:33.33%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/es/cellxgene.md"><strong>gget cellxgene</strong></a><br>Obtener matrices de conteo de ARN de células individuales listas para usar para ciertos tejidos/enfermedades/etc.</td>
    <td style="width:33.33%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/es/cosmic.md"><strong>gget cosmic</strong></a><br>Buscar genes, mutaciones y otros factores asociados con ciertos cánceres.</td>
    <td style="width:33.33%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/es/diamond.md"><strong>gget diamond</strong></a><br>Alinear secuencias de aminoácidos a una referencia.</td>
  </tr>
  <tr>
    <td style="width:33.33%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/es/elm.md"><strong>gget elm</strong></a><br>Encontrar dominios y funciones de interacción de proteínas en una secuencia de aminoácidos.</td>
    <td style="width:33.33%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/es/enrichr.md"><strong>gget enrichr</strong></a><br>Verificar si una lista de genes está asociada con un tipo celular específico/ vía/ enfermedad/ etc.</td>
    <td style="width:33.33%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/es/info.md"><strong>gget info</strong></a><br>Recuperar toda la información asociada con un ID de Ensembl.</td>
  </tr>
  <tr>
    <td style="width:33.33%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/es/muscle.md"><strong>gget muscle</strong></a><br>Alinear múltiples secuencias de nucleótidos o aminoácidos entre sí.</td>
    <td style="width:33.33%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/es/mutate.md"><strong>gget mutate</strong></a><br>Mutar secuencias de nucleótidos según mutaciones específicas.</td>
    <td style="width:33.33%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/es/opentargets.md"><strong>gget opentargets</strong></a><br>Explorar qué enfermedades y medicamentos están asociados con un gen.</td>
  </tr>
  <tr>
    <td style="width:33.33%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/es/pdb.md"><strong>gget pdb</strong></a><br>Recuperar datos de la Base de Datos de Proteínas (PDB) según un ID de PDB.</td>
    <td style="width:33.33%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/es/ref.md"><strong>gget ref</strong></a><br>Obtener genomas de referencia de Ensembl.</td>
    <td style="width:33.33%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/es/search.md"><strong>gget search</strong></a><br>Encontrar IDs de Ensembl asociados con la palabra de búsqueda especificada.</td>
  </tr>
  <tr>
    <td style="width:33.33%; padding:20px; text-align:center; vertical-align:top;"><a href="/gget/es/seq.md"><strong>gget seq</strong></a><br>Recuperar la secuencia de nucleótidos o aminoácidos de un gen.</td>
  </tr>
</table>

<br>  

Si usa `gget` en una publicación, por favor [cite*](cite.md):    
```
Luebbert, L., & Pachter, L. (2023). Efficient querying of genomic reference databases with gget. Bioinformatics. https://doi.org/10.1093/bioinformatics/btac836
```
Lea el artículo aquí: [https://doi.org/10.1093/bioinformatics/btac836](https://doi.org/10.1093/bioinformatics/btac836)

Gracias a [Victor Garcia-Ruiz](https://github.com/victorg775) y [Anna Karen Orta](https://www.linkedin.com/in/akorta/) por su ayuda con la traduccion del sitio web.  

<br>  
<iframe width="560" height="315" src="https://www.youtube.com/embed/cVR0k6Mt97o?si=BJwRyaymmxF9w65f" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen></iframe>

![logo-bmbf](https://github.com/user-attachments/assets/68449ba0-09d1-4f1e-b747-5ae2fec08a21)
![logo-okfn](https://github.com/user-attachments/assets/452ae8d8-69f0-4d0d-848c-ddfb40357eb2)
