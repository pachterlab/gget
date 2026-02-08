[<kbd> Ver el codigo fuente de la pagina en GitHub </kbd>](https://github.com/pachterlab/gget/blob/main/docs/src/es/opentargets.md)

> Parámetros de Python són iguales a los parámetros largos (`--parámetro`) de Terminal, si no especificado de otra manera. Banderas son parámetros de verdadero o falso (True/False) en Python. El manuál para cualquier modulo de gget se puede llamar desde la Terminal con la bandera `-h` `--help`.  
# gget opentargets 🎯
**Obtenga enfermedades o fármacos asociados con ciertos genes desde [OpenTargets](https://platform.opentargets.org/).**  
Formato de salida: JSON/CSV (línea de comandos) o marco de datos (Python).  

Este módulo fue escrito por [Sam Wagenaar](https://github.com/techno-sam).  

**Argumento posicional**  
`ens_id`  
ID de gen Ensembl, por ejemplo, ENSG00000169194.

**Argumentos opcionales**  
`-r` `--resource`   
Define el tipo de información a devolver en la salida. Predeterminado: 'diseases' (enfermedades).   
Los recursos posibles son:

| Recurso            | Valor devuelto                                                    | Filtros válidos                                   | Fuentes                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
|--------------------|-------------------------------------------------------------------|---------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `diseases`         | Enfermedades asociadas                                            | Ninguno                                           | Varias:<ul><li>[Open&nbsp;Targets](https://genetics.opentargets.org/)</li><li>[ChEMBL](https://www.ebi.ac.uk/chembl/)</li><li>[Europe&nbsp;PMC](http://europepmc.org/)</li></ul>etc.                                                                                                                                                                                                                                                                                                               |
| `drugs`            | Fármacos asociados                                                | `disease_id`                                      | [ChEMBL](https://www.ebi.ac.uk/chembl/)                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| `tractability`     | Datos de tractabilidad                                            | Ninguno                                           | [Open&nbsp;Targets](https://platform-docs.opentargets.org/target/tractability)                                                                                                                                                                                                                                                                                                                                                                                                                     |
| `pharmacogenetics` | Respuestas farmacogenéticas                                       | `drug_id`                                         | [PharmGKB](https://www.pharmgkb.org/)                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| `expression`       | Datos de expresión génica (por tejidos, órganos y sistemas anatómicos) | `tissue_id`<br/>`anatomical_system`<br/>`organ`   | <ul><li>[ExpressionAtlas](https://www.ebi.ac.uk/gxa/home)</li><li>[HPA](https://www.proteinatlas.org/)</li><li>[GTEx](https://www.gtexportal.org/home/)</li></ul>                                                                                                                                                                                                                                                                                                                                  |
| `depmap`           | Datos de efecto gen&rarr;enfermedad en DepMap.                    | `tissue_id`                                       | [DepMap&nbsp;Portal](https://depmap.org/portal/)                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| `interactions`     | Interacciones proteína&rlarr;proteína                             | `protein_a_id`<br/>`protein_b_id`<br/>`gene_b_id` | <ul><li>[Open&nbsp;Targets](https://platform-docs.opentargets.org/target/molecular-interactions)</li><li>[IntAct](https://platform-docs.opentargets.org/target/molecular-interactions#intact)</li><li>[Signor](https://platform-docs.opentargets.org/target/molecular-interactions#signor)</li><li>[Reactome](https://platform-docs.opentargets.org/target/molecular-interactions#reactome)</li><li>[String](https://platform-docs.opentargets.org/target/molecular-interactions#string)</li></ul> |

`-l` `--limit`  
Limitar el número de resultados, por ejemplo, 10. Predeterminado: Sin límite.     
Nota: No es compatible con los recursos `tractability` y `depmap`.

`-o` `--out`    
Ruta al archivo JSON donde se guardarán los resultados, por ejemplo, path/to/directory/results.json. Predeterminado: Salida estándar.  
Python: `save=True` guardará la salida en el directorio de trabajo actual.

**Argumentos opcionales de filtrado**

`-fd` `--filter_disease` `disease_id`  
Filtrar por ID de enfermedad, por ejemplo, 'EFO_0000274'. *Válido solo para el recurso `drugs`.*

`-fc` `--filter_drug` `drug_id`  
Filtrar por ID de fármaco, por ejemplo, 'CHEMBL1743081'. *Válido solo para el recurso `pharmacogenetics`.*

`-ft` `--filter_tissue` `tissue_id`  
Filtrar por ID de tejido, por ejemplo, 'UBERON_0000473'. *Válido solo para los recursos `expression` y `depmap`.*

`-fa` `--filter_anat_sys`  
Filtrar por sistema anatómico, por ejemplo, 'sistema nervioso'. *Válido solo para el recurso `expression`.*

`-fo` `--filter_organ` `anatomical_system`  
Filtrar por órgano, por ejemplo, 'cerebro'. *Válido solo para el recurso `expression`.*

`-fpa` `--filter_protein_a` `protein_a_id`  
Filtrar por ID de la primera proteína en la interacción, por ejemplo, 'ENSP00000304915'. *Válido solo para el recurso `interactions`.*

`-fpb` `--filter_protein_b` `protein_b_id`  
Filtrar por ID de la segunda proteína en la interacción, por ejemplo, 'ENSP00000379111'. *Válido solo para el recurso `interactions`.*

`-fgb` `--filter_gene_b` `gene_b_id`  
Filtrar por ID de gen de la segunda proteína en la interacción, por ejemplo, 'ENSG00000077238'. *Válido solo para el recurso `interactions`.*

`filters`  
Solo para Python. Un diccionario de filtros, por ejemplo:
```python
{'disease_id': ['EFO_0000274', 'HP_0000964']}

`filter_mode`  
Solo para Python. `filter_mode='or'` combina filtros de diferentes IDs con lógica OR.  
`filter_mode='and'` combina filtros de diferentes IDs con lógica AND (predeterminado).

**Banderas**   
`-csv` `--csv`  
Solo en línea de comandos. Devuelve la salida en formato CSV, en lugar de formato JSON.  
Python: Use `json=True` para devolver la salida en formato JSON.

`-q` `--quiet`   
Solo en línea de comandos. Evita que se muestre la información de progreso.  
Python: Use `verbose=False` para evitar que se muestre la información de progreso. 

`-or` `--or`  
Solo en línea de comandos. Los filtros se combinan con lógica OR. Predeterminado: lógica AND.

`wrap_text`  
Solo para Python. `wrap_text=True` muestra el marco de datos con texto ajustado para facilitar la lectura (predeterminado: False).

### Ejemplos

**Obtenga enfermedades asociadas a un gen específico:**   
```bash
gget opentargets ENSG00000169194 -r diseases -l 1
```
```python
# Python
import gget
gget.opentargets('ENSG00000169194', resource='diseases', limit=1)
```
&rarr; Devuelve la principal enfermedad asociada con el gen ENSG00000169194.

| id            | name               | description                                                           | score            |
|---------------|--------------------|-----------------------------------------------------------------------|------------------|
| EFO_0000274   | atopic&nbsp;eczema | A chronic inflammatory genetically determined disease of the skin ... | 0.66364347241831 |

<br/><br/>

**Obtener medicamentos asociados para un gen específico:**   
```bash
gget opentargets ENSG00000169194 -r drugs -l 2
```
```python
# Python
import gget
gget.opentargets('ENSG00000169194', resource='drugs', limit=2)
```

&rarr; Devuelve los 2 principales medicamentos asociados con el gen ENSG00000169194.

| id            | name         | type     | action_mechanism                    | description                                                  | synonyms                                           | trade_names           | disease_id  | disease_name                  | trial_phase | trial_status | trial_ids     | approved |
|---------------|--------------|----------|-------------------------------------|--------------------------------------------------------------|----------------------------------------------------|-----------------------|-------------|-------------------------------|-------------|--------------|---------------|----------|
| CHEMBL1743081 | TRALOKINUMAB | Antibody | Interleukin&#8209;13&nbsp;inhibitor | Antibody drug with a maximum clinical trial phase of IV ...  | ['CAT-354', 'Tralokinumab']                        | ['Adbry', 'Adtralza'] | EFO_0000274 | atopic&nbsp;eczema            | 4           |              | []            | True     |
| CHEMBL4297864 | CENDAKIMAB   | Antibody | Interleukin&#8209;13&nbsp;inhibitor | Antibody drug with a maximum clinical trial phase of III ... | [ABT-308, Abt-308, CC-93538, Cendakimab, RPC-4046] | []                    | EFO_0004232 | eosinophilic&nbsp;esophagitis | 3           | Recruiting   | [NCT04991935] | False    |

*Note: Los `trial_ids` devueltos son identificadores de [ClinicalTrials.gov](https://clinicaltrials.gov)*

<br/><br/>

**Obtenga datos de trazabilidad para un gen específico:**   
```bash
gget opentargets ENSG00000169194 -r tractability
```
```python
# Python
import gget
gget.opentargets('ENSG00000169194', resource='tractability')
```

&rarr; Devuelve datos de trazabilidad para el gen ENSG00000169194.

| label                 | modality       |
|-----------------------|----------------|
| High-Quality Pocket   | Small molecule |
| Approved Drug         | Antibody       |
| GO CC high conf       | Antibody       |
| UniProt loc med conf  | Antibody       |
| UniProt SigP or TMHMM | Antibody       |

<br/><br/>

**Obtenga respuestas farmacogenéticas para un gen específico:**
```bash
gget opentargets ENSG00000169194 -r pharmacogenetics -l 1
```
```python
# Python
import gget
gget.opentargets('ENSG00000169194', resource='pharmacogenetics', limit=1)
```

&rarr; Devuelve respuestas farmacogenéticas para el gen ENSG00000169194.

| rs_id     | genotype_id       | genotype | variant_consequence_id | variant_consequence_label | drugs                                                                                                                                   | phenotype                                                               | genotype_annotation                                                                                          | response_category | direct_target | evidence_level | source   | literature |
|-----------|-------------------|----------|------------------------|---------------------------|-----------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------|-------------------|---------------|----------------|----------|------------|
| rs1295686 | 5_132660151_T_T,T | TT       | SO:0002073             | no_sequence_alteration    | &nbsp;&nbsp;&nbsp;&nbsp;id&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;name<br/>0&nbsp;&nbsp;None&nbsp;&nbsp;hepatitis&nbsp;vaccines | increased risk for non&#8209;immune response to the hepatitis B vaccine | Patients with the TT genotype may be at increased risk for non-immune response to the hepatitis B vaccine... | efficacy          | False         | 3              | pharmgkb | [21111021] |

*Note: Los identificadores de `literature` devueltos son identificadores de [PMC de Europa](https://europepmc.org/article/med/)*

<br/><br/>

**Obtenga tejidos donde un gen se expresa más:**
```bash
gget opentargets ENSG00000169194 -r expression -l 2
```
```python
# Python
import gget
gget.opentargets('ENSG00000169194', resource='expression', limit=2)
```

&rarr; Devuelve los 2 tejidos principales donde se expresa más el gen ENSG00000169194.

| tissue_id      | tissue_name                           | rna_zscore | rna_value | rna_unit | rna_level | anatomical_systems                                                   | organs                                                 |
|----------------|---------------------------------------|------------|-----------|----------|-----------|----------------------------------------------------------------------|--------------------------------------------------------|
| UBERON_0000473 | testis                                | 5          | 1026      |          | 3         | [reproductive&nbsp;system]                                           | [reproductive&nbsp;organ, reproductive&nbsp;structure] |
| CL_0000542     | EBV&#8209;transformed&nbsp;lymphocyte | 1          | 54        |          | 2         | [hemolymphoid&nbsp;system, immune&nbsp;system, lymphoid&nbsp;system] | [immune organ]                                         |

<br/><br/>

**Obtenga datos sobre el efecto de la enfermedad genética de DepMap para un gen específico:**
```bash
gget opentargets ENSG00000169194 -r depmap
```
```python
# Python
import gget
gget.opentargets('ENSG00000169194', resource='depmap')
```

&rarr; Devuelve datos del efecto de la enfermedad del gen DepMap para el gen ENSG00000169194.

| depmap_id        | expression | effect   | tissue_id      | tissue_name | cell_line_name | disease_cell_line_id | disease_name         | mutation |
|------------------|------------|----------|----------------|-------------|----------------|----------------------|----------------------|----------|
| ACH&#8209;001532 | 0.176323   | 0.054950 | UBERON_0002113 | kidney      | JMU-RTK-2      | None                 | Rhabdoid&nbsp;Cancer | None     |

<br/><br/>

**Obtener interacciones proteína-proteína para un gen específico:**
```bash
gget opentargets ENSG00000169194 -r interactions -l 2
```
```python
# Python
import gget
gget.opentargets('ENSG00000169194', resource='interactions', limit=2)
```

&rarr; Devuelve las 2 interacciones proteína-proteína principales para el gen ENSG00000169194.

| evidence_score | evidence_count | source_db | protein_a_id    | gene_a_id       | gene_a_symbol | role_a                | taxon_a | protein_b_id    | gene_b_id       | gene_b_symbol | role_b                | taxon_b |
|----------------|----------------|-----------|-----------------|-----------------|---------------|-----------------------|---------|-----------------|-----------------|---------------|-----------------------|---------|
| 0.999          | 3              | string    | ENSP00000304915 | ENSG00000169194 | IL13          | unspecified&nbsp;role | 9606    | ENSP00000379111 | ENSG00000077238 | IL4R          | unspecified&nbsp;role | 9606    |
| 0.999          | 3              | string    | ENSP00000304915 | ENSG00000169194 | IL13          | unspecified&nbsp;role | 9606    | ENSP00000360730 | ENSG00000131724 | IL13RA1       | unspecified&nbsp;role | 9606    |

<br/><br/>

**Obtenga interacciones proteína-proteína para un gen específico, filtrando por ID de proteínas y genes:**
```bash
gget opentargets ENSG00000169194 -r interactions -fpa P35225 --filter_gene_b ENSG00000077238
```
```python
# Python
import gget
gget.opentargets('ENSG00000169194', resource='interactions', filters={'protein_a_id': 'P35225', 'gene_b_id': 'ENSG00000077238'})
```

&rarr; Devuelve interacciones proteína-proteína para el gen ENSG00000169194, donde la primera proteína es P35225 **y** el segundo gen es ENSG00000077238:

| evidence_score | evidence_count | source_db | protein_a_id | gene_a_id       | gene_a_symbol | role_a                | taxon_a | protein_b_id | gene_b_id       | gene_b_symbol | role_b                | taxon_b |
|----------------|----------------|-----------|--------------|-----------------|---------------|-----------------------|---------|--------------|-----------------|---------------|-----------------------|---------|
| None           | 3              | reactome  | P35225       | ENSG00000169194 | IL13          | unspecified&nbsp;role | 9606    | P24394       | ENSG00000077238 | IL4R          | unspecified&nbsp;role | 9606    |
| None           | 2              | signor    | P35225       | ENSG00000169194 | IL13          | regulator             | 9606    | P24394       | ENSG00000077238 | IL4R          | regulator&nbsp;target | 9606    |

<br/><br/>

**Obtenga interacciones proteína-proteína para un gen específico, filtrando por ID de proteína o gen:**
```bash
gget opentargets ENSG00000169194 -r interactions -fpa P35225 --filter_gene_b ENSG00000077238 ENSG00000111537 --or -l 5
```
```python
# Python
import gget
gget.opentargets(
    'ENSG00000169194',
    resource='interactions',
    filters={'protein_a_id': 'P35225', 'gene_b_id': ['ENSG00000077238', 'ENSG00000111537']},
    filter_mode='or',
    limit=5
)
```

&rarr; Devuelve interacciones proteína-proteína para el gen ENSG00000169194, donde la primera proteína es P35225 **o** el segundo gen es ENSG00000077238 o ENSG00000111537.
| evidence_score | evidence_count | source_db | protein_a_id    | gene_a_id       | gene_a_symbol | role_a                | taxon_a | protein_b_id    | gene_b_id       | gene_b_symbol | role_b                | taxon_b |
|----------------|----------------|-----------|-----------------|-----------------|---------------|-----------------------|---------|-----------------|-----------------|---------------|-----------------------|---------|
| 0.999          | 3              | string    | ENSP00000304915 | ENSG00000169194 | IL13          | unspecified&nbsp;role | 9606    | ENSP00000379111 | ENSG00000077238 | IL4R          | unspecified&nbsp;role | 9606    |
| 0.961          | 2              | string    | ENSP00000304915 | ENSG00000169194 | IL13          | unspecified&nbsp;role | 9606    | ENSP00000229135 | ENSG00000111537 | IFNG          | unspecified&nbsp;role | 9606    |
| 0.800          | 9              | intact    | P35225          | ENSG00000169194 | IL13          | unspecified&nbsp;role | 9606    | Q14627          | ENSG00000123496 | IL13RA2       | unspecified&nbsp;role | 9606    |
| 0.740          | 6              | intact    | P35225          | ENSG00000169194 | IL13          | unspecified&nbsp;role | 9606    | P78552          | ENSG00000131724 | IL13RA1       | unspecified&nbsp;role | 9606    |
| 0.400          | 1              | intact    | P35225          | ENSG00000169194 | IL13          | unspecified&nbsp;role | 9606    | Q86XT9          | ENSG00000149932 | TMEM219       | stimulator            | 9606    |


    
#### [Más ejemplos](https://github.com/pachterlab/gget_examples)

# Citar    
Si utiliza `gget opentargets` en una publicación, favor de citar los siguientes artículos:

- Luebbert, L., & Pachter, L. (2023). Efficient querying of genomic reference databases with gget. Bioinformatics. [https://doi.org/10.1093/bioinformatics/btac836](https://doi.org/10.1093/bioinformatics/btac836)

- Ochoa D, Hercules A, Carmona M, Suveges D, Baker J, Malangone C, Lopez I, Miranda A, Cruz-Castillo C, Fumis L, Bernal-Llinares M, Tsukanov K, Cornu H, Tsirigos K, Razuvayevskaya O, Buniello A, Schwartzentruber J, Karim M, Ariano B, Martinez Osorio RE, Ferrer J, Ge X, Machlitt-Northen S, Gonzalez-Uriarte A, Saha S, Tirunagari S, Mehta C, Roldán-Romero JM, Horswell S, Young S, Ghoussaini M, Hulcoop DG, Dunham I, McDonagh EM. The next-generation Open Targets Platform: reimagined, redesigned, rebuilt. Nucleic Acids Res. 2023 Jan 6;51(D1):D1353-D1359. doi: [10.1093/nar/gkac1046](https://doi.org/10.1093/nar/gkac1046). PMID: 36399499; PMCID: PMC9825572.
