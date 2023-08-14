> Parámetros de Python són iguales a los parámetros largos (`--parámetro`) de Terminal, si no especificado de otra manera. Banderas son parámetros de verdadero o falso (True/False) en Python. El manuál para cualquier modulo de gget se puede llamar desde la Terminal con la bandera `-h` `--help`.  
## gget cellxgene 🍱  
Regresa data de la base de datos [CZ CELLxGENE Discover](https://cellxgene.cziscience.com/) usando [CZ CELLxGENE Discover Census](https://github.com/chanzuckerberg/cellxgene-census).  
Regresa: Un objeto del tipo AnnData que contiene la matriz de recuentos de genes y los metadatos de los datos de single-cell RNA sequencing de los tejidos/genes/etc. definidos.  

Antes de usar `gget cellxgene` por primera vez, corre `gget setup cellxgene` / `gget.setup("cellxgene")` (ver también [`gget setup`](setup.md)).  

**Parámetros optionales**  
`-s` `--species`  
'homo_sapiens' o 'mus_musculus'. Por defecto: 'homo_sapiens'.  

`-g` `--gene`  
Str o lista de genes de interés o ID(s) del tipo Ensembl. Por defecto: None (ninguno).  
Atención: Utilice la bandera `-e / --ensembl` (Python: `ensembl=True`) cuando ingresas ID(s) del tipo Ensembl.    
Ver https://cellxgene.cziscience.com/gene-expression para ejemplos de genes.  

`-cv` `--census_version`  
Versión del CZ CELLxGENE Discover Census (str), p. "2023-05-15", o "latest" (la ultima) o "stable" (estable). Por defecto: "stable" (estable).  

`-cn` `--column_names`  
Lista de columnas de metadatos a regresar (almacenadas en AnnData.obs).  
Por defecto: ['dataset_id', 'assay', 'suspension_type', 'sex', 'tissue_general', 'tissue', 'cell_type']   
Para más opciones, ver: https://api.cellxgene.cziscience.com/curation/ui/#/ -> 'Schemas' -> 'dataset'  

`-o` `--out`   
Ruta al archivo para guardar el objeto AnnData de tipo .h5ad (o .csv con bandera `-mo / --meta_only`).  
¡Requerido cuando se usa desde la Terminal!  

**Banderas**  
`-e` `--ensembl`  
Usa esta bandera si `gene` se ingresa como ID del tipo Ensembl.    

`-mo` `--meta_only`  
Solo regresa la tabla (Dataframe) con metadatos (corresponde a AnnData.obs).  

`-q` `--quiet`   
Solo para la Terminal. Impide la informacion de progreso de ser exhibida durante la corrida.  
Para Python, usa `verbose=False` para imipidir la informacion de progreso de ser exhibida durante la corrida.  

**Parámetros optionales correspondientes a los atributos de metadatos de CZ CELLxGENE Discover**  
`--tissue`  
Str o lista de tejido(s), p. ['lung', 'blood']. Por defecto: None.  
Ver https://cellxgene.cziscience.com/gene-expression para ejemplos de tejidos.  

`--cell_type`  
Str o lista de tipo(s) de célula(s), p. ['mucus secreting cell', 'neuroendocrine cell']. Por defecto: None.  
Ver https://cellxgene.cziscience.com/gene-expression y seleccione un tejido para ejemplos de tipos de células.  

`--development_stage`  
Str o lista de etapa(s) de desarrollo. Por defecto: None.  

`--disease`  
Str o lista de enfermedade(s). Por defecto: None.  

`--sex`  
Str o lista de sexo(s), p. 'female' (femenina). Por defecto: None.  

`--dataset_id`  
Str o lista de CELLxGENE ID(s). Por defecto: None.  

`--tissue_general_ontology_term_id`  
Str o lista de tejido(s) del tipo high-level UBERON ID. Por defecto: None.  
Tejidos y sus IDs de UBERON se enumeran [aquí](https://github.com/chanzuckerberg/single-cell-data-portal/blob/9b94ccb0a2e0a8f6182b213aa4852c491f6f6aff/backend/wmg/data/tissue_mapper.py).  

`--tissue_general`  
Str o lista de tejido(s) del tipo high-level. Por defecto: None.  
Tejidos y sus IDs de UBERON se enumeran [aquí](https://github.com/chanzuckerberg/single-cell-data-portal/blob/9b94ccb0a2e0a8f6182b213aa4852c491f6f6aff/backend/wmg/data/tissue_mapper.py).  

`--tissue_ontology_term_id`  
Str o lista de ID(s) de 'tissue ontology term' como están definidos en el [esquema de datos del CELLxGENE](https://github.com/chanzuckerberg/single-cell-curation/tree/main/schema). Por defecto: None.   

`--assay_ontology_term_id`  
Str o lista de ID(s) de 'assay ontology term' como están definidos en el [esquema de datos del CELLxGENE](https://github.com/chanzuckerberg/single-cell-curation/tree/main/schema). Por defecto: None.  

`--assay`  
Str o lista de 'assays' (métodos) como están definidos en el [esquema de datos del CELLxGENE](https://github.com/chanzuckerberg/single-cell-curation/tree/main/schema). Por defecto: None.  

`--cell_type_ontology_term_id`  
Str o lista de ID(s) de 'celltype ontology term' como están definidos en el [esquema de datos del CELLxGENE](https://github.com/chanzuckerberg/single-cell-curation/tree/main/schema). Por defecto: None.  

`--development_stage_ontology_term_id`   
Str o lista de ID(s) de 'development stage ontology term' como están definidos en el [esquema de datos del CELLxGENE](https://github.com/chanzuckerberg/single-cell-curation/tree/main/schema). Por defecto: None.  

`--disease_ontology_term_id`  
Str o lista de ID(s) de 'disease ontology term' como están definidos en el [esquema de datos del CELLxGENE](https://github.com/chanzuckerberg/single-cell-curation/tree/main/schema). Por defecto: None.  

`--donor_id`  
Str o lista de ID(s) de 'donor' (donador) como están definidos en el [esquema de datos del CELLxGENE](https://github.com/chanzuckerberg/single-cell-curation/tree/main/schema). Por defecto: None.  

`--self_reported_ethnicity_ontology_term_id`  
Str o lista de ID(s) de 'self-reported ethnicity ontology' como están definidos en el [esquema de datos del CELLxGENE](https://github.com/chanzuckerberg/single-cell-curation/tree/main/schema). Por defecto: None.  

`--self_reported_ethnicity`  
Str o lista de etnias autoinformadas como están definidas en el [esquema de datos del CELLxGENE](https://github.com/chanzuckerberg/single-cell-curation/tree/main/schema). Por defecto: None.  

`--sex_ontology_term_id`  
Str o lista de ID(s) de 'sex ontology' como están definidos en el [esquema de datos del CELLxGENE](https://github.com/chanzuckerberg/single-cell-curation/tree/main/schema). Por defecto: None.  

`--suspension_type`  
Str o lista de tipo(s) de suspensión como están definidos en el [esquema de datos del CELLxGENE](https://github.com/chanzuckerberg/single-cell-curation/tree/main/schema). Por defecto: None.  

  
### Por ejemplo
```bash
gget cellxgene --gene ACE2 ABCA1 SLC5A1 --tissue lung --cell_type 'mucus secreting cell' 'neuroendocrine cell' -o example_adata.h5ad
```
```python
# Python
adata = gget.cellxgene(
    gene = ["ACE2", "ABCA1", "SLC5A1"],
    tissue = "lung",
    cell_type = ["mucus secreting cell", "neuroendocrine cell"]
)
adata
```
&rarr; Regresa un objeto AnnData que contiene la matriz de recuentos de scRNAseq de los genes ACE2, ABCA1 y SLC5A1 en 3322 células secretoras de moco y neuroendocrinas del pulmón humano y sus metadatos correspondientes.

<br/><br/>

Obtene solo los metadatos (corresponde a AnnData.obs):  
```bash
gget cellxgene --meta_only --gene ENSMUSG00000015405 --ensembl --tissue lung --species mus_musculus -o example_meta.csv
```
```python
# Python
df = gget.cellxgene(
    meta_only = True,
    gene = "ENSMUSG00000015405",
    ensembl = True,
    tissue = "lung",  
    species = "mus_musculus"
)
df
```
&rarr; Regresa solo los metadatos de conjuntos de datos de células de pulmón de ratón que contienen recuentos para ENSMUSG00000015405 (ACE2).  
