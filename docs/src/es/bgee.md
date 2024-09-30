
> Parámetros de Python són iguales a los parámetros largos (`--parámetro`) de Terminal, si no especificado de otra manera. Banderas son parámetros de verdadero o falso (True/False) en Python. El manuál para cualquier modulo de gget se puede llamar desde la Terminal con la bandera `-h` `--help`.  
# gget bgee 🐝

Obtenga datos de ortología y expresión genética de [Bgee](https://www.bgee.org/) utilizando IDs de Ensembl.  
Resultado: JSON/CSV (línea de comandos) o marco de datos (Python).  

> Si estás interesado específicamente en datos de expresión génica humana, considera usar [gget opentargets](./opentargets.md) o [gget archs4](./archs4.md) en su lugar. **gget bgee** tiene menos datos, pero admite más especies.

Este módulo fue escrito por [Sam Wagenaar](https://github.com/techno-sam).

**Argumento posicional**  
`ens_id`  
ID de gen Ensembl, por ejemplo, ENSG00000169194 o ENSSSCG00000014725.

**Argumentos requeridos**  
`-t` `--type`  
Tipo de datos a obtener. Opciones: `orthologs`, `expression`.  

**Argumentos opcionales**  
`-o` `--out`  
Ruta al archivo JSON donde se guardarán los resultados, por ejemplo, path/to/directory/results.json. Por defecto: Salida estándar.

**Banderas**  
`-csv` `--csv`  
Solo en línea de comandos. Devuelve la salida en formato CSV, en lugar de formato JSON.  
Python: Usa `json=True` para devolver la salida en formato JSON.

`-q` `--quiet`  
Solo en línea de comandos. Evita que se muestre la información de progreso.  
Python: Usa `verbose=False` para evitar que se muestre la información de progreso.
  
### Ejemplos

**Obtener ortólogos para un gen**

```bash
gget bgee ENSSSCG00000014725 -t orthologs
```
```python
import gget
gget.bgee("ENSSSCG00000014725", type="orthologs")
```

&rarr; Devuelve ortólogos para el gen con el ID de Ensembl ENSSSCG00000014725.

| gene_id            | gene_name    | species_id | genus   | species    |
|--------------------|--------------|------------|---------|------------|
| 734881             | hbb1         | 8355       | Xenopus | laevis     |
| ENSFCAG00000038029 | LOC101098159 | 9685       | Felis   | catus      |
| ENSBTAG00000047356 | LOC107131172 | 9913       | Bos     | taurus     |
| ENSOARG00000019163 | LOC101105437 | 9940       | Ovis    | aries      |
| ENSXETG00000025667 | hbg1         | 8364       | Xenopus | tropicalis |
| ...                | ...          | ...        | ...     | ...        |

<br/><br/>

**Obtener datos de expresión génica para un gen**

```bash
gget bgee ENSSSCG00000014725 -t expression
```
```python
import gget
gget.bgee("ENSSSCG00000014725", type="expression")
```

&rarr; Devuelve datos de expresión génica para el gen con el ID de Ensembl ENSSSCG00000014725.

| anat_entity_id | anat_entity_name            | score | score_confidence | expression_state |
|----------------|-----------------------------|-------|------------------|------------------|
| UBERON:0000178 | blood                       | 99.98 | high             | expressed        |
| UBERON:0002106 | spleen                      | 99.96 | high             | expressed        |
| UBERON:0002190 | subcutaneous adipose tissue | 99.70 | high             | expressed        |
| UBERON:0005316 | endocardial endothelium     | 99.61 | high             | expressed        |
| UBERON:0002107 | liver                       | 99.27 | high             | expressed        |
| ...            | ...                         | ...   | ...              | ...              |

    
#### [Más ejemplos](https://github.com/pachterlab/gget_examples)

# Citar    
Si utiliza `gget bgee` en una publicación, favor de citar los siguientes artículos:

- Luebbert, L., & Pachter, L. (2023). Efficient querying of genomic reference databases with gget. Bioinformatics. [https://doi.org/10.1093/bioinformatics/btac836](https://doi.org/10.1093/bioinformatics/btac836)

- Frederic B Bastian, Julien Roux, Anne Niknejad, Aurélie Comte, Sara S Fonseca Costa, Tarcisio Mendes de Farias, Sébastien Moretti, Gilles Parmentier, Valentine Rech de Laval, Marta Rosikiewicz, Julien Wollbrett, Amina Echchiki, Angélique Escoriza, Walid H Gharib, Mar Gonzales-Porta, Yohan Jarosz, Balazs Laurenczy, Philippe Moret, Emilie Person, Patrick Roelli, Komal Sanjeev, Mathieu Seppey, Marc Robinson-Rechavi (2021). The Bgee suite: integrated curated expression atlas and comparative transcriptomics in animals. Nucleic Acids Research, Volume 49, Issue D1, 8 January 2021, Pages D831–D847, [https://doi.org/10.1093/nar/gkaa793](https://doi.org/10.1093/nar/gkaa793)
