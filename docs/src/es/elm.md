> Parámetros de Python són iguales a los parámetros largos (`--parámetro`) de Terminal, si no especificado de otra manera. Banderas son parámetros de verdadero o falso (True/False) en Python. El manuál para cualquier modulo de gget se puede llamar desde la Terminal con la bandera `-h` `--help`.  
## gget elm 🎭
Prediga localmente motivos lineales eucarióticos (ELMs) a partir de una secuencia de aminoácidos o UniProt ID utilizando datos de la [base de datos ELM](http://elm.eu.org/).  
Produce: Resultados en formato JSON (Terminal) o Dataframe/CSV (Python). Este módulo devuelve dos tipos de resultados (ver ejemplos).   

Antes de usar `gget elm` por primera vez, ejecute `gget setup elm` / `gget.setup("elm")` una vez (consulte también [`gget setup`](setup.md)).   

**Parámetro posicional**  
`sequence`  
Secuencia de aminoácidos o ID de Uniprot (str).  
Al proporcionar una ID de Uniprot, use la bandera `--uniprot` (Python: `uniprot==True`).  

**Parámetros optionales**  
`-s` `sensitivity`  
Sensibilidad de la alineación DIAMOND (str). Por defecto: "very-sensitive" (muy sensible).  
Uno de los siguientes: fast, mid-sensitive, sensitive, more-sensitive, very-sensitive, or ultra-sensitive.  

`-t` `threads`  
Número de hilos de procesamiento utilizados en la alineación de secuencias con DIAMOND (int). Por defecto: 1.  

`-bin` `diamond_binary`  
Ruta al binario DIAMOND (str). Por defecto: None -> Utiliza el binario DIAMOND instalado automáticamente con `gget`.  

`-o` `--out`   
Ruta al archivo en el que se guardarán los resultados (str), p. ej. "ruta/al/directorio". Por defecto: salida estándar (STDOUT); los archivos temporales se eliminan.  

**Banderas**  
`-u` `--uniprot`  
Use esta bandera cuando `sequence` es un ID de Uniprot en lugar de una secuencia de aminoácidos.      

`-csv` `--csv`  
Solo para Terminal. Produce los resultados en formato CSV.    
Para Python, usa `json=True` para producir los resultados en formato JSON.  

`-q` `--quiet`   
Solo para Terminal. Impide la información de progreso de ser exhibida durante la ejecución del programa.  
Para Python, usa `verbose=False` para impedir la información de progreso de ser exhibida durante la ejecución del programa.  

### Ejemplo
Encuentre ELM en una secuencia de aminoácidos:  
```bash
gget setup elm          # Descarga/actualiza la base de datos ELM local
gget elm -o gget_elm_results LIAQSIGQASFV
```
```python
# Python
gget.setup(“elm”)      # Descarga/actualiza la base de datos ELM local
ortholog_df, regex_df = gget.elm("LIAQSIGQASFV")
```
  
Encuentre ELM que proporcionen un ID de UniProt como entrada: 
```bash
gget setup elm          # Descarga/actualiza la base de datos ELM local
gget elm -o gget_elm_results --uniprot Q02410
```
```python
# Python
gget.setup(“elm”)      # Descarga/actualiza la base de datos ELM local
ortholog_df, regex_df = gget.elm("Q02410", uniprot=True)
```
&rarr; Produce dos resultados con información extensa sobre ELMs asociados con proteínas ortólogas y motivos encontrados en la secuencia de entrada directamente en función de sus expresiones regex:  

Ortholog data frame:  
|Ortholog_UniProt_ID|ProteinName|class_accession|ELMIdentifier  |FunctionalSiteName                   |Description                                                                                                                              |Probability|Methods                                                                                                             |Organism    |motif_start_in_subject|motif_end_in_subject|…  |
|:-----------------:|:---------:|:-------------:|:-------------:|:-----------------------------------:|:---------------------------------------------------------------------------------------------------------------------------------------:|:---------:|:------------------------------------------------------------------------------------------------------------------:|:----------:|:--------------------:|:------------------:|:-:|
|Q02410             |APBA1_HUMAN|ELME000357     |LIG_CaMK_CASK_1|CASK CaMK domain binding ligand motif|Motif that mediates binding to the calmodulin-dependent protein kinase (CaMK) domain of the peripheral plasma membrane protein CASK/Lin2.|1,12E+06   |comigration in gel electrophoresis; far western blotting; glutathione s tranferase tag; mutation analysis; pull down|Homo sapiens|376                   |382                 |…  |
|Q02410             |APBA1_HUMAN|ELME000091     |LIG_PDZ_Class_2|PDZ domain ligands                   |The C-terminal class 2 PDZ-binding motif is classically represented by a pattern such as                                                 |7,89E+06   |mutation analysis; nuclear magnetic resonance; two hybrid                                                           |Homo sapiens|832                   |837                 |…  |

Regex data frame:  
|Instance_accession|ELMIdentifier     |FunctionalSiteName             |ELMType|Description                                                                                                                                            |Instances (Matched Sequence)|motif_start_in_query|motif_end_in_query|ProteinName|Organism                      |…  |
|:----------------:|:----------------:|:-----------------------------:|:-----:|:-----------------------------------------------------------------------------------------------------------------------------------------------------:|:--------------------------:|:------------------:|:----------------:|:---------:|:----------------------------:|:-:|
|ELME000321        |CLV_C14_Caspase3-7|Caspase cleavage motif         |CLV    |Caspase-3 and Caspase-7 cleavage site.                                                                                                                 |ERSDG                       |239                 |244               |IL16_MOUSE |Mus musculus                  |…  |
|ELME000102        |CLV_NRD_NRD_1     |NRD cleavage site              |CLV    |N-Arg dibasic convertase (NRD/Nardilysin) cleavage site.                                                                                               |RRA                         |142                 |145               |PDYN_RAT   |Rattus norvegicus             |…  |
|ELME000100        |CLV_PCSK_PC1ET2_1 |PCSK cleavage site             |CLV    |NEC1/NEC2 cleavage site.                                                                                                                               |KRD                         |341                 |344               |COLI_MOUSE |Mus musculus                  |…  |
|ELME000146        |CLV_PCSK_SKI1_1   |PCSK cleavage site             |CLV    |Subtilisin/kexin isozyme-1 (SKI1) cleavage site.                                                                                                       |RLLTA                       |825                 |830               |SRBP2_HUMAN|Homo sapiens                  |…  |
|ELME000231        |DEG_APCC_DBOX_1   |APCC-binding Destruction motifs|DEG    |An RxxL-based motif that binds to the Cdh1 and Cdc20 components of APC/C thereby targeting the protein for destruction in a cell cycle dependent manner|SRVKLNIVR                   |732                 |741               |ACM1_YEAST |Saccharomyces cerevisiae S288c|…  |
|…                 |…                 |…                              |…      |…                                                                                                                                                      |…                           |…                   |…                 |…          |…                             |…  |


(Los motivos que aparecen en muchas especies diferentes pueden parecer repetidos, pero todas las filas deben ser únicas.)

#### [Màs ejemplos](https://github.com/pachterlab/gget_examples)  
