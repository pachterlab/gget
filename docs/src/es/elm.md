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

`-db` `diamond_binary`  
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
```bash
gget setup elm          # Descarga/actualiza la base de datos ELM local
gget elm -o gget_elm_results LIAQSIGQASFV
gget elm -o gget_elm_results --uniprot Q02410
```
```python
# Python
gget.setup(“elm”)      # Descarga/actualiza la base de datos ELM local
ortholog_df, regex_df = gget.elm(“LIAQSIGQASFV”)
ortholog_df, regex_df = gget.elm(“Q02410”, uniprot=True)
```
&rarr; Produce dos resultados con información extensa sobre ELMs asociados con proteínas ortólogas y motivos encontrados en la secuencia de entrada directamente en función de sus expresiones regex.  

#### [Màs ejemplos](https://github.com/pachterlab/gget_examples)  
