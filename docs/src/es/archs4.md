> Parámetros de Python són iguales a los parámetros largos (`--parámetro`) de Terminal, si no es especificado de otra manera. Las banderas son designadas como cierto o falso (True/False) en Python. El manuál para cualquier modulo de gget se puede obtener desde Terminal con la bandera `-h` `--help`.  
## gget archs4 🐁
Encuentra los genes más correlacionados a un gen de interés, o bién, encuentra los tejidos donde un gen se expresa usando la base de datos [ARCHS4](https://maayanlab.cloud/archs4/).  
Produce: Resultados en formato JSON (Terminal) o Dataframe/CSV (Python).  

**Parámetro posicional**  
`gene`  
Nombre corto (símbolo del gen) del gen de interés, p. ej. STAT4.  
Alternativamente: usa la bandera `--ensembl` para ingresar un ID del tipo Ensembl, p. ej. ENSG00000138378.  

**Parámetros optionales**  
 `-w` `--which`  
'correlation' (correlación; esto se usa por defecto) o 'tissue' (tejido).  
'correlation' regresa una tabla que contiene los 100 genes más correlacionados con el gen de interés. La correlación de Pearson se calcula sobre todas las muestras y tejidos en [ARCHS4](https://maayanlab.cloud/archs4/).  
'tissue' regresa un atlas de expresión tisular calculado sobre todas las muestras humanas o de ratón (según lo definido usando el parámetro `--species` (especies)) en [ARCHS4](https://maayanlab.cloud/archs4/).  

`-s` `--species`  
'human' (humano; esto se usa por defecto) o 'mouse' (ratón).   
Define si usar muestras humanas o de ratón de [ARCHS4](https://maayanlab.cloud/archs4/).  
(Solo aplicable para el atlas de expresión tisular.)  

`-o` `--out`   
Ruta al archivo en el que se guardarán los resultados, p. ej. ruta/al/directorio/resultados.csv (o .json). Por defecto: salida estándar (STDOUT).  
Para Python, usa `save=True` para guardar los resultados en el directorio de trabajo actual.  
  
**Banderas**   
`-e` `--ensembl`  
Usa esta bandera si `gene` se ingresa como ID del tipo Ensembl.   

`-csv` `--csv`  
Solo para la Terminal. Regresa los resultados en formato CSV.    
Para Python, usa `json=True` para regresar los resultados en formato JSON.    

`-q` `--quiet`   
Solo para la Terminal. Impide la informacion de progreso de ser exhibida durante la corrida.  
Para Python, usa `verbose=False` para imipidir la informacion de progreso de ser exhibida durante la corrida.  
  
  
### Por ejemplo
```bash
gget archs4 ACE2
```
```python
# Python
gget.archs4("ACE2")
```
&rarr; Regresa los 100 genes más correlacionados con el gen ACE2:  

| gene_symbol     | pearson_correlation     |
| -------------- |-------------------------| 
| SLC5A1 | 0.579634 | 	
| CYP2C18 | 0.576577 | 	
| . . . | . . . | 	

<br/><br/>

```bash
gget archs4 -w tissue ACE2
```
```python
# Python
gget.archs4("ACE2", which="tissue")
```
&rarr; Regresa la expresión tisular de ACE2 (por defecto, se utilizan datos humanos):  

| id     | min     | q1 |  median | q3 | max |
| ------ |--------| ------ |--------| ------ |--------| 
| System.Urogenital/Reproductive System.Kidney.RENAL CORTEX | 0.113644 | 8.274060 | 9.695840 | 10.51670 | 11.21970 |
| System.Digestive System.Intestine.INTESTINAL EPITHELIAL CELL | 0.113644 | 	5.905560 | 9.570450 | 13.26470 | 13.83590 | 
| . . . | . . . | . . . | . . . | . . . | . . . |

<br/><br/>
Consulte [este tutorial](https://davetang.org/muse/2023/05/16/check-where-a-gene-is-expressed-from-the-command-line/) de Dave Tang, quien escribió un script R para crear esta visualización con los resultados de `gget archs4` en formato JSON:  

![image](https://github.com/pachterlab/gget/assets/56094636/f2a34a9e-beaa-45a5-a678-d38399dd3017)


#### [Más ejemplos](https://github.com/pachterlab/gget_examples)  
