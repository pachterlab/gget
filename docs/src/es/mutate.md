> Parámetros de Python són iguales a los parámetros largos (`--parámetro`) de Terminal, si no especificado de otra manera. Banderas son parámetros de verdadero o falso (True/False) en Python. El manuál para cualquier modulo de gget se puede llamar desde la Terminal con la bandera `-h` `--help`.   
## gget mutate 🧟
Recibe secuencias de nucleótidos y mutaciones (en [anotación de mutación estándar](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1867422/)) y devuelve versiones mutadas de las secuencias según las mutaciones proporcionadas.  
Formato de devolución: Guarda las secuencias mutadas en formato FASTA (o devuelve una lista que contiene las secuencias mutadas si `out=None`).  

**Parámetro posicional**  
`sequences`   
Ruta al archivo FASTA que contiene las secuencias a mutar, por ejemplo, 'path/to/seqs.fa'.  
Los identificadores de secuencia que siguen al carácter '>' deben corresponder a los identificadores en la columna seq_ID de mutations.  
NOTA: Solo se utilizará la cadena que sigue al '>' hasta el primer espacio o punto como identificador de secuencia. -> Se ignorarán los números de versión de los IDs de Ensembl.

Ejemplo de formato del archivo FASTA:  
```
>seq1 (o ENSG00000106443)  
ACTGCGATAGACT  
>seq2  
AGATCGCTAG
```

Alternativamente: Proporcione secuencia(s) como una cadena o lista, por ejemplo, 'AGCTAGCT'.

**Otros parámetros requeridos**  
`-m` `--mutations`  
Ruta al archivo csv o tsv (por ejemplo, 'path/to/mutations.csv') o marco de datos (objeto DataFrame) que contiene información sobre las mutaciones en el siguiente formato (la columna 'notas' no es necesaria):  

| mutation         | mut_ID | seq_ID | notas |
|------------------|--------|--------|-|
| c.2C>T           | mut1   | seq1   | -> Aplicar la mutación 1 a la secuencia 1 |
| c.9_13inv        | mut2   | seq2   | -> Aplicar la mutación 2 a la secuencia 2 |
| c.9_13inv        | mut2   | seq4   | -> Aplicar la mutación 2 a la secuencia 4 |
| c.9_13delinsAAT  | mut3   | seq4   | -> Aplicar la mutación 3 a la secuencia 4 |
| ...              | ...    | ...    | |

'mutation' = Columna que contiene las mutaciones a realizar, escritas en anotación de mutación estándar  
'mut_ID' = Columna que contiene el identificador de cada mutación  
'seq_ID' = Columna que contiene los identificadores de las secuencias a mutar (deben corresponder a la cadena que sigue al carácter '>' en el archivo FASTA de 'sequences'; NO incluir espacios ni puntos)  

Alternativamente: Mutación(es) de entrada como una cadena o lista, por ejemplo, 'c.2C>T'.  
Si se proporciona una lista, el número de mutaciones debe ser igual al número de secuencias de entrada.  
Para uso desde el terminal (bash): Encierre las anotaciones de mutación individuales entre comillas para evitar errores.  

**Parámetros opcionales**  
`-k` `--k`  
Longitud de las secuencias que flanquean la mutación. Por defecto: 30.  
Si k > longitud total de la secuencia, se mantendrá toda la secuencia.

`-mc` `--mut_column`  
Nombre de la columna que contiene las mutaciones a realizar en `mutations`. Por defecto: 'mutation'.  

`-mic` `--mut_id_column`  
Nombre de la columna que contiene los IDs de cada mutación en `mutations`. Por defecto: 'mut_ID'.  

`-sic` `--seq_id_column`  
Nombre de la columna que contiene los IDs de las secuencias a mutar en `mutations`. Por defecto: 'seq_ID'.  

`-o` `--out`   
Ruta al archivo FASTA de salida que contiene las secuencias mutadas, por ejemplo, 'path/to/output_fasta.fa'.  
Por defecto: `None` -> devuelve una lista de las secuencias mutadas a la salida estándar.  
Los identificadores (después del '>') de las secuencias mutadas en el FASTA de salida serán '>[seq_ID]_[mut_ID]'.

**Flags**  
`-q` `--quiet`   
Solo para Terminal. Impide la información de progreso de ser exhibida durante la ejecución del programa.  
Para Python, usa `verbose=False` para imipidir la información de progreso de ser exhibida durante la ejecución del programa.  

### Ejemplos
```bash
gget mutate ATCGCTAAGCT -m 'c.4G>T'
```
```python
# Python
gget.mutate("ATCGCTAAGCT", "c.4G>T")
```
&rarr; Devuelve ATCTCTAAGCT.  

<br/><br/>

**Lista de secuencias con una mutación para cada secuencia proporcionada en una lista:**  
```bash
gget mutate ATCGCTAAGCT TAGCTA -m 'c.4G>T' 'c.1_3inv' -o mut_fasta.fa
```
```python
# Python
gget.mutate(["ATCGCTAAGCT", "TAGCTA"], ["c.4G>T", "c.1_3inv"], out="mut_fasta.fa")
```
&rarr; Guarda el archivo 'mut_fasta.fa' que contiene:  
```
>seq1_mut1  
ATCTCTAAGCT  
>seq2_mut2  
GATCTA
```

<br/><br/>

**Una mutación aplicada a varias secuencias con k ajustado:**  
```bash
gget mutate ATCGCTAAGCT TAGCTA -m 'c.1_3inv' -k 3
```
```python
# Python
gget.mutate(["ATCGCTAAGCT", "TAGCTA"], "c.1_3inv", k=3)
```
&rarr; Devuelve ['CTAGCT', 'GATCTA'].  
