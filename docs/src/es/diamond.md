[<kbd> Ver el codigo fuente de la pagina en GitHub </kbd>](https://github.com/pachterlab/gget/blob/main/docs/src/es/diamond.md)

> Parámetros de Python són iguales a los parámetros largos (`--parámetro`) de Terminal, si no especificado de otra manera. Las banderas son parámetros de verdadero o falso (True/False) en Python. El manuál para cualquier modulo de gget se puede llamar desde la Terminal con la bandera `-h` `--help`.  
# gget diamond 💎
Alinee múltiples proteínas o secuencias de ADN traducidas usando [DIAMOND](https://www.nature.com/articles/nmeth.3176) (DIAMOND es similar a BLAST, pero este es un cálculo local).       
Produce: Resultados en formato JSON (Terminal) o Dataframe/CSV (Python).  

**Parámetro posicional**  
`query`  
Secuencia(s) (str o lista) de aminoácidos, o una ruta a un archivo tipo FASTA.    

**Parámetro requerido**  
`-ref` `--reference`  
Secuencias de aminoácidos de referencia (str o lista), o una ruta a un archivo tipo FASTA.  

**Parámetros optionales**  
`-db` `--diamond_db`  
Ruta para guardar la base de datos DIAMOND creada a partir de `reference` (str).  
Por defecto: None -> El archivo de base de datos DIAMOND temporal se eliminará después de la alineación o se guardará en `out` si se proporciona `out`.  

`-s` `--sensitivity`  
Sensibilidad de la alineación (str). Por defecto: "very-sensitive" (muy sensible).  
Uno de los siguientes: fast, mid-sensitive, sensitive, more-sensitive, very-sensitive, or ultra-sensitive.    

`-t` `--threads`  
Número de hilos de procesamiento utilizados (int). Por defecto: 1.  

`-db` `--diamond_binary`  
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
# !!! Asegúrese de enumerar primero el argumento posicional aquí para que no se agregue como secuencia de referencia
gget diamond GGETISAWESQME ELVISISALIVE LQVEFRANKLIN PACHTERLABRQCKS -ref GGETISAWESQMEELVISISALIVELQVEFRANKLIN PACHTERLABRQCKS
```
```python
# Python
gget.diamond(["GGETISAWESQME", "ELVISISALIVE", "LQVEFRANKLIN", "PACHTERLABRQCKS"], reference=["GGETISAWESQMEELVISISALIVELQVEFRANKLIN", "PACHTERLABRQCKS"])
```
&rarr; Produce los resultados de la alineación en formato JSON (Terminal) o Dataframe/CSV:  

|query_accession|subject_accession|identity_percentage|query_seq_length|subject_seq_length|length|mismatches|gap_openings|query_start|query_end|subject_start|subject_end|e-value |bit_score|
|---------------|-----------------|-------------------|----------------|------------------|------|----------|------------|-----------|---------|-------------|-----------|--------|---------|
|Seq0           |Seq0             |100                |13              |37                |13    |0         |0           |1          |13       |1            |13         |2.82e-09|30.8     |
|Seq2           |Seq0             |100                |12              |37                |12    |0         |0           |1          |12       |26           |37         |4.35e-08|27.7     |
|Seq3           |Seq1             |100                |15              |15                |15    |0         |0           |1          |15       |1            |15         |2.01e-11|36.2     |


#### [Màs ejemplos](https://github.com/pachterlab/gget_examples)

# Citar    
Si utiliza `gget diamond` en una publicación, favor de citar los siguientes artículos:

- Luebbert, L., & Pachter, L. (2023). Efficient querying of genomic reference databases with gget. Bioinformatics. [https://doi.org/10.1093/bioinformatics/btac836](https://doi.org/10.1093/bioinformatics/btac836)

- Buchfink, B., Xie, C. & Huson, D. Fast and sensitive protein alignment using DIAMOND. Nat Methods 12, 59–60 (2015). [https://doi.org/10.1038/nmeth.3176](https://doi.org/10.1038/nmeth.3176)
