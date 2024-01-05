> Parámetros de Python són iguales a los parámetros largos (`--parámetro`) de Terminal, si no especificado de otra manera. Banderas son parámetros de verdadero o falso (True/False) en Python. El manuál para cualquier modulo de gget se puede llamar desde la Terminal con la bandera `-h` `--help`.  
## gget ref 📖
Obtenga enlaces FTP y sus respectivos metadatos (o use la bandera `ftp` para regresar solo los enlaces) para referenciar genomas y anotaciones de [Ensembl](https://www.ensembl.org/).  
Regresa: Resultados en formato JSON.  

**Mientras Ensembl está en el proceso de actualizar su base de datos a una nueva versión, es posible que reciba un ERROR 404.**   
Si este es el caso, especifique una versión anterior de Ensembl usando el argumento `release`.  
Ejemplo: `gget ref -r 110 human` (Python: `gget.ref("human", release=110)`)

**Parámetro posicional**  
`species`  
La especie por la cual que se buscará los FTP en el formato género_especies, p. ej. homo_sapiens.  
Nota: No se requiere cuando se llama a la bandera `--list_species`.    
Accesos directos: 'human', 'mouse'  

**Parámetros optionales**  
`-w` `--which`  
Define qué resultados devolver. Por defecto: 'all' -> Regresa todos los resultados disponibles.   
Las entradas posibles son uno solo o una combinación de las siguientes (como lista separada por comas):   
'gtf' - Regresa la anotación (GTF).  
'cdna' - Regresa el transcriptoma (cDNA).  
'dna' - Regresa el genoma (DNA).  
'cds' - Regresa las secuencias codificantes correspondientes a los genes Ensembl. (No contiene UTR ni secuencia intrónica).  
'cdrna' - Regresa secuencias de transcripción correspondientes a genes de ARN no codificantes (ncRNA).   
'pep' - Regresa las traducciones de proteínas de los genes Ensembl.    

`-r` `--release`  
Define el número de versión de Ensembl desde el que se obtienen los archivos, p. ej. 104. Default: latest Ensembl release.  

`-o` `--out`   
Ruta al archivo en el que se guardarán los resultados, p. ej. ruta/al/directorio/resultados.json. Por defecto: salida estándar (STDOUT).  
Para Python, usa `save=True` para guardar los resultados en el directorio de trabajo actual.  

**Banderas**  
`-l` `--list_species`   
Enumera todas las especies disponibles. (Para Python: combina con `species=None`.)  

`-ftp` `--ftp`   
Regresa solo los enlaces FTP solicitados.  

`-d` `--download`   
Solo para Terminal. Descarga los FTP solicitados al directorio actual (requiere [curl](https://curl.se/docs/) para ser instalado).  

`-q` `--quiet`   
Solo para la Terminal. Impide la informacion de progreso de ser exhibida durante la corrida.  
Para Python, usa `verbose=False` para imipidir la informacion de progreso de ser exhibida durante la corrida.  
  
### Por ejemplo
**Use `gget ref` en combinación con [kallisto | bustools](https://www.kallistobus.tools/kb_usage/kb_ref/) para construir un índice de referencia:**
```bash
kb ref -i INDEX -g T2G -f1 FASTA $(gget ref --ftp -w dna,gtf homo_sapiens)
```
&rarr; kb ref crea un índice de referencia utilizando los últimos archivos de ADN y GTF de especies **Homo sapiens** que le ha pasado `gget ref`.  

<br/><br/>

**Enumere todos los genomas disponibles de la versión 103 de Ensembl:**  
```bash
gget ref --list_species -r 103
```
```python
# Python
gget.ref(species=None, list_species=True, release=103)
```
&rarr; Regresa una lista con todos los genomas disponibles (`gget ref` verifica si GTF y FASTA están disponibles) de la versión 103 de Ensembl.  
(Si no se especifica ninguna versión, `gget ref` siempre devolverá información de la última versión de Ensembl).  

<br/><br/>
  
**Obtenga la referencia del genoma para una especie específica:**   
```bash
gget ref -w gtf,dna homo_sapiens
```
```python
# Python
gget.ref("homo_sapiens", which=["gtf", "dna"])
```
&rarr; Regresa un JSON con los últimos FTP humanos GTF y FASTA, y sus respectivos metadatos, en el formato:
```
{
    "homo_sapiens": {
        "annotation_gtf": {
            "ftp": "http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz",
            "ensembl_release": 106,
            "release_date": "28-Feb-2022",
            "release_time": "23:27",
            "bytes": "51379459"
        },
        "genome_dna": {
            "ftp": "http://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
            "ensembl_release": 106,
            "release_date": "21-Feb-2022",
            "release_time": "09:35",
            "bytes": "881211416"
        }
    }
}
```

#### [Más ejemplos](https://github.com/pachterlab/gget_examples)
