> Parámetros de Python són iguales a los parámetros largos (`--parámetro`) de Terminal, si no especificado de otra manera. Banderas son parámetros de verdadero o falso (True/False) en Python. El manuál para cualquier modulo de gget se puede llamar desde la Terminal con la bandera `-h` `--help`.  
## gget alphafold 🪢
Predice la estructura en 3D de cualquier proteína derivada de su secuencia de aminoácidos usando una versión simplificada del algoritmo [AlphaFold2](https://github.com/deepmind/alphafold) de [DeepMind](https://www.deepmind.com/), originalmente producido y publicado para [AlphaFold Colab](https://colab.research.google.com/github/deepmind/alphafold/blob/main/notebooks/AlphaFold.ipynb).  
Resultado: Predicción de la estructura (en formato PDB) y el errór de alineación (en formato json).  

Antes de usar `gget alphafold` por primera vez:
1. Instale openmm v7.5.1 (o v7.7.0 para Python >= 3.10) ejecutando el siguiente comando desde la línea de comando:  
   `conda install -qy conda==4.13.0 && conda install -qy -c conda-forge openmm=7.5.1`  
   (reemplazar con `openmm=7.7.0` para Python >= 3.10)  
   Recomendación: siga con `conda update -qy conda` para actualizar _conda_ a la última versión.  
3. Corre `gget setup alphafold` / `gget.setup("alphafold")` (ver también [`gget setup`](setup.md)). Al ejecutar `gget setup alphafold` / `gget.setup("alphafold")` se descargará e instalará la última versión de AlphaFold2 alojada en el [AlphaFold GitHub Repo](https://github.com/deepmind/alphafold). Puede volver a ejecutar este comando en cualquier momento para actualizar el software cuando hay una nueva versión de AlphaFold.    

**Parámetro posicional**  
`sequence`  
Secuencia de aminoácidos (str), o una lista de secuencias (*gget alphafold automaticamente usa el algoritmo del multímero si múltiples secuencias son ingresadas*), o una ruta a un archivo formato FASTA.  

**Parámetros optionales**  
`-mr` `--multimer_recycles`  
El algoritmo de multímero se reciclara hasta que las predicciones dejen de cambiar, el limite de ciclos esta indicado aqui. Por defecto: 3  
Para obtener más exactitud, ajusta este limite a 20 (al costo de ejecuciones mas tardadas).  

`-o` `--out`   
Ruta a la carpeta para guardar los resultados de la predicción (str). Por defecto: "./[fecha_tiempo]_gget_alphafold_prediction".  
   
**Banderas**   
`-mfm` `--multimer_for_monomer`  
Usa el algoritmo de multímero para un monómero.  

`-r` `--relax`   
Relaja el mejor modelo con el algoritmo AMBER.  

`-q` `--quiet`   
Uso limitado para Terminal. Impide la información de progreso de ser exhibida durante la ejecución del programa.  
Para Python, usa `verbose=False`.  

`plot`  
Solo para Python. `plot=True` provée una visualización interactiva de la predicción con el errór de alineación en 3D con [py3Dmol](https://pypi.org/project/py3Dmol/) y [matplotlib](https://matplotlib.org/) (por defecto: True).  

`show_sidechains`  
Solo para Python. `show_sidechains=True` incluye las cadenas laterales de proteínas en el esquema (por defecto: True).  
  
  
### Ejemplo
```bash
# Predice la estructura de una proteína derivada de su secuencia de aminoácidos
gget alphafold MAAHKGAEHHHKAAEHHEQAAKHHHAAAEHHEKGEHEQAAHHADTAYAHHKHAEEHAAQAAKHDAEHHAPKPH

# Encuentra secuencias similares previamente depositadas en el PDB para análisis comparativo
gget blast --database pdbaa MAAHKGAEHHHKAAEHHEQAAKHHHAAAEHHEKGEHEQAAHHADTAYAHHKHAEEHAAQAAKHDAEHHAPKPH

# Busca los archivos PDB de estructuras similares resultantes de gget blast para comparar y obtener una medida de calidad del modelo predecido.
gget pdb 3UQ3 -o 3UQ3.pdb
gget pdb 2K42 -o 2K42.pdb
```
```python
# Python
# Predice la estructura de una proteína derivada de su secuencia de aminoácidos
gget.alphafold("MAAHKGAEHHHKAAEHHEQAAKHHHAAAEHHEKGEHEQAAHHADTAYAHHKHAEEHAAQAAKHDAEHHAPKPH")

# Encuentra secuencias similares previamente depositadas en el PDB para análisis comparativo
gget.blast("MAAHKGAEHHHKAAEHHEQAAKHHHAAAEHHEKGEHEQAAHHADTAYAHHKHAEEHAAQAAKHDAEHHAPKPH", database="pdbaa")

# Busca los archivos PDB de estructuras similares resultantes de gget blast para comparar y obtener una medida de calidad del modelo predecido.
gget.pdb("3UQ3", save=True)
gget.pdb("2K42", save=True)
```
&rarr; `gget alphafold` produce la estructura predecida (en formato PDB) y el errór de alineación (en formato json) en una nueva carpeta ("./[fecha_tiempo]_gget_alphafold_prediction"). Este ejemplo demuestra como usar [`gget blast`](blast.md) y [`gget pdb`](pdb.md) para correr un análisis comparativo. Los archivos PDB se pueden ver en 3D con [RCSB 3D view](https://rcsb.org/3d-view), o usando programas como [PyMOL](https://pymol.org/) o [Blender](https://www.blender.org/). Para comparar múltiples archivos PDB, use [RCSB alignment](https://rcsb.org/alignment). Python también produce [esquemas interactivos](https://twitter.com/NeuroLuebbert/status/1555968042948915200), los cuales se pueden generar de los archivos PDB y JSON, como es describido en [gget alphafold FAQ](https://github.com/pachterlab/gget/discussions/39) Q4.

<iframe width="560" height="315" src="https://www.youtube.com/embed/4qxGF1tbZ3I?si=mEqQ5oSnDYtg2OP7" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen></iframe>

### [Ejemplo en Google Colab](https://github.com/pachterlab/gget_examplblob/main/gget_alphafold.ipynb)
### [gget alphafold - preguntas más frecuentes](https://github.com/pachterlab/gget/discussions/39)
