> Parámetros de Python són iguales a los parámetros largos (`--parámetro`) de Terminal, si no especificado de otra manera. Banderas son parámetros de verdadero o falso (True/False) en Python. El manuál para cualquier modulo de gget se puede llamar desde el Terminal con la bandera `-h` `--help`.  
## gget alphafold 🪢
Predice la structura en 3D de cualquier proteína basada sobre su secuencia de aminoácidos usando una versión simpleficada del algoritmo [AlphaFold2](https://github.com/deepmind/alphafold) de [DeepMind](https://www.deepmind.com/), originalmente producido i hecho público para [AlphaFold Colab](https://colab.research.google.com/github/deepmind/alphafold/blob/main/notebooks/AlphaFold.ipynb).  
Regresa: La structura pedicada (en formato PDB) i el errór de alineación (en formato json).  

Antes de usar `gget alphafold` por primera vez, corre `gget setup alphafold` / `gget.setup("alphafold")` (ver también [`gget setup`](setup.md)).  

**Parámetro posicional**  
`sequence`  
Secuencia de aminoácidos (str), o una lista de secuencias (*gget alphafold automaticamente usa el algoritmo de multímero si múltiple secuencias son ingresadas*), o una ruta a un archivo tipo FASTA.  

**Parámetros optionales**  
`-mr` `--multimer_recycles`  
El algoritmo de multímero continuara a recyclar hasta que los prediciones paren de cambiar, hasta el limite indicado aquí. Por defecto: 3  
Para obtener más exactitude ajusta a 20 (al costo de corridas más largas).  

`-o` `--out`   
Ruta a una carpeta para guardar los resultados de la predicción (str). Por defecto: "./[fecha_tiempo]_gget_alphafold_prediction".  
   
**Banderas**   
`-mfm` `--multimer_for_monomer`  
Usa el algoritmo de multímero para un monómero.  

`-r` `--relax`   
Relaja el mejor modelo con el algoritmo AMBER.  

`-q` `--quiet`   
Solo para la Terminal. Impide la informacion de progreso de ser exhibida durante la corrida.  
Para Python, usa `verbose=False` para imipidir la informacion de progreso de ser exhibida durante la corrida.  

`plot`  
Solo para Python. `plot=True` provee una visualición interactiva de la predicción con el errór de alineación en 3D con [py3Dmol](https://pypi.org/project/py3Dmol/) i [matplotlib](https://matplotlib.org/) (por defecto: True).  

`show_sidechains`  
Solo para Python. `show_sidechains=True` incluye las cadenas laterales de proteínas en la visualición (por defecto: True).  
  
  
### Por ejemplo
```bash
# Predice la structura de una proteína basada sobre su secuencia de aminoácidos
gget alphafold MAAHKGAEHHHKAAEHHEQAAKHHHAAAEHHEKGEHEQAAHHADTAYAHHKHAEEHAAQAAKHDAEHHAPKPH

# Encuentra secuencias similares depositadas en el PDB para análisis comparativo
gget blast --database pdbaa MAAHKGAEHHHKAAEHHEQAAKHHHAAAEHHEKGEHEQAAHHADTAYAHHKHAEEHAAQAAKHDAEHHAPKPH

# Busca los archivos del PDB por las proteínas encuentradas con gget blast para tener algo con que comparar la predicción
gget pdb 3UQ3 -o 3UQ3.pdb
gget pdb 2K42 -o 2K42.pdb
```
```python
# Python
# Predice la structura de una proteína basada sobre su secuencia de aminoácidos
gget.alphafold("MAAHKGAEHHHKAAEHHEQAAKHHHAAAEHHEKGEHEQAAHHADTAYAHHKHAEEHAAQAAKHDAEHHAPKPH")

# Encuentra secuencias similares depositadas en el PDB para análisis comparativo
gget.blast("MAAHKGAEHHHKAAEHHEQAAKHHHAAAEHHEKGEHEQAAHHADTAYAHHKHAEEHAAQAAKHDAEHHAPKPH", database="pdbaa")

# Busca los archivos del PDB por las proteínas encuentradas con gget blast para tener algo con que comparar la predicción
gget.pdb("3UQ3", save=True)
gget.pdb("2K42", save=True)
```
&rarr; `gget alphafold` devuelve la structura pedicada (en formato PDB) i el errór de alineación (en formato json) en una carpeta nueva ("./[fecha_tiempo]_gget_alphafold_prediction"). Este ejemplo demuestra como usar [`gget blast`](blast.md) i [`gget pdb`](pdb.md) para correr un análisis comparativo. Los archivos PDB se pueden ver en 3D con [RCSB 3D view](https://rcsb.org/3d-view), o usando programas como [PyMOL](https://pymol.org/) o [Blender](https://www.blender.org/). Para comparar múltiple archivos PDB, usen [RCSB alignment](https://rcsb.org/alignment). Python también devuelve [visualiciónes interactivas](https://twitter.com/NeuroLuebbert/status/1555968042948915200), que también se pueden generar con los archivos PDB i JSON como describido en [gget alphafold FAQ](https://github.com/pachterlab/gget/discussions/39) Q4.

### [Ejemplo en Google Colab](https://github.com/pachterlab/gget_examples/blob/main/gget_alphafold.ipynb)
### [gget alphafold FAQ](https://github.com/pachterlab/gget/discussions/39)
