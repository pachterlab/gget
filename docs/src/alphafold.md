> Python arguments are equivalent to long-option arguments (`--arg`), unless otherwise specified. Flags are True/False arguments in Python. The manual for any gget tool can be called from the command-line using the `-h` `--help` flag.  
## gget alphafold 🪢
Predict the 3D structure of a protein from its amino acid sequence using a simplified version of [DeepMind](https://www.deepmind.com/)’s [AlphaFold2](https://github.com/deepmind/alphafold) originally released and benchmarked for [AlphaFold Colab](https://colab.research.google.com/github/deepmind/alphafold/blob/main/notebooks/AlphaFold.ipynb).  
Returns: Predicted structure (PDB) and alignment error (json).  

Before using `gget alphafold` for the first time, run `gget setup alphafold` / `gget.setup("alphafold")` once (also see [`gget setup`](setup.md)).  

**Positional argument**  
`sequence`  
Amino acid sequence (str), list of sequences (for multimers), or path to FASTA file.

**Optional arguments**  
`-mr` `--multimer_recycles`  
The multimer model will continue recycling until the predictions stop changing, up to the limit set here. Default: 3.  
For higher accuracy, at the potential cost of longer inference times, set this to 20.  

`-o` `--out`   
Path to folder to save prediction results in (str). Default: "./[date_time]_gget_alphafold_prediction".  
  
**Flags**   
`-mfm` `--multimer_for_monomer`  
Use multimer model for a monomer.  

`-r` `--relax`   
AMBER relax the best model. 

`plot`  
Python only. `plot=True` provides an interactive, 3D graphical overview of the predicted structure and alignment quality using [py3Dmol](https://pypi.org/project/py3Dmol/) and [matplotlib](https://matplotlib.org/) (default: True).  

`show_sidechains`  
Python only. `show_sidechains=True` includes side chains in the plot (default: True).  

`-q` `--quiet`   
Command-line only. Prevents progress information from being displayed.  
Python: Use `verbose=False` to prevent progress information from being displayed. 
  
  
### Example
```bash
# Generate new prediction from amino acid sequence
gget alphafold MAAHKGAEHHHKAAEHHEQAAKHHHAAAEHHEKGEHEQAAHHADTAYAHHKHAEEHAAQAAKHDAEHHAPKPH

# Find similar sequences deposited on the PDB for comparative analysis
gget blast --database pdbaa MAAHKGAEHHHKAAEHHEQAAKHHHAAAEHHEKGEHEQAAHHADTAYAHHKHAEEHAAQAAKHDAEHHAPKPH

# Fetch the PDB files of similar structures returned by gget blast for comparison, to get a measure for model quality
gget pdb 3UQ3 -o 3UQ3.pdb
gget pdb 2K42 -o 2K42.pdb
```
```python
# Python
gget.alphafold("MAAHKGAEHHHKAAEHHEQAAKHHHAAAEHHEKGEHEQAAHHADTAYAHHKHAEEHAAQAAKHDAEHHAPKPH")

# Find similar sequences deposited on the PDB for comparative analysis
gget.blast("MAAHKGAEHHHKAAEHHEQAAKHHHAAAEHHEKGEHEQAAHHADTAYAHHKHAEEHAAQAAKHDAEHHAPKPH", database="pdbaa")

# Fetch the PDB files of similar structures returned by gget blast for comparison, to get a measure for model quality
gget.pdb("3UQ3", save=True)
gget.pdb("2K42", save=True)
```
&rarr; `gget alphafold` returns the predicted structure (PDB) and predicted alignment error (.json) in a new folder ("./[date_time]_gget_alphafold_prediction"). The use case above exemplifies how to use [`gget blast`](blast.md) and [`gget pdb`](pdb.md) for a comparative analysis of the new prediction. PDB files can be viewed interactively in 3D [online](https://rcsb.org/3d-view), or using programs like [PyMOL](https://pymol.org/) or [Blender](https://www.blender.org/). To compare two PDB files, you can use [this website](https://rcsb.org/alignment). The Python interface also returns [interactive plots](https://twitter.com/NeuroLuebbert/status/1555968042948915200), which can be generated from the PDB and JSON as described in the [gget alphafold FAQ](https://github.com/pachterlab/gget/discussions/39) Q4.

#### [Example in Google Colab](https://github.com/pachterlab/gget_examples/blob/main/gget_alphafold.ipynb)
#### Also see: [gget alphafold FAQ](https://github.com/pachterlab/gget/discussions/39)
