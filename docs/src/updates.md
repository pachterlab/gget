## ✨ What's new  
**Version ≥ 0.27.3:**  
- [`gget info`](./info.md) excludes PDB IDs by default to increase speed (PDB results can be included using flag `--pdb` / `pdb=True`).

**Version ≥ 0.27.2:**  
- Updated [`gget alphafold`](./alphafold.md) to [DeepMind's AlphaFold v2.3.0](https://github.com/deepmind/alphafold/releases/tag/v2.3.0) (including new arguments `multimer_for_monomer` and `multimer_recycles`)

**Version ≥ 0.27.0:**  
- Updated [`gget alphafold`](./alphafold.md) to match recent changes by DeepMind  
- Updated version number to match gget's creator's age following a long-standing Pachter lab tradition  

**Version ≥ 0.3.13:**  
- Reduced runtime for [`gget enrichr`](./enrichr.md) and [`gget archs4`](./archs4.md) when used with Ensembl IDs  

**Version ≥ 0.3.12:**  
- [`gget info`](./info.md) now also returns subcellular localisation data from UniProt
- New [`gget info`](./info.md) flag `ensembl_only` returns only Ensembl results
- Reduced runtime for [`gget info`](./info.md) and [`gget seq`](./seq.md)

**Version ≥ 0.3.11: [`gget pdb`](./pdb.md)**  

**NOTE:** The [Ensembl FTP site](http://ftp.ensembl.org/pub/) changed its structure on August 8, 2022. Please upgrade to `gget` version ≥ 0.3.7 if you use `gget ref`.  

**Version ≥ 0.3.0: [`gget alphafold`](./alphafold.md)**  

**NOTE:** [UniProt](https://www.uniprot.org/) changed the structure of their API on June 28, 2022. Please upgrade to `gget` version ≥ 0.2.5 if you use any of the modules querying data from UniProt (`gget info` and `gget seq`).

**Version ≥ 0.2.0:**  
- JSON is now the default output format for the command-line interface for modules that previously returned data frame (CSV) format by default (the output can be converted to data frame/CSV using flag `[-csv][--csv]`). Data frame/CSV remains the default output for Jupyter Lab / Google Colab (and can be converted to JSON with `json=True`).
- For all modules, the first required argument was converted to a positional argument and should not be named anymore in the command-line, e.g. `gget ref -s human` &rarr; `gget ref human`.
- `gget info`: `[--expand]` is deprecated. The module will now always return all of the available information.
- Slight changes to the output returned by `gget info`, including the return of versioned Ensembl IDs.
- `gget info` and `gget seq` now support 🪱 WormBase and 🪰 FlyBase IDs.
- `gget archs4` and `gget enrichr` now also take Ensembl IDs as input with added flag `[-e][--ensembl]` (`ensembl=True` in Jupyter Lab / Google Colab).
- `gget seq` argument `seqtype` was replaced by flag `[-t][--translate]` (`translate=True/False` in Jupyter Lab / Google Colab) which will return either nucleotide (`False`) or amino acid (`True`) sequences.
- `gget search` argument `seqtype` was renamed to `id_type` for clarity (still taking the same arguments 'gene' or 'transcript').
- Version ≥ 0.2.6: `gget ref` supports plant genomes! 🌱  
