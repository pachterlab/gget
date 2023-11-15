## ✨ What's new  
**Version ≥ 0.28.2** (November 15, 2023):  
- [`gget info`](./info.md): Return a logging error message when the NCBI server fails for a reason other than a fetch fail (this is an error on the server side rather than an error with `gget`)
- Replace deprecated 'text' argument to find()-type methods whenever used with dependency `BeautifulSoup`
- [`gget elm`](elm.md): Remove false positive and true negative instances from returned results
  
**Version ≥ 0.28.0** (November 5, 2023):  
- Updated documentation of [`gget muscle`](./muscle.md) to add a tutorial on how to visualize sequences with sequence name lengths + slight change to returned visualization so it's a bit more robust to varying sequence names
- [`gget muscle`](./muscle.md) now also allows a list of sequences as input (as an alternative to providing the path to a FASTA file)
- Allow missing gene filter for [`gget cellxgene`](cellxgene.md)  (fixes [bug](https://github.com/pachterlab/gget/issues/110))
- [`gget seq`](./seq.md): Allow missing gene names (fixes [https://github.com/pachterlab/gget/issues/107](https://github.com/pachterlab/gget/issues/107))
- New arguments for [`gget enrichr`](enrichr.md): Use arguments `kegg_out` and `kegg_rank` to create an image of the KEGG pathway with the genes from the enrichment analysis highlighted (thanks to [this PR](https://github.com/pachterlab/gget/pull/106) by [Noriaki Sato](https://github.com/noriakis))  
- New modules: [`gget elm`](elm.md) and [`gget diamond`](diamond.md)  
  
**Version ≥ 0.27.9** (August 7, 2023):  
- New arguments for [`gget enrichr`](enrichr.md): Use argument `background_list` to provide a list of background genes
- [`gget search`](search.md) now also searches [Ensembl](https://ensembl.org/) synonyms (in addition to gene descriptions and names) to return more comprehensive search results (thanks to [Samuel Klein](https://github.com/KleinSamuel) for the [suggestion](https://github.com/pachterlab/gget/issues/90))

**Version ≥ 0.27.8** (July 12, 2023):  
- New argument for [`gget search`](search.md): Specify the Ensembl release from which information is fetched with `-r` `--release`
- Fixed [bug](https://github.com/pachterlab/gget/issues/91) in [`gget pdb`](pdb.md) (this bug was introduced in version 0.27.5)

**Version ≥ 0.27.7** (May 15, 2023):  
- Moved dependencies for modules [`gget gpt`](gpt.md) and [`gget cellxgene`](cellxgene.md) from automatically installed requirements to [`gget setup`](setup.md).  
- Updated [`gget alphafold`](alphafold.md) dependencies for compatibility with Python >= 3.10.  
- Added `census_version` argument to [`gget cellxgene`](cellxgene.md).

**Version ≥ 0.27.6** (May 1, 2023) (YANKED due to problems with dependencies -> replaced with version 0.27.7):  
- Thanks to PR by [Tomás Di Domenico](https://github.com/tdido): [`gget search`](search.md) can now also query plant 🌱 Ensembl IDs.  
- New module: [`gget cellxgene`](cellxgene.md)  

**Version ≥ 0.27.5** (April 6, 2023):  
- Updated [`gget search`](search.md) to function correctly with new [Pandas](https://pypi.org/project/pandas/2.0.0/) version 2.0.0 (released on April 3rd, 2023) as well as older versions of Pandas
- Updated [`gget info`](info.md) with new flags `uniprot` and `ncbi` which allow turning off results from these databases independently to save runtime (note: flag `ensembl_only` was deprecated)
- All gget modules now feature a `-q / --quiet` (Python: `verbose=False`) flag to turn off progress information

**Version ≥ 0.27.4** (March 19, 2023):  
- New module: [`gget gpt`](gpt.md)  

**Version ≥ 0.27.3** (March 11, 2023):  
- [`gget info`](info.md) excludes PDB IDs by default to increase speed (PDB results can be included using flag `--pdb` / `pdb=True`).  

**Version ≥ 0.27.2** (January 1, 2023):    
- Updated [`gget alphafold`](alphafold.md) to [DeepMind's AlphaFold v2.3.0](https://github.com/deepmind/alphafold/releases/tag/v2.3.0) (including new arguments `multimer_for_monomer` and `multimer_recycles`)  

**Version ≥ 0.27.0** (December 10, 2022):  
- Updated [`gget alphafold`](alphafold.md) to match recent changes by DeepMind  
- Updated version number to match [gget's creator](https://github.com/lauraluebbert)'s age following a long-standing Pachter lab tradition  

**Version ≥ 0.3.13** (November 11, 2022):  
- Reduced runtime for [`gget enrichr`](enrichr.md) and [`gget archs4`](archs4.md) when used with Ensembl IDs  

**Version ≥ 0.3.12** (November 10, 2022):  
- [`gget info`](info.md) now also returns subcellular localisation data from UniProt
- New [`gget info`](info.md) flag `ensembl_only` returns only Ensembl results
- Reduced runtime for [`gget info`](info.md) and [`gget seq`](seq.md)

**Version ≥ 0.3.11** (September 7, 2022):  
- New module: [`gget pdb`](pdb.md)  

**Version ≥ 0.3.10** (September 2, 2022):  
- [`gget alphafold`](alphafold.md) now also returns pLDDT values for generating plots from output without rerunning the program (also see the [gget alphafold FAQ](https://github.com/pachterlab/gget/discussions/39))

**Version ≥ 0.3.9** (August 25, 2022):  
- Updated openmm installation instructions for [`gget alphafold`](alphafold.md)  

**Version ≥ 0.3.8** (August 12, 2022):  
- Fixed mysql-connector-python version requirements

**Version ≥ 0.3.7** (August 9, 2022):  
- **NOTE:** The [Ensembl FTP site](http://ftp.ensembl.org/pub/) changed its structure on August 8, 2022. Please upgrade to `gget` version ≥ 0.3.7 if you use [`gget ref`](ref.md)  

**Version ≥ 0.3.5** (August 6, 2022):  
- New module: [`gget alphafold`](alphafold.md)  

**Version ≥ 0.2.6** (July 7, 2022):  
- [`gget ref`](ref.md) now supports plant genomes! 🌱  

**Version ≥ 0.2.5** (June 30, 2022):  
- **NOTE:** [UniProt](https://www.uniprot.org/) changed the structure of their API on June 28, 2022. Please upgrade to `gget` version ≥ 0.2.5 if you use any of the modules querying data from UniProt ([`gget info`](info.md) and [`gget seq`](seq.md)).

**Version ≥ 0.2.3:** (June 26, 2022):  
- JSON is now the default output format for the command-line interface for modules that previously returned data frame (CSV) format by default (the output can be converted to data frame/CSV using flag `[-csv][--csv]`). Data frame/CSV remains the default output for Jupyter Lab / Google Colab (and can be converted to JSON with `json=True`).
- For all modules, the first required argument was converted to a positional argument and should not be named anymore in the command-line, e.g. `gget ref -s human` &rarr; `gget ref human`.
- `gget info`: `[--expand]` is deprecated. The module will now always return all of the available information.
- Slight changes to the output returned by `gget info`, including the return of versioned Ensembl IDs.
- `gget info` and `gget seq` now support 🪱 WormBase and 🪰 FlyBase IDs.
- `gget archs4` and `gget enrichr` now also take Ensembl IDs as input with added flag `[-e][--ensembl]` (`ensembl=True` in Jupyter Lab / Google Colab).
- `gget seq` argument `seqtype` was replaced by flag `[-t][--translate]` (`translate=True/False` in Jupyter Lab / Google Colab) which will return either nucleotide (`False`) or amino acid (`True`) sequences.
- `gget search` argument `seqtype` was renamed to `id_type` for clarity (still taking the same arguments 'gene' or 'transcript').
