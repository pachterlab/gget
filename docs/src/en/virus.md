[<kbd> View page source on GitHub </kbd>](https://github.com/pachterlab/gget/blob/main/docs/src/en/virus.md)

> Python arguments are equivalent to long-option arguments (`--arg`), unless otherwise specified. Flags are True/False arguments in Python. The manual for any gget tool can be called from the command-line using the `-h` `--help` flag.  
# gget virus 🦠  

Download viral nucleotide sequences, along with rich, linked metadata, from across the International Nucleotide Sequence Database Collaboration ([INSDC](https://www.insdc.org/)), including NCBI, [ENA](ebi.ac.uk/ena/browser/), and [DDBJ](https://www.ddbj.nig.ac.jp/index-e.html) (accessed via [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/)), with the option to further enrich results using metadata from NCBI GenBank (e.g. gene and protein annotations, amino acid sequences, and more). `gget virus` applies sequential server-side and local filters to efficiently download customized datasets.  

Return format: FASTA, CSV, and JSONL files saved to an output folder.  

[No-code, shareable Google Colab notebook for querying viral sequences.](https://colab.research.google.com/github/pachterlab/gget_examples/blob/main/gget_virus/gget_virus_colab.ipynb)

This module was written by [Ferdous Nasri](https://github.com/ferbsx).

**Note**: For SARS-CoV-2 and Alphainfluenza (Influenza A) queries, `gget virus` uses NCBI's optimized cached data packages via the [NCBI datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/). The datasets CLI binary is bundled with gget for all major platforms—no additional installation required. If you already have the NCBI datasets CLI installed on your system, gget will automatically use your existing installation.

**Positional argument**  
`virus`  
Virus taxon name (e.g. 'Zika virus'), taxon ID (e.g. '2697049'), NCBI accession number (e.g. 'NC\_045512.2'), space-separated list of accessions (e.g. 'NC\_045512.2 MN908947.3 MT020781.1'), or path to a text file containing accession numbers (one per line) (e.g. 'path/to/text.txt').  

Add `--is_accession` when passing an NCBI accession number. Add `--is_sars_cov2` or `is_alphainfluenza` for optimized download of SARS-CoV2 or Alphainfluenza sequences, respectively.  

For SARS-CoV-2 and Alphainfluenza cached downloads, supports:  
  - Single accession: `NC_045512.2`  
  - Space-separated list: `NC_045512.2 MN908947.3 MT020781.1`   
  - Text file path: `accessions.txt` (one accession per line)  

Use flag `--download_all_accessions` to apply filters without searching for a specific virus.  

**Optional arguments**   

_Host filters_  

`--host`  
Filter by host organism name or NCBI Taxonomy ID (e.g., 'human', 'Aedes aegypti', `1335626`).

_Sequence & Gene filters_  

`--nuc_completeness`  
Filter by nucleotide completeness. One of: 'complete' or 'partial'.  
Set to 'complete' to only return nucleotide sequences marked as complete; set to 'partial' to only return sequences that are marked as partial.  

`--min_seq_length`  
Filter by minimum sequence length.

`--max_seq_length`  
Filter by maximum sequence length.

`--min_gene_count`  
Filter by minimum number of genes.

`--max_gene_count`  
Filter by maximum number of genes.

`--min_protein_count`  
Filter by minimum number of proteins.

`--max_protein_count`  
Filter by maximum number of proteins.

`--min_mature_peptide_count`  
Filter by minimum number of mature peptides.

`--max_mature_peptide_count`  
Filter by maximum number of mature peptides.

`--max_ambiguous_chars`  
Filter by maximum number of ambiguous nucleotide characters (N's).

`--has_proteins`  
Filter for sequences containing specific proteins or genes (e.g., 'spike', 'ORF1ab'). Can be a single protein name or a list of protein names.  
Python: `has_proteins="spike"` or `has_proteins=["spike", "ORF1ab"]`

`--annotated`  
'true' or 'false'. Filter for sequences that have been annotated with gene/protein information.  
Command line: `--annotated true` to fetch only that have been annotated with gene/protein information, or `--annotated false` to exclude them.  
Python: `annotated=True` or `annotated=False` (`annotated=None` for no filter).

`--lab_passaged`  
'true' or 'false'. Filter for or against lab-passaged samples.   
Command line: `--lab_passaged true` to fetch only lab-passaged samples, or `--lab_passaged false` to exclude them.  
Python: `lab_passaged=True` or `lab_passaged=False` (`lab_passaged=None` for no filter).

`--vaccine_strain`  
Filter for or against vaccine strain sequences.  
Command line: `--vaccine_strain true` to fetch only vaccine strains, or `--vaccine_strain false` to exclude them.  
Python: `vaccine_strain=True` or `vaccine_strain=False` (`vaccine_strain=None` for no filter).

`--segment`  
Filter for sequences with specific segment(s) (e.g. 'HA', 'NA'). Can be a single segment name or a list of segment names.
Python: `segment="HA"` or `segment=["HA", "NA", "PB1"]`

_Date filters_  

`--min_collection_date`  
Filter by minimum sample collection date (YYYY-MM-DD).

`--max_collection_date`  
Filter by maximum sample collection date (YYYY-MM-DD).

`--min_release_date`  
Filter by minimum sequence release date (YYYY-MM-DD).

`--max_release_date`  
Filter by maximum sequence release date (YYYY-MM-DD).

_Location & Submitter filters_

`--geographic_location`  
Filter by geographic location of sample collection (e.g., 'USA', 'Asia').

`--submitter_country`  
Filter by the country of the sequence submitter. Can be a single country or a comma-separated list.

`--source_database`  
Filter by source database. One of: 'genbank' or 'refseq'.

_SARS-CoV-2 specific filters_

`--lineage`  
Filter by SARS-CoV-2 lineage (e.g. 'B.1.1.7', 'P.1'). Can be a single lineage or a list of lineages.
Python: `lineage="B.1.1.7"` or `lineage=["B.1.1.7", "P.1"]`

_Workflow configurations_

`--genbank_batch_size`  
Batch size for GenBank metadata API requests. Default: 200. Larger batches are faster but may be more prone to timeouts.  

`-o` `--out`  
Path to the folder where results will be saved. Default: current working directory.  
Python: `outfolder="path/to/folder"`

**Flags**  
`-a` `--is_accession`  
Flag to indicate that the `virus` positional argument is an accession number, a space-separated list of accessions, or a path to a text file containing accession numbers (one per line).  

`--download_all_accessions`   
Use this flag when applying filters without searching for a specific virus (leave `virus` argument empty).     
⚠️ **WARNING**: If you do not specify additional filters, this flag downloads ALL available viral sequences from NCBI (entire Viruses taxonomy, taxon ID 10239). This is an extremely large dataset that can take many hours to download and require significant disk space. Use with caution and ensure you have adequate storage and bandwidth. When this flag is set, the `virus` argument is ignored.

`--is_sars_cov2`  
Use NCBI's optimized cached data packages for a SARS-CoV-2 query. This provides faster and more reliable downloads. The system can auto-detect SARS-CoV-2 taxon-name queries, but for accession-based queries, you must set this flag explicitly.

`--is_alphainfluenza`  
Use NCBI's optimized cached data packages for an Alphainfluenza (Influenza A virus) query. This provides faster and more reliable downloads for large Influenza A datasets. The system can auto-detect Alphainfluenza taxon-name queries, but for accession-based queries, you must set this flag explicitly.

`-g` `--genbank_metadata`  
Fetch and save additional detailed metadata from GenBank, including collection dates, host details, and publication references, in a separate `{virus}_genbank_metadata.csv` file (plus full XML/CSV dumps).

`--proteins_complete`  
Flag to only include sequences where all annotated proteins are complete.  

`-kt` `--keep_temp`  
Flag to keep all intermediate/temporary files generated during processing. By default, only final output files are retained.

`-q` `--quiet`  
Command-line only. Prevents progress information from being displayed.  
Python: Use `verbose=False` to prevent progress information from being displayed. 

### Example

```bash
gget virus "Zika virus" --nuc_completeness complete --host human --out zika_data
```

```python
# Python
import gget

gget.virus(
  "Zika virus",
  nuc_completeness="complete",
  host="human",
  outfolder="zika_data"
)
```

→ Downloads complete Zika virus genomes from human hosts. Results are saved in the `zika_data` folder as `Zika_virus_sequences.fasta`, `Zika_virus_metadata.csv`, `Zika_virus_metadata.jsonl`, and `command_summary.txt`.


<br><br>
**Download a specific SARS-CoV-2 reference genome using its accession number:**

```bash
gget virus NC_045512.2 --is_accession --is_sars_cov2
```

```python
# Python
import gget

gget.virus("NC_045512.2", is_accession=True, is_sars_cov2=True)
```

→ Uses the optimized download method for SARS-CoV-2 to fetch the reference genome and its metadata.

<br><br>
**Download SARS-CoV-2 sequences with cached optimization AND GenBank metadata:**

```bash
gget virus "SARS-CoV-2" --host human --nuc_completeness complete --min_seq_length 29000 --genbank_metadata
```

```python
# Python
import gget

gget.virus(
  "SARS-CoV-2", 
  host="human", 
  nuc_completeness="complete",
  min_seq_length=29000,
  genbank_metadata=True,
  is_sars_cov2=True,
  outfolder="covid_data"
)
```

→ Uses cached download for speed (via NCBI's SARS-CoV-2 data packages when available), applies the sequence length filter post-download, and fetches detailed GenBank metadata for all filtered sequences.

<br><br>
**Download Influenza A virus sequences with optimized caching and post-download filtering:**

```bash
gget virus "Influenza A virus" --host human --nuc_completeness complete --max_seq_length 15000 --genbank_metadata --is_alphainfluenza
```

```python
# Python
import gget

gget.virus(
  "Influenza A virus", 
  host="human", 
  nuc_completeness="complete",
  max_seq_length=15000,
  genbank_metadata=True,
  is_alphainfluenza=True,
  outfolder="influenza_a_data"
)
```

→ Uses NCBI's cached data packages for Alphainfluenza to download complete Influenza A genomes from human hosts much faster than the standard API method, then applies the sequence length filter and fetches GenBank metadata.

#### [More examples](https://github.com/pachterlab/gget_examples/tree/main/gget_virus)

# References

If you use `gget virus` in a publication, please cite the following articles:

  - Nasri, F. et al (2026). Coming soon.

  - Luebbert, L., & Pachter, L. (2023). Efficient querying of genomic reference databases with gget. Bioinformatics. [https://doi.org/10.1093/bioinformatics/btac836](https://doi.org/10.1093/bioinformatics/btac836)

  - O’Leary, N.A., Cox, E., Holmes, J.B. et al (2024). Exploring and retrieving sequence and metadata for species across the tree of life with NCBI Datasets. Sci Data 11, 732. [https://doi.org/10.1038/s41597-024-03571-y](https://doi.org/10.1038/s41597-024-03571-y)


