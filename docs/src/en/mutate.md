> Python arguments are equivalent to long-option arguments (`--arg`), unless otherwise specified. Flags are True/False arguments in Python.  The manual for any gget tool can be called from the command-line using the `-h` `--help` flag.  
## gget mutate 🧟
Takes in nucleotide sequences and mutations (in [standard mutation annotation](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1867422/) and returns mutated versions of the input sequences according to the provided mutations.  
Return format: Saves mutated sequences in FASTA format (or returns a list containing the mutated sequences if `out=None`).  

This module was written by [Joseph Rich](https://github.com/josephrich98).

**Positional argument**  
`sequences`   
Path to the FASTA file containing the sequences to be mutated, e.g., 'path/to/seqs.fa'.  
Sequence identifiers following the '>' character must correspond to the identifiers in the seq_ID column of `mutations`.  

Example format of the FASTA file:  
```
>seq1 (or ENSG00000106443)  
ACTGCGATAGACT  
>seq2  
AGATCGCTAG
```

Alternatively: Input sequence(s) as a string or list, e.g. 'AGCTAGCT'.

NOTE: Only the letters until the first space or dot will be used as sequence identifiers - Version numbers of Ensembl IDs will be ignored.  
NOTE: When the `sequences` input is a genome fasta file, also see the `gtf` argument below.

**Required arguments**  
`-m` `--mutations`  
Path to the csv or tsv file (e.g., 'path/to/mutations.csv') or data frame (DataFrame object) containing information about the mutations in the following format (the 'notes' and 'mut_ID' columns are optional):  

| mutation         | mut_ID | seq_ID | notes |
|------------------|--------|--------|-------|
| c.2C>T           | mut1   | seq1   | -> Apply mutation 1 to sequence 1 |
| c.9_13inv        | mut2   | seq2   | -> Apply mutation 2 to sequence 2 |
| c.9_13inv        | mut2   | seq4   | -> Apply mutation 2 to sequence 4 |
| c.9_13delinsAAT  | mut3   | seq4   | -> Apply mutation 3 to sequence 4 |
| ...              | ...    | ...    |                                   |

'mutation' = Column containing the mutations to be performed written in standard mutation annotation  
'mut_ID' = Column containing the identifier for each mutation  
'seq_ID' = Column containing the identifiers of the sequences to be mutated (must correspond to the string following the '>' character in the 'sequences' FASTA file; do NOT include spaces or dots)  

Alternatively: Input mutation(s) as a string or list, e.g., 'c.2C>T'.  
If a list is provided, the number of mutations must equal the number of input sequences.  

For use from the terminal (bash): Enclose individual mutation annotations in quotation marks to prevent parsing errors.  

**Optional input-related arguments**  
`-mc` `--mut_column`  
Name of the column containing the mutations to be performed in `mutations`. Default: 'mutation'.  

`-sic` `--seq_id_column`  
Name of the column containing the IDs of the sequences to be mutated in `mutations`. Default: 'seq_ID'.

`-mic` `--mut_id_column`  
Name of the column containing the IDs of each mutation in `mutations`. Default: Same as `mut_column`.

`-gtf` `--gtf`  
Path to a .gtf file. When providing a genome fasta file as input for 'sequences', you can provide a .gtf file here and the input sequences will be defined according to the transcript boundaries, e.g. 'path/to/genome_annotation.gtf'. Default: None

`-gtic` `--gtf_transcript_id_column`  
Column name in the input `mutations` file containing the transcript ID. In this case, column `seq_id_column` should contain the chromosome number.  
Required when `gtf` is provided. Default: None  
  
**Optional mutant sequence generation/filtering arguments**  
`-k` `--k`  
Length of sequences flanking the mutation. Default: 30.  
If k > total length of the sequence, the entire sequence will be kept.  

`-msl` `--min_seq_len`  
Minimum length of the mutant output sequence, e.g. 100. Mutant sequences smaller than this will be dropped. Default: None

`-ma` `--max_ambiguous`                
Maximum number of 'N' (or 'n') characters allowed in the output sequence, e.g. 10. Default: None (no ambiguous character filter will be applied)

**Optional mutant sequence generation/filtering flags**  
`-ofr` `--optimize_flanking_regions`  
Removes nucleotides from either end of the mutant sequence to ensure (when possible) that the mutant sequence does not contain any k-mers also found in the wildtype/input sequence. 

`-rswk` `--remove_seqs_with_wt_kmers`  
Removes output sequences where at least one k-mer is also present in the wildtype/input sequence in the same region.  
When used with `--optimize_flanking_regions`, only sequences for which a wildtype k-mer is still present after optimization will be removed.

`-mio` `--merge_identical_off`          
Do not merge identical mutant sequences in the output (by default, identical sequences will be merged by concatenating the sequence headers for all identical sequences).

**Optional arguments to generate additional output**   
This output is activated using the  `--update_df` flag and will be stored in a copy of the `mutations` DataFrame.  

`-udf_o` `--update_df_out`               
Path to output csv file containing the updated DataFrame, e.g. 'path/to/mutations_updated.csv'. Only valid when used with `--update_df`.  
Default: None -> the new csv file will be saved in the same directory as the `mutations` DataFrame with appendix '_updated'  

`-ts` `--translate_start`              
(int or str) The position in the input nucleotide sequence to start translating, e.g. 5. If a string is provided, it should correspond to a column name in `mutations` containing the open reading frame start positions for each sequence/mutation. Only valid when used with `--translate`.  
Default: translates from the beginning of each sequence  

`-te` `--translate_end`                
(int or str) The position in the input nucleotide sequence to end translating, e.g. 35. If a string is provided, it should correspond to a column name in `mutations` containing the open reading frame end positions for each sequence/mutation. Only valid when used with `--translate`.  
Default: translates until the end of each sequence  

**Optional flags to modify additional output**  
`-udf` `--update_df`   
Updates the input `mutations` DataFrame to include additional columns with the mutation type, wildtype nucleotide sequence, and mutant nucleotide sequence (only valid if `mutations` is a .csv or .tsv file).  

`-sfs` `--store_full_sequences`         
Includes the complete wildtype and mutant sequences in the updated `mutations` DataFrame (not just the sub-sequence with k-length flanks). Only valid when used with `--update_df`.   

`-tr` `--translate`                  
Adds additional columns to the updated `mutations` DataFrame containing the wildtype and mutant amino acid sequences. Only valid when used with `--store_full_sequences`.   
                                  
**Optional general arguments**  
`-o` `--out`   
Path to output FASTA file containing the mutated sequences, e.g., 'path/to/output_fasta.fa'.  
Default: None -> returns a list of the mutated sequences to standard out.    
The identifiers (following the '>') of the mutated sequences in the output FASTA will be '>[seq_ID]_[mut_ID]'. 

**Optional general flags**  
`-q` `--quiet`   
Command-line only. Prevents progress information from being displayed.  
Python: Use `verbose=False` to prevent progress information from being displayed.  

### Examples
```bash
gget mutate ATCGCTAAGCT -m 'c.4G>T'
```
```python
# Python
gget.mutate("ATCGCTAAGCT", "c.4G>T")
```
&rarr; Returns ATCTCTAAGCT.  

<br/><br/>

**List of sequences with a mutation for each sequence provided in a list:**  
```bash
gget mutate ATCGCTAAGCT TAGCTA -m 'c.4G>T' 'c.1_3inv' -o mut_fasta.fa
```
```python
# Python
gget.mutate(["ATCGCTAAGCT", "TAGCTA"], ["c.4G>T", "c.1_3inv"], out="mut_fasta.fa")
```
&rarr; Saves 'mut_fasta.fa' file containing: 
```
>seq1_mut1  
ATCTCTAAGCT  
>seq2_mut2  
GATCTA
```

<br/><br/>

**One mutation applied to several sequences with adjusted `k`:**  
```bash
gget mutate ATCGCTAAGCT TAGCTA -m 'c.1_3inv' -k 3
```
```python
# Python
gget.mutate(["ATCGCTAAGCT", "TAGCTA"], "c.1_3inv", k=3)
```
&rarr; Returns ['CTAGCT', 'GATCTA'].  


<br/><br/>

**Add mutations to an entire genome with extended output**  
Main input:   
- mutation information as a `mutations` CSV (by having `seq_id_column` contain chromosome information, and `mut_column` contain mutation information with respect to genome coordinates)  
- the genome as the `sequences` file  

Since we are passing the path to a gtf file to the `gtf` argument, transcript boundaries will be respected (the genome will be split into transcripts). `gtf_transcript_id_column` specifies the name of the column in `mutations` containing the transcript IDs corresponding to the transcript IDs in the `gtf` file.  

The `optimize_flanking_regions` argument maximizes the length of the resulting mutation-containing sequences while maintaining specificity (no wildtype k-mer will be retained).

`update_df` activates the creation of a new CSV file with updated information about each input and output sequence. This new CSV file will be saved as `update_df_out`. Since `store_full_sequences` is activated, this new CSV file will not only contain the output sequences (restricted in size by flanking regiong of size `k`), but also the complete input and output sequences. This allows us to observe the mutation in the context of the entire sequence. Lastly, we are also adding the translated versions of the complete sequences by adding the with the `translate` flag, so we can observe how the resulting amino acid sequence is changed. The `translate_start` and `translate_end` arguments specify the names of the columns in `mutations` that contain the start and end positions of the open reading frame (start and end positions for translating the nucleotide sequence to an amino acid sequence), respectively.  

```bash
gget mutate \
  -m mutations_input.csv \
  -o mut_fasta.fa \
  -k 4 \
  -sic Chromosome \
  -mic Mutation \
  -gtf genome_annotation.gtf \
  -gtic Ensembl_Transcript_ID \
  -ofr \
  -update_df \
  -udf_o mutations_updated.csv \
  -sfs \
  -tr \
  -ts Translate_Start \
  -te Translate_End \
  genome_reference.fa
```
```python
# Python
gget.mutate(sequences="genome_reference.fa", mutations="mutations_input.csv", out="mut_fasta.fa", k=4, seq_id_column="Chromosome", mut_column="Mutation", gtf="genome_annotation.gtf", gtf_transcript_id_column="Ensembl_Transcript_ID", optimize_flanking_regions=True, update_df=True, update_df_out="mutations_updated.csv", store_full_sequences=True, translate=True, translate_start="Translate_Start", translate_end="Translate_End")
```
&rarr; Takes as input 'mutations_input.csv' file containing: 
```
| Chromosome | Mutation          | Ensembl_Transcript_ID  | Translate_Start | Translate_End |
|------------|-------------------|------------------------|-----------------|---------------|
| 1          | g.224411A>C       | ENST00000193812        | 0               | 100           |
| 8          | g.25111del        | ENST00000174411        | 0               | 294           |
| X          | g.1011_1012insAA  | ENST00000421914        | 9               | 1211          |
``` 
&rarr; Saves 'mut_fasta.fa' file containing: 
```
>1:g.224411A>C  
TGCTCTGCT  
>8:g.25111del  
GAGTCGAT
>X:g.1011_1012insAA
TTAGAACTT
``` 
&rarr; Saves 'mutations_updated.csv' file containing: 
```

| Chromosome | Mutation          | Ensembl_Transcript_ID  | mutation_type | wt_sequence | mutant_sequence | wt_sequence_full  | mutant_sequence_full | wt_sequence_aa_full | mutant_sequence_aa_full |
|------------|-------------------|------------------------|---------------|-------------|-----------------|-------------------|----------------------|---------------------|-------------------------|
| 1          | g.224411A>C       | ENSMUST00000193812     | Substitution  | TGCTATGCT   | TGCTCTGCT       | ...TGCTATGCT...   | ...TGCTCTGCT...      | ...CYA...           | ...CSA...               |
| 8          | g.25111del        | ENST00000174411        | Deletion      | GAGTCCGAT   | GAGTCGAT        | ...GAGTCCGAT...   | ...GAGTCGAT...       | ...ESD...           | ...ES...                |
| X          | g.1011_1012insAA  | ENST00000421914        | Insertion     | TTAGCTT     | TTAGAACTT       | ...TTAGCTT...     | ...TTAGAACTT...      | ...A...             | ...EL...                |

```

# References
If you use `gget mutate` in a publication, please cite the following articles:   

- Luebbert, L., & Pachter, L. (2023). Efficient querying of genomic reference databases with gget. Bioinformatics. [https://doi.org/10.1093/bioinformatics/btac836](https://doi.org/10.1093/bioinformatics/btac836)
