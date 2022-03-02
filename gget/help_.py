def help_():
    print("""
## gget ref
Fetch links to GTF and FASTA files from the Ensembl FTP site.

-l --list
List all available species.

-s --species
Species for which the FTPs will be fetched in the format genus_species, e.g. homo_sapiens.

-w --which
Defines which results to return. Possible entries are: 'all' - Returns GTF, cDNA, and DNA links and associated info (default). Or one or a combination of the following: 'gtf' - Returns the GTF FTP link and associated info. 'cdna' - Returns the cDNA FTP link and associated info. 'dna' - Returns the DNA FTP link and associated info.

-r --release
Ensemble release the FTPs will be fetched from, e.g. 104 (default: None → uses latest Ensembl release).

-ftp --ftp
If True: returns only a list containing the requested FTP links (default: False).

-d --download
Download the requested FTPs to the current directory.

-o --out
Path to the file the results will be saved in, e.g. path/to/directory/results.json (default: None → just prints results).
For Jupyter Lab / Google Colab: save=True will save the output to the current working directory.


## gget search
Query Ensembl for genes using free form search words.

-sw --searchwords
One or more free form searchwords for the query, e.g. gaba, nmda. Searchwords are not case-sensitive.

-s --species
Species or database to be searched.
Species can be passed in the format 'genus_species', e.g. 'homo_sapiens'. To pass a specific CORE database (e.g. a specific mouse strain), enter the name of the CORE database, e.g. 'mus_musculus_dba2j_core_105_1'. All availabale species databases can be found here: http://ftp.ensembl.org/pub/release-105/mysql/

-t --d_type
Possible entries: 'gene' (default), 'transcript' Returns either genes or transcripts, respectively, which match the searchwords.

-ao --andor
Possible entries: 'or', 'and' 'or': ID descriptions must include at least one of the searchwords (default). 'and': Only return IDs whose descriptions include all searchwords.

-l --limit
Limits the number of search results to the top [limit] genes found (default: None).

-o --out
Path to the file the results will be saved in, e.g. path/to/directory/results.csv (default: None → just prints results).
For Jupyter Lab / Google Colab: save=True will save the output to the current working directory.


## gget info
Look up gene or transcript Ensembl IDs. 

-id --ens_ids
One or more Ensembl IDs.

-e --expand
Expand returned information (default: False). For genes: add isoform information. For transcripts: add translation and exon information.

-H --homology
Returns homology information of ID (default: False).

-x --xref
Returns information from external references (default: False).

-o --out
Path to the file the results will be saved in, e.g. path/to/directory/results.json (default: None → just prints results).
For Jupyter Lab / Google Colab: save=True will save the output to the current working directory.


## gget seq
Fetch DNA sequences from gene or transcript Ensembl IDs.

-id --ens_ids
One or more Ensembl IDs.

-i --isoforms
If a gene Ensembl ID is passed, this returns sequences of all known transcript isoforms.

-o --out
Path to the file the results will be saved in, e.g. path/to/directory/results.fa (default: None → just prints results).
For Jupyter Lab / Google Colab: save=True will save the output FASTA to the current working directory.

Author: Laura Luebbert
""")