# Copyright 2022 Laura Luebbert

## Packages for gget search
import pandas as pd
import mysql.connector as sql
import time

## Packages for gget ref/info
from bs4 import BeautifulSoup
import requests
import re
import numpy as np
import json

## Packages for gget blast
import logging
# Add and format time stamp in logging messages
logging.basicConfig(format="%(asctime)s %(message)s", datefmt="%d %b %Y %H:%M:%S")
# Using urllib instead of requests here because requests does not 
# allow long queries (queries very long here due to input sequence)
from urllib.request import urlopen
from urllib.parse import urlencode
from urllib.request import Request

## Packages for gget seq
import urllib
from io import StringIO 

## Packages for gget muscle
import os

## To write standard error output
import sys

## Custom functions
from .utils import (
    rest_query,
    ref_species_options,
    gget_species_options,
    parse_blast_ref_page
)
from .compile import (
    compile_muscle,
    MUSCLE_PATH
)
# Constants
from .constants import (
    BLAST_URL,
    BLAST_CLIENT,
    ENSEMBL_REST_API,
    MUSCLE_GITHUB_LINK,
    UNIPROT_REST_API
)

## gget muscle
def muscle(fasta, 
           super5=False, 
           out="default"
          ):
    """
    Align multiple nucleotide or amino acid sequences against each other (using the Muscle v5 algorithm).
    
    Args:
    - fasta
    Path to fasta file containing the sequences to be aligned.
    - super5
    True/False (default: False). 
    If True, align input using Super5 algorithm instead of PPP algorithm to decrease time and memory.
    Use for large inputs (a few hundred sequences).
    - out
    Path to the 'aligned FASTA' (.afa) file the results will be saved in (default: "muscle_results.afa").
        
    Returns alignment results in an "aligned FASTA" (.afa) file.
    """
    # Muscle v5 documentation: https://drive5.com/muscle5

    # Get absolute path to input fasta file
    abs_fasta_path = os.path.abspath(fasta)
    
    # Get absolute path to output .afa file
    if out == "default":
        abs_out_path = os.path.join(os.getcwd(), "muscle_results.afa")
    else:
        abs_out_path = os.path.abspath(out)
    
    # Compile muscle if it is not already compiled
    if os.path.isfile(MUSCLE_PATH) == False:
        # Compile muscle
        compile_muscle()
        
    else:
        logging.warning(
            "MUSCLE already compiled. "
        )
        
    # Define muscle command
    if super5:
        command = f"{MUSCLE_PATH} -super5 {abs_fasta_path} -output {abs_out_path}"
    else:
        command = f"{MUSCLE_PATH} -align {abs_fasta_path} -output {abs_out_path}"
     
    # Run align command  
    logging.warning(
            "MUSCLE aligning... "
        )
    start_time = time.time()
    os.system(command)
    
    logging.warning(
        f"MUSCLE alignment complete. Alignment time: {round(time.time() - start_time, 2)} seconds."
    )
    

## gget info
def info(
    ens_ids, 
    expand=False, 
    homology=False, 
    xref=False, 
    save=False
):
    """
    Fetch gene and transcript metadata using Ensembl IDs.

    Args:
    - ens_ids
    One or more Ensembl IDs to look up (string or list of strings).
    - expand
    Expand returned information (only for genes and transcripts) (default: False). 
    For genes, this adds isoform information. 
    For transcripts, this adds translation and exon information.
    - homology 
    If True, returns homology information of ID (default: False).
    - xref
    If True, returns information from external references (default: False).
    - save
    If True, saves json with query results in current working directory.

    Returns a dictionary/json file containing the requested information about the Ensembl IDs.
    """
    # Define Ensembl REST API server
    server = ENSEMBL_REST_API
    # Define type of returned content from REST
    content_type = "application/json"

    ## Clean up Ensembl IDs
    # If single Ensembl ID passed as string, convert to list
    if type(ens_ids) == str:
        ens_ids = [ens_ids]
    # Remove Ensembl ID version if passed
    ens_ids_clean = []
    for ensembl_ID in ens_ids:
        ens_ids_clean.append(ensembl_ID.split(".")[0])
        
    # Initiate dictionary to save results for all IDs in
    master_dict = {}

    # Query REST APIs from https://rest.ensembl.org/
    for ensembl_ID in ens_ids_clean:
        # Create dict to save query results
        results_dict = {ensembl_ID:{}}

        ## lookup/id/ query: Find the species and database for a single identifier 
        # Define the REST query
        if expand == True:
            query = "lookup/id/" + ensembl_ID + "?" + "expand=1"
        else:
            query = "lookup/id/" + ensembl_ID + "?"
        # Submit query
        try:
            df_temp = rest_query(server, query, content_type)
        # If query returns in an error:
        except:
            # Try submitting query without expand (expand does not work for exons and translation IDs)
            try:
                query = "lookup/id/" + ensembl_ID + "?"
                df_temp = rest_query(server, query, content_type)
            # Raise error if this also did not work
            except:
                sys.stderr.write(
                    f"Ensembl ID {ensembl_ID} not found. "
                    "Please double-check spelling/arguments and try again.\n"
                )
                sys.exit()
            
        ## Delete superfluous entries
        # Delete superfluous entries in general info
        keys_to_delete = ["version", "source", "db_type", "logic_name"]
        for key in keys_to_delete:
            # Pop keys, None -> do not raise an error if key to delete not found
            df_temp.pop(key, None)
        
        # If looking up gene, delete superfluous entries in transcript isoforms info    
        if "Transcript" in df_temp.keys():
            transcript_keys_to_delete = ["assembly_name", "start", "is_canonical", "seq_region_name", "db_type", "source", "strand", "end", "Parent", "species", "version", "logic_name", "Exon", "Translation", "object_type"]
        
            try:
                # More than one isoform present
                for isoform in np.arange(len(df_temp["Transcript"])):
                    for key in transcript_keys_to_delete:
                        df_temp["Transcript"][isoform].pop(key, None)
            except:           
                # Just one isoform present
                for key in transcript_keys_to_delete:
                    df_temp["Transcript"].pop(key, None)
                    
        # If looking up transcript, delete superfluous entries in translation and exon info              
        if "Translation" in df_temp.keys():
            # Delete superfluous entries in Translation info
            translation_keys_to_delete = ["Parent", "species", "db_type", "object_type", "version"]
            
            try:
                # More than one translation present
                for transl in np.arange(len(df_temp["Translation"])):
                    for key in translation_keys_to_delete:
                        df_temp["Translation"][transl].pop(key, None)
            except:
                # Just one translation present
                for key in translation_keys_to_delete:
                    df_temp["Translation"].pop(key, None)
        
        if "Exon" in df_temp.keys():        
            # Delete superfluous entries in Exon info
            exon_keys_to_delete = ["version", "species", "object_type", "db_type", "assembly_name", "seq_region_name", "strand"]
            
            try:
                # More than one exon present
                for exon in np.arange(len(df_temp["Exon"])):
                    for key in exon_keys_to_delete:
                        df_temp["Exon"][exon].pop(key, None)
            except:
                # Just one exon present
                for key in translation_keys_to_delete:
                    df_temp["Exon"].pop(key, None)
            
        ## Add results to main dict
        results_dict[ensembl_ID].update(df_temp)

        ## homology/id/ query: Retrieves homology information (orthologs) by Ensembl gene id
        if homology == True:
            # Define the REST query
            query = "homology/id/" + ensembl_ID + "?"
                
            try:
                # Submit query
                df_temp = rest_query(server, query, content_type)
                # Add results to main dict
                results_dict[ensembl_ID].update({"homology":df_temp["data"][0]["homologies"]})
            except:
                sys.stderr.write(f"No homology information found for {ensembl_ID}.\n")

        ## xrefs/id/ query: Retrieves external reference information by Ensembl gene id
        if xref == True:
            # Define the REST query
            query = "xrefs/id/" + ensembl_ID + "?"

            try:
                # Submit query
                df_temp = rest_query(server, query, content_type)
                # Add results to main dict
                results_dict[ensembl_ID].update({"xrefs":df_temp})
            except:
                sys.stderr.write(f"No external reference information found for {ensembl_ID}.\n")
    
        # Add results to master dict
        master_dict.update(results_dict)

    # Save
    if save == True:
        with open('info_results.json', 'w', encoding='utf-8') as f:
            json.dump(master_dict, f, ensure_ascii=False, indent=4)

    ## Sort nested master_dict alphabetically at all levels
    # Sort IDs keys alphabetically 
    master_dict = {key: value for key, value in sorted(master_dict.items())}
    # Sort ID info level keys alphabetically
    for dict_ens_id in master_dict.keys():
        master_dict[dict_ens_id] = {
            key: value for key, value in sorted(master_dict[dict_ens_id].items())
        }
        # Sort transcript/translation/exon level keys alphabetically
        for trans_id in master_dict[dict_ens_id].keys():
            # Sort if entry is a dict
            if type(master_dict[dict_ens_id][trans_id]) == dict:
                master_dict[dict_ens_id][trans_id] = {
                    key: value
                    for key, value in sorted(master_dict[dict_ens_id][trans_id].items())
                }
            # Sort if entry is a list of dicts
            if type(master_dict[dict_ens_id][trans_id]) == list:
                for index, list_item in enumerate(master_dict[dict_ens_id][trans_id]):
                    master_dict[dict_ens_id][trans_id][index] = {
                        key: value
                        for key, value in sorted(
                            master_dict[dict_ens_id][trans_id][index].items()
                        )
                    }
    
    # Return sorted dictionary containing results
    return master_dict   


## gget search
def search(searchwords, 
           species, 
           seqtype="gene", 
           andor="or", 
           limit=None, 
           save=False
          ):
    """
    Function to query Ensembl for genes based on species and free form search terms. 
    
    Args:
    - searchwords
    The parameter "searchwords" is a list of one or more strings containing free form search terms 
    (e.g.searchwords = ["GABA", "gamma-aminobutyric"]).
    All results that contain at least one of the search terms are returned.
    The search is not case-sensitive.
    - species 
    Species or database. 
    Species can be passed in the format "genus_species", e.g. "homo_sapiens".
    To pass a specific database (e.g. specific mouse strain),
    enter the name of the core database without "/", e.g. 'mus_musculus_dba2j_core_105_1'. 
    All availabale species databases can be found here: http://ftp.ensembl.org/pub/release-105/mysql/
    - seqtype
    Possible entries: "gene" (default) or "transcript". 
    Defines whether genes or transcripts matching the searchwords are returned.
    - andor
    Possible entries: "and", "or"  
    "or": Returns all genes that include at least one of the searchwords in their description (default)  
    "and": Returns only genes that include all of the searchwords in their description 
    - limit
    "Limit" limits the number of search results to the top {limit} genes found.
    - save
    If "save=True", the data frame is saved as a csv in the current directory.
    
    Returns a data frame with the query results.
    """
    start_time = time.time()
    
    ## Check validity or arguments
    # Check if seqtype is valid
    seqtypes = ["gene", "transcript"]
    seqtype = seqtype.lower()
    if seqtype not in seqtypes:
        raise ValueError(
            f"Sequence type specified is {seqtype}. Expected one of {', '.join(seqtypes)}"
        )
        
    # Check if 'andor' arg is valid
    andors = ["and", "or"]
    andor = andor.lower()
    if andor not in andors:
        raise ValueError(
            f"'andor' argument specified as {andor}. Expected one of {', '.join(andors)}"
        )

    ## Get database for specified species
    # Species shortcuts
    if species == "human":
        species = "homo_sapiens"
    if species == "mouse":
        species = "mus_musculus"
    
    # If a specific database is passed with the "/" at the end, remove it
    if "/" in species:
        species = species.split("/")[0]
    
    # In case species was passed with upper case letters
    species = species.lower()

    # Fetch all available databases
    databases = gget_species_options(release=105)
    db = []
    for datab in databases:
        if species in datab:
            db.append(datab)
    
    # Unless an unambigious mouse database is specified, 
    # the standard core database will be used
    if len(db) > 1 and "mus_musculus" in species:
        db = "mus_musculus_core_105_39"

    # Check for ambigious species matches in species other than mouse
    elif len(db) > 1 and "mus_musculus" not in species:
        raise ValueError(
            "Species matches more than one database.\n"
            "Please double-check spelling or pass specific CORE database.\n" 
            "All available databases can be found here:\n"
            "http://ftp.ensembl.org/pub/release-105/mysql/"
            )
    # Raise error if no matching database was found 
    elif len(db) == 0:
        raise ValueError(
            "Species not found in database.\n"
            "Please double-check spelling or pass specific CORE database.\n" 
            "All available databases can be found here:\n"
            "http://ftp.ensembl.org/pub/release-105/mysql/"
            )

    else:
        db = db[0]
        
    logging.warning(f"Fetching results from database: {db}")

    ## Connect to data base
    db_connection = sql.connect(host="ensembldb.ensembl.org", 
                                database=db, 
                                user="anonymous", 
                                password="")

    ## Clean up list of searchwords
    # If single searchword passed as string, convert to list
    if type(searchwords) == str:
        searchwords = [searchwords]
    
    ## Find genes
    for i, searchword in enumerate(searchwords):
        if seqtype == "gene":
            # If limit is specified, fetch only the first {limit} genes for which the searchword appears in the description
            if limit != None:
                query = f"""
                SELECT gene.stable_id, xref.display_label, gene.description, xref.description, gene.biotype
                FROM gene
                LEFT JOIN xref ON gene.display_xref_id = xref.xref_id
                WHERE (gene.description LIKE '%{searchword}%' OR xref.description LIKE '%{searchword}%' OR xref.display_label LIKE '%{searchword}%')
                LIMIT {limit}
                """
            # Else, fetch all genes for which the searchword appears in the description
            else:
                query = f"""
                SELECT gene.stable_id, xref.display_label, gene.description, xref.description, gene.biotype
                FROM gene
                LEFT JOIN xref ON gene.display_xref_id = xref.xref_id
                WHERE (gene.description LIKE '%{searchword}%' OR xref.description LIKE '%{searchword}%' OR xref.display_label LIKE '%{searchword}%')
                """

            # Fetch the search results from the host using the specified query
            df_temp = pd.read_sql(query, con=db_connection)
            # Order by ENSEMBL ID (I am using pandas for this instead of SQL to increase speed)
            df_temp = df_temp.sort_values("stable_id").reset_index(drop=True)

            # If andor="or", keep all results
            if andor == "or":
                # In the first iteration, make the search results equal to the master data frame
                if i == 0:
                    df = df_temp.copy()
                # Add new search results to master data frame
                else:
                    df = pd.concat([df, df_temp])

            # If andor="and", only keep overlap between results
            if andor == "and":
                # In the first iteration, make the search results equal to the master data frame
                if i == 0:
                    df = df_temp.copy()
                # Only keep overlapping results in master data frame
                else:
                    val = np.intersect1d(df["stable_id"], df_temp["stable_id"])
                    df = df[df.stable_id.isin(val)]
        
        if seqtype == "transcript":
            # If limit is specified, fetch only the first {limit} transcripts for which the searchword appears in the description
            if limit != None:
                query = f"""
                SELECT transcript.stable_id, xref.display_label, transcript.description, xref.description, transcript.biotype
                FROM transcript
                LEFT JOIN xref ON transcript.display_xref_id = xref.xref_id
                WHERE (transcript.description LIKE '%{searchword}%' OR xref.description LIKE '%{searchword}%' OR xref.display_label LIKE '%{searchword}%')
                LIMIT {limit}
                """
            # Else, fetch all transcripts for which the searchword appears in the description
            else:
                query = f"""
                SELECT transcript.stable_id, xref.display_label, transcript.description, xref.description, transcript.biotype
                FROM transcript
                LEFT JOIN xref ON transcript.display_xref_id = xref.xref_id
                WHERE (transcript.description LIKE '%{searchword}%' OR xref.description LIKE '%{searchword}%' OR xref.display_label LIKE '%{searchword}%')
                """

            # Fetch the search results from the host using the specified query
            df_temp = pd.read_sql(query, con=db_connection)
            # Order by ENSEMBL ID (I am using pandas for this instead of SQL to increase speed)
            df_temp = df_temp.sort_values("stable_id").reset_index(drop=True)

            # If andor="or", keep all results
            if andor == "or":
                # In the first iteration, make the search results equal to the master data frame
                if i == 0:
                    df = df_temp.copy()
                # Add new search results to master data frame
                else:
                    df = pd.concat([df, df_temp])

            # If andor="and", only keep overlap between results
            if andor == "and":
                # In the first iteration, make the search results equal to the master data frame
                if i == 0:
                    df = df_temp.copy()
                # Only keep overlapping results in master data frame
                else:
                    val = np.intersect1d(df["stable_id"], df_temp["stable_id"])
                    df = df[df.stable_id.isin(val)]

    # Rename columns
    df = df.rename(columns={"stable_id": "ensembl_id", 
                            "display_label":"gene_name",
                            "biotype": "biotype"})
    # Changing description columns name by column index since they were returned with the same name ("description")
    df.columns.values[2] = "ensembl_description"
    df.columns.values[3] = "ext_ref_description"

    # Remove any duplicate search results from the master data frame and reset the index
    df = df.drop_duplicates().reset_index(drop=True)

    # Add URL to gene summary on Ensembl
    df["url"] = "https://uswest.ensembl.org/" + "_".join(db.split("_")[:2]) + "/Gene/Summary?g=" + df["ensembl_id"]
    
    # Print query time and number of genes/transcripts fetched
    logging.warning(f"Query time: {round(time.time() - start_time, 2)} seconds")
    logging.warning(f"Matches found: {len(df)}")
    
    # Save
    if save == True:
        df.to_csv("gget_results.csv", index=False)
    
    # Return data frame
    return df

## gget ref
def ref(species, which="all", release=None, ftp=False, save=False):
    """
    Fetch FTPs for reference genomes and annotations by species.
    
    Args:
    - species
    Defines the species for which the files should be fetched in the format "<genus>_<species>", 
    e.g.species = "homo_sapiens".
    - which
    Defines which results to return. Possible entries are:
    "all" - Returns all links (default).
    Or one or a combination (as a list of strings) of the following:  
    "gtf" - Returns the GTF FTP link and associated information.
    "cdna" - Returns the cDNA FTP link and associated information.
    "dna" - Returns the DNA FTP link and associated information.
    "cds" - Returns the coding sequences corresponding to Ensembl genes. (Does not contain UTR or intronic sequence.)
    "cdrna" - Returns transcript sequences corresponding to non-coding RNA genes (ncRNA).
    "pep" - Returns the protein translations of Ensembl genes.
    - release
    Defines the Ensembl release number from which the files are fetched, e.g. release = 104.
    (Ensembl releases earlier than release 48 are not suupported.)
    By default, the latest Ensembl release is used.
    - FTP
    If True, returns a list containing only the requested URLs instead of the comprehensive json/dictionary.
    - save
    If "save=True", the json containing all results is saved in the current directory. Only works if "returnval='json'".

    Returns a dictionary containing the requested URLs with their respective Ensembl version and release date and time.
    (If FTP=True, returns a list containing only the URLs.)
    """
    
    # Species shortcuts
    if species == "human":
        species = "homo_sapiens"
    if species == "mouse":
        species = "mus_musculus"
        
    # In case species was passed with upper case letters
    species = species.lower()

    ## Find latest Ensembl release
    url = "http://ftp.ensembl.org/pub/"
    html = requests.get(url)
    
    # Raise error if status code not "OK" Response
    if html.status_code != 200:
        raise RuntimeError(f"HTTP response status code {html.status_code}. Please try again.")
    
    soup = BeautifulSoup(html.text, "html.parser")
    # Find all releases
    releases = soup.body.findAll(text=re.compile("release-"))
    # Get release numbers
    rels = []
    for rel in releases:
        rels.append(rel.split("/")[0].split("-")[-1])

    # Find highest release number (= latest release)
    ENS_rel = np.array(rels).astype(int).max()
        
    # If release != None, use user-defined Ensembl release    
    if release != None:
        # Do not allow user-defined release if it is higher than the latest release
        if release > ENS_rel:
            raise ValueError("Defined Ensembl release number cannot be greater than latest release.")
        else:
            ENS_rel = release

    ## Raise error if species not found
    # Find all available species for GTFs for this Ensembl release
    species_list_gtf = ref_species_options('gtf', release=ENS_rel)
    # Find all available species for FASTAs for this Ensembl release
    species_list_dna = ref_species_options('dna', release=ENS_rel) 

    # Find intersection of the two lists 
    # (Only species which have GTF and FASTAs available can continue)
    species_list = list(set(species_list_gtf) & set(species_list_dna))

    if species not in species_list:
        raise ValueError(
            f"Species does not match any available species for Ensembl release {ENS_rel}.\n"
            "Please double-check spelling.\n"
            "$ gget ref --list -> lists out all available species.\n"
            "Combine with [-r] to define specific release (default: latest).\n"
            )
    
    ## Get GTF link for this species and release
    url = f"http://ftp.ensembl.org/pub/release-{ENS_rel}/gtf/{species}/"
    html = requests.get(url)
    
    # Raise error if status code not "OK" Response
    if html.status_code != 200:
        raise RuntimeError(f"HTTP response status code {html.status_code}. Please try again.")
    
    soup = BeautifulSoup(html.text, "html.parser")
    
    # The url can be found under an <a> object tag in the html, 
    # but the date and size do not have an object tag (element=None)
    nones = []
    a_elements = []
    pre = soup.find('pre')
    for element in pre.descendants:
        if element.name == "a":
            a_elements.append(element)
        elif element.name != "a":
            nones.append(element)
    
    # Find the <a> element containing the url
    for i, string in enumerate(a_elements):
        if f"{ENS_rel}.gtf.gz" in string['href']:
            gtf_str = string
            # Get release date and time from <None> elements (since there are twice as many, 2x and +1 to move from string to date)
            gtf_date_size = nones[i*2+1]
            
    gtf_url = f"http://ftp.ensembl.org/pub/release-{ENS_rel}/gtf/{species}/{gtf_str['href']}"
            
    gtf_date = gtf_date_size.strip().split("  ")[0]
    gtf_size = gtf_date_size.strip().split("  ")[-1]

    ## Get cDNA FASTA link for this species and release
    url = f"http://ftp.ensembl.org/pub/release-{ENS_rel}/fasta/{species}/cdna"
    html = requests.get(url)
    
    # Raise error if status code not "OK" Response
    if html.status_code != 200:
        raise RuntimeError(f"HTTP response status code {html.status_code}. Please try again.")
    
    soup = BeautifulSoup(html.text, "html.parser")
    
    # The url can be found under an <a> object tag in the html, 
    # but the date and size do not have an object tag (element=None)
    nones = []
    a_elements = []
    pre = soup.find('pre')
    for element in pre.descendants:
        if element.name == "a":
            a_elements.append(element)
        elif element.name != "a":
            nones.append(element)
            
    # Find the <a> element containing the url       
    for i, string in enumerate(a_elements):
        if "cdna.all.fa" in string['href']:
            cdna_str = string
            # Get release date and time from <None> elements (since there are twice as many, 2x and +1 to move from string to date)
            cdna_date_size = nones[i*2+1]
            
    cdna_url = f"http://ftp.ensembl.org/pub/release-{ENS_rel}/fasta/{species}/cdna/{cdna_str['href']}"
            
    cdna_date = cdna_date_size.strip().split("  ")[0]
    cdna_size = cdna_date_size.strip().split("  ")[-1]
    
    ## Get DNA FASTA link for this species and release
    url = f"http://ftp.ensembl.org/pub/release-{ENS_rel}/fasta/{species}/dna"
    html = requests.get(url)
    
    # Raise error if status code not "OK" Response
    if html.status_code != 200:
        raise RuntimeError(f"HTTP response status code {html.status_code}. Please try again.")
    
    soup = BeautifulSoup(html.text, "html.parser")
    
    # The url can be found under an <a> object tag in the html, 
    # but the date and size do not have an object tag (element=None)
    nones = []
    a_elements = []
    pre = soup.find('pre')
    
    for element in pre.descendants:
        if element.name == "a":
            a_elements.append(element)
        elif element.name != "a":
            nones.append(element)

    # Get primary assembly if available, otherwise toplevel assembly
    dna_str = None
    for i, string in enumerate(a_elements):
        if ".dna.primary_assembly.fa" in string['href']:
            dna_str = string
            dna_search = ".dna.primary_assembly.fa"
            # Get date from non-assigned values (since there are twice as many, 2x and +1 to move from string to date)
            dna_date_size = nones[i*2+1]
            
    # Find the <a> element containing the url        
    if dna_str == None:
        for i, string in enumerate(a_elements):
            if ".dna.toplevel.fa" in string['href']:
                dna_str = string
                dna_search = ".dna.toplevel.fa"
                # Get date from non-assigned values (since there are twice as many, 2x and +1 to move from string to date)
                dna_date_size = nones[i*2+1]
    
    dna_url = f"http://ftp.ensembl.org/pub/release-{ENS_rel}/fasta/{species}/dna/{dna_str['href']}"
      
    dna_date = dna_date_size.strip().split("  ")[0]
    dna_size = dna_date_size.strip().split("  ")[-1]
    # Strip again to remove any extra spaces
    dna_date = dna_date.strip()
    dna_size = dna_size.strip()
    
    ## Get CDS FASTA link for this species and release
    url = f"http://ftp.ensembl.org/pub/release-{ENS_rel}/fasta/{species}/cds"
    html = requests.get(url)
    
    # Raise error if status code not "OK" Response
    if html.status_code != 200:
        raise RuntimeError(f"HTTP response status code {html.status_code}. Please try again.")
    
    soup = BeautifulSoup(html.text, "html.parser")
    
    # The url can be found under an <a> object tag in the html, 
    # but the date and size do not have an object tag (element=None)
    nones = []
    a_elements = []
    pre = soup.find('pre')
    for element in pre.descendants:
        if element.name == "a":
            a_elements.append(element)
        elif element.name != "a":
            nones.append(element)
            
    # Find the <a> element containing the url       
    for i, string in enumerate(a_elements):
        if "cds.all.fa" in string['href']:
            cds_str = string
            # Get release date and time from <None> elements (since there are twice as many, 2x and +1 to move from string to date)
            cds_date_size = nones[i*2+1]
            
    cds_url = f"http://ftp.ensembl.org/pub/release-{ENS_rel}/fasta/{species}/cds/{cds_str['href']}"
            
    cds_date = cds_date_size.strip().split("  ")[0]
    cds_size = cds_date_size.strip().split("  ")[-1]
    
    ## Get ncRNA FASTA link for this species and release
    url = f"http://ftp.ensembl.org/pub/release-{ENS_rel}/fasta/{species}/ncrna"
    html = requests.get(url)
    
    # Raise error if status code not "OK" Response
    if html.status_code != 200:
        raise RuntimeError(f"HTTP response status code {html.status_code}. Please try again.")
    
    soup = BeautifulSoup(html.text, "html.parser")
    
    # The url can be found under an <a> object tag in the html, 
    # but the date and size do not have an object tag (element=None)
    nones = []
    a_elements = []
    pre = soup.find('pre')
    for element in pre.descendants:
        if element.name == "a":
            a_elements.append(element)
        elif element.name != "a":
            nones.append(element)
            
    # Find the <a> element containing the url       
    for i, string in enumerate(a_elements):
        if ".ncrna.fa" in string['href']:
            ncrna_str = string
            # Get release date and time from <None> elements (since there are twice as many, 2x and +1 to move from string to date)
            ncrna_date_size = nones[i*2+1]
            
    ncrna_url = f"http://ftp.ensembl.org/pub/release-{ENS_rel}/fasta/{species}/ncrna/{ncrna_str['href']}"
            
    ncrna_date = ncrna_date_size.strip().split("  ")[0]
    ncrna_size = ncrna_date_size.strip().split("  ")[-1]
    
    ## Get pep FASTA link for this species and release
    url = f"http://ftp.ensembl.org/pub/release-{ENS_rel}/fasta/{species}/pep"
    html = requests.get(url)
    
    # Raise error if status code not "OK" Response
    if html.status_code != 200:
        raise RuntimeError(f"HTTP response status code {html.status_code}. Please try again.")
    
    soup = BeautifulSoup(html.text, "html.parser")
    
    # The url can be found under an <a> object tag in the html, 
    # but the date and size do not have an object tag (element=None)
    nones = []
    a_elements = []
    pre = soup.find('pre')
    for element in pre.descendants:
        if element.name == "a":
            a_elements.append(element)
        elif element.name != "a":
            nones.append(element)
            
    # Find the <a> element containing the url       
    for i, string in enumerate(a_elements):
        if ".pep.all.fa" in string['href']:
            pep_str = string
            # Get release date and time from <None> elements (since there are twice as many, 2x and +1 to move from string to date)
            pep_date_size = nones[i*2+1]
            
    pep_url = f"http://ftp.ensembl.org/pub/release-{ENS_rel}/fasta/{species}/pep/{pep_str['href']}"
            
    pep_date = pep_date_size.strip().split("  ")[0]
    pep_size = pep_date_size.strip().split("  ")[-1]
    
    ## Return results
    # If single which passed as string, convert to list
    if type(which) == str:
        which = [which]

    # Raise error if several values are passed and 'all' is included
    if len(which) > 1 and "all" in which:
        raise ValueError("Parameter 'which' must be 'all', or any one or a combination of the following: 'gtf', 'cdna', 'dna', 'cds', 'ncrna' 'pep'.")

    # If FTP=False, return dictionary/json of specified results
    if ftp == False:
        ref_dict = {species:{}}
        for return_val in which:
            if return_val == "all":
                ref_dict = {
                    species: {
                        "transcriptome_cdna": {
                            "ftp":cdna_url,
                            "ensembl_release": int(ENS_rel),
                            "release_date": cdna_date.split(" ")[0],
                            "release_time": cdna_date.split(" ")[1],
                            "bytes": cdna_size
                        },
                        "genome_dna": {
                            "ftp":dna_url,
                            "ensembl_release": int(ENS_rel),
                            "release_date": dna_date.split(" ")[0],
                            "release_time": dna_date.split(" ")[1],
                            "bytes": dna_size
                        },
                        "annotation_gtf": {
                            "ftp":gtf_url,
                            "ensembl_release": int(ENS_rel),
                            "release_date": gtf_date.split(" ")[0],
                            "release_time": gtf_date.split(" ")[1],
                            "bytes": gtf_size
                        },
                        "coding_seq_cds": {
                            "ftp":cds_url,
                            "ensembl_release": int(ENS_rel),
                            "release_date": cds_date.split(" ")[0],
                            "release_time": cds_date.split(" ")[1],
                            "bytes": cds_size
                        },
                        "non-coding_seq_ncRNA": {
                            "ftp":ncrna_url,
                            "ensembl_release": int(ENS_rel),
                            "release_date": ncrna_date.split(" ")[0],
                            "release_time": ncrna_date.split(" ")[1],
                            "bytes": ncrna_size
                        },
                        "protein_translation_pep": {
                            "ftp":pep_url,
                            "ensembl_release": int(ENS_rel),
                            "release_date": pep_date.split(" ")[0],
                            "release_time": pep_date.split(" ")[1],
                            "bytes": pep_size
                        }
                    }
                }
            elif return_val == "gtf":
                dict_temp = {
                        "annotation_gtf": {
                            "ftp":gtf_url,
                            "ensembl_release": int(ENS_rel),
                            "release_date": gtf_date.split(" ")[0],
                            "release_time": gtf_date.split(" ")[1],
                            "bytes": gtf_size
                    },
                }
                ref_dict[species].update(dict_temp)
            elif return_val == "cdna":
                dict_temp = {
                        "transcriptome_cdna": {
                            "ftp":cdna_url,
                            "ensembl_release": int(ENS_rel),
                            "release_date": cdna_date.split(" ")[0],
                            "release_time": cdna_date.split(" ")[1],
                            "bytes": cdna_size
                    },
                }
                ref_dict[species].update(dict_temp)
            elif return_val == "dna":
                dict_temp = {
                        "genome_dna": {
                            "ftp":dna_url,
                            "ensembl_release": int(ENS_rel),
                            "release_date": dna_date.split(" ")[0],
                            "release_time": dna_date.split(" ")[1],
                            "bytes": dna_size
                    },
                }
                ref_dict[species].update(dict_temp)
            elif return_val == "cds":
                dict_temp = {
                        "coding_seq_cds": {
                            "ftp":cds_url,
                            "ensembl_release": int(ENS_rel),
                            "release_date": cds_date.split(" ")[0],
                            "release_time": cds_date.split(" ")[1],
                            "bytes": cds_size
                    },
                }
                ref_dict[species].update(dict_temp)
            elif return_val == "ncrna":
                dict_temp = {
                        "non-coding_seq_ncRNA": {
                            "ftp":ncrna_url,
                            "ensembl_release": int(ENS_rel),
                            "release_date": ncrna_date.split(" ")[0],
                            "release_time": ncrna_date.split(" ")[1],
                            "bytes": ncrna_size
                    },
                }
                ref_dict[species].update(dict_temp)
            elif return_val == "pep":
                dict_temp = {
                        "protein_translation_pep": {
                            "ftp":pep_url,
                            "ensembl_release": int(ENS_rel),
                            "release_date": pep_date.split(" ")[0],
                            "release_time": pep_date.split(" ")[1],
                            "bytes": pep_size
                    },
                }
                ref_dict[species].update(dict_temp)
            else:
                raise ValueError("Parameter 'which' must be 'all', or any one or a combination of the following: 'gtf', 'cdna', 'dna'.")

        if save == True:
            with open('ref_results.json', 'w', encoding='utf-8') as f:
                json.dump(ref_dict, f, ensure_ascii=False, indent=4)

        sys.stderr.write(f"Fetching from Ensembl release: {ENS_rel}\n")
        return ref_dict
        
    # If FTP==True, return only the specified URLs as a list 
    if ftp == True:
        results = []
        for return_val in which:
            if return_val == "all":
                results.append(gtf_url)
                results.append(cdna_url)
                results.append(dna_url)
            elif return_val == "gtf":
                results.append(gtf_url)
            elif return_val == "cdna":
                results.append(cdna_url)
            elif return_val == "dna":
                results.append(dna_url)
            elif return_val == "cds":
                results.append(cds_url)
            elif return_val == "ncrna":
                results.append(ncrna_url)
            elif return_val == "pep":
                results.append(pep_url)
            else:
                raise ValueError("Parameter 'which' must be 'all', or any one or a combination of the following: 'gtf', 'cdna', 'dna', 'cds', 'ncrna', 'pep'.")

        if save == True:
            file = open("ref_results.txt", "w")
            for element in results:
                file.write(element + "\n")
            file.close()

        return results
    
# gget seq helper function (here and not in .utils because it uses gget info)
def get_uniprot_seqs(server, ensembl_ids):
    """
    Retrieve UniProt sequences based on Ensemsbl identifiers.

    Args:
    - server
    Link to UniProt REST API server.
    - ensembl_ids: 
    One or more Ensembl IDs (string or list of strings).
    All passed Ensembl IDs must be of the same object type (gene or transcript).

    Returns data frame with UniProt ID, gene name, organism, sequence, sequence length, and query ID.
    """
    
    # If a single UniProt ID is passed as string, convert to list
    if type(ensembl_ids) == str:
        ensembl_ids = [ensembl_ids]
        
    # Check whether first Ensembl ID is gene or transcript using gget info
    # The first Ensembl ID defines the input type for all IDs
    id_ = ensembl_ids[0]
    id_type = info(id_)[id_]["object_type"]
    if id_type == "Gene":
        input_id = "ENSEMBL_ID"
        # sys.stderr.write("Ensembl ID recognized as gene.\n")
    elif id_type == "Transcript":
        input_id = "ENSEMBL_TRS_ID"
        # sys.stderr.write("Ensembl ID recognized as transcript.\n")
    else: 
        sys.stderr.write(f"Object type of Ensembl ID {id_} not recognized.")
        sys.exit()
    
    # Define query arguments
    # Columns documentation: https://www.uniprot.org/help/uniprotkb%5Fcolumn%5Fnames
    # from/to IDs documentation: https://www.uniprot.org/help/api_idmapping
    query_args = {
        "from": input_id,
        "to": "ACC",
        "format": "tab",
        "query": " ".join(ensembl_ids),
        "columns": "id,genes,organism,sequence,length",
    }
    # Reformat query arguments
    query_args = urllib.parse.urlencode(query_args)
    query_args = query_args.encode("ascii")
    
    # Submit query to UniProt server
    request = urllib.request.Request(server, query_args)
    
    # Read and clean up results
    with urllib.request.urlopen(request) as response:
        res = response.read()
        
    try:
        df = pd.read_csv(StringIO(res.decode("utf-8")), sep="\t")
        # Rename columns
        df.columns = ["UniProt ID", "Gene name", "Organism", "Sequence", "Sequence Length", "Query"]
        # Split rows if two different UniProt IDs for a single query ID are returned
        df = df.assign(Query=df["Query"].str.split(",")).explode("Query")

        return df
    
    except:
        sys.stderr.write(f"No UniProt transcript sequences were found for these Ensembl ID(s).\n")
        sys.exit()

# gget seq
def seq(ens_ids,
        seqtype="gene",
        isoforms=False, 
        save=False,
       ):
    """
    Fetch nucleotide or amino acid sequence (FASTA) of a gene 
    (and all isoforms) or transcript by Ensembl ID. 

    Args:
    - ens_ids
    One or more Ensembl IDs (passed as string or list of strings).
    - seqtype
    'gene' (default) or 'transcript'. 
    Defines whether nucleotide or amino acid sequences are returned.
    Nucleotide sequences are fetched from the Ensembl REST API server.
    Amino acid sequences are fetched from the UniProt REST API server.
    - isoforms
    If True, returns the sequences of all known transcripts (for gene IDs only)
    (default: False).
    - save
    If True: Save output FASTA to current directory.
    
    Returns a list (or FASTA file if 'save=True') containing the sequences of the Ensembl IDs.
    """
    
    ## Clean up arguments
    # Check if seqtype is valid
    seqtypes = ["gene", "transcript"]
    seqtype = seqtype.lower()
    if seqtype not in seqtypes:
        raise ValueError(
            f"Sequence type specified is {seqtype}. Expected one of {', '.join(seqtypes)}"
        )
        
    # Clean up Ensembl IDs
    # If single Ensembl ID passed as string, convert to list
    if type(ens_ids) == str:
        ens_ids = [ens_ids]
    # Remove Ensembl ID version if passed
    ens_ids_clean = []
    for ensembl_ID in ens_ids:
        ens_ids_clean.append(ensembl_ID.split(".")[0])
    
    ## Fetch nucleotide sequences
    if seqtype == "gene":
        # Define Ensembl REST API server
        server = ENSEMBL_REST_API
        # Define type of returned content from REST
        content_type = "application/json"

        ## Clean up Ensembl IDs
        # If single Ensembl ID passed as string, convert to list
        if type(ens_ids) == str:
            ens_ids = [ens_ids]
        # Remove Ensembl ID version if passed
        ens_ids_clean = []
        for ensembl_ID in ens_ids:
            ens_ids_clean.append(ensembl_ID.split(".")[0])

        # Initiate dictionary to save results for all IDs in
        master_dict = {}

        # Query REST APIs from https://rest.ensembl.org/
        for ensembl_ID in ens_ids_clean:
            # Create dict to save query results
            results_dict = {ensembl_ID:{}}

            ## Fetch nucleotide sequece
            # sequence/id/ query: Request sequence by stable identifier

            # If isoforms false, just fetch sequences of passed Ensembl ID
            if isoforms == False:
                # Define the REST query
                query = "sequence/id/" + ensembl_ID + "?"
                # Submit query
                df_temp = rest_query(server, query, content_type)

                # Delete superfluous entries
                keys_to_delete = ["query", "id", "version", "molecule"]
                for key in keys_to_delete:
                    # Pop keys, None -> do not raise an error if key to delete not found
                    df_temp.pop(key, None)

                # Add results to main dict
                results_dict[ensembl_ID].update({"seq":df_temp})

            # If isoforms true, fetch sequences of isoforms instead
            if isoforms == True:
                # Get transcripts using gget info
                info_dict = info(ensembl_ID, expand=True)

                # If this is a gene, get the sequence of all isoforms using gget info
                try:
                    info_dict[ensembl_ID]["Transcript"]

                    # If only one transcript present
                    try:
                        transcipt_id = info_dict[ensembl_ID]["Transcript"]["id"]

                        # Define the REST query
                        query = "sequence/id/" + transcipt_id + "?"
                        # Submit query
                        df_temp = rest_query(server, query, content_type)

                        # Delete superfluous entries
                        keys_to_delete = ["query", "version", "molecule"]
                        for key in keys_to_delete:
                            # Pop keys, None -> do not raise an error if key to delete not found
                            df_temp.pop(key, None)

                        # Add results to main dict
                        results_dict[ensembl_ID].update({"transcript":df_temp})

                    # If more than one transcript present    
                    except:
                        for isoform in np.arange(len(info_dict[ensembl_ID]["Transcript"])):
                            transcipt_id = info_dict[ensembl_ID]["Transcript"][isoform]["id"]

                            # Define the REST query
                            query = "sequence/id/" + transcipt_id + "?"
                            # Submit query
                            df_temp = rest_query(server, query, content_type)

                            # Delete superfluous entries
                            keys_to_delete = ["query", "version", "molecule"]
                            for key in keys_to_delete:
                                # Pop keys, None -> do not raise an error if key to delete not found
                                df_temp.pop(key, None)

                            # Add results to main dict
                            results_dict[ensembl_ID].update({f"transcript{isoform}":df_temp})
                except:
                    # If the gget info results do not include "Trancript" 
                    # (i.e. because input is not a gene), raise value error
                    sys.stderr.write("Isoform option not availabale for non-gene Ensembl IDs.")
                    sys.exit()

            # Add results to master dict
            master_dict.update(results_dict)

        # Build FASTA file
        fasta = []
        for ens_ID in master_dict:
            for key in master_dict[ens_ID].keys():
                if key == 'seq':
                    fasta.append(">" + ens_ID + " " + master_dict[ens_ID][key]['desc'])
                    fasta.append(master_dict[ens_ID][key]['seq'])
                else:
                    fasta.append(">" + master_dict[ens_ID][key]['id'] + " " + master_dict[ens_ID][key]['desc'])
                    fasta.append(master_dict[ens_ID][key]['seq'])
    
    ## Fetch amino acid sequences from UniProt
    if seqtype == "transcript":
        if isoforms == False:
            # Fetch amino acid sequences of the passed Ensembl iD from UniProt REST API
            df_uniprot = get_uniprot_seqs(UNIPROT_REST_API, ens_ids_clean)
        
        if isoforms == True:
            # Get the transcript isoforms using gget info
            info_dict = info(ensembl_ID, expand=True)

            # If this is a gene, get the sequence of all isoforms
            try:
                info_dict[ensembl_ID]["Transcript"]

                # If only one transcript present
                try:
                    # Get the ID of the transcript from the gget info results
                    transcipt_id = info_dict[ensembl_ID]["Transcript"]["id"]
                    # Fetch amino acid sequence from the UniProt REST API
                    df_uniprot = get_uniprot_seqs(UNIPROT_REST_API, transcipt_id)

                # If more than one transcript present    
                except:
                    # Get the ID of all transcripts from the gget info results
                    transcipt_ids = []
                    for isoform_idx in np.arange(len(info_dict[ensembl_ID]["Transcript"])):
                        transcipt_ids.append(info_dict[ensembl_ID]["Transcript"][isoform_idx]["id"])
                    
                    # Fetch amino acid sequences of all isoforms from the UniProt REST API
                    df_uniprot = get_uniprot_seqs(UNIPROT_REST_API, transcipt_ids)
                
            # If the gget info results do not include "Trancript" 
            # (i.e. because input is not a gene), raise value error
            except:
                sys.stderr.write("Isoform option not availabale for non-gene Ensembl IDs.")
                sys.exit()
        
        # Build UniProt results FASTA file
        fasta = []
        for uniprot_id, query_ensembl_id, gene_name, organism, uniprot_seq in zip(
            df_uniprot["UniProt ID"].values,
            df_uniprot["Query"].values,
            df_uniprot["Gene name"].values,
            df_uniprot["Organism"].values,
            df_uniprot["Sequence"].values
        ):
            fasta.append(
                ">"
                + "uniprot_id: " + uniprot_id
                + " ensembl_id: " + query_ensembl_id 
    #                 + " gene_name:" + gene_name 
    #                 + " organism:" + organism
            )
            fasta.append(uniprot_seq)
            
    # Save
    if save == True:
        file = open("seq_results.fa", "w")
        for element in fasta:
            file.write(element + "\n")
        file.close()
    
    return fasta   

## gget blast
def blast(
    sequence,
    program="default",
    database="default",
    ncbi_gi=False,
    descriptions=500,
    alignments=500,
    hitlist_size=50,
    expect=10.0,
    low_comp_filt=False,
    megablast=True,
    verbose=True,
):
    """
    BLAST a nucleotide or amino acid sequence against any BLAST DB.
    Args:
     - sequence       Sequence (str) or path to fasta file containing one sequence.
     - program        'blastn', 'blastp', 'blastx', 'tblastn', or 'tblastx'. 
                      Default: 'blastn' for nucleotide sequences; 'blastp' for amino acid sequences.
     - database       'nt', 'nr', 'refseq_rna', 'refseq_protein', 'swissprot', 'pdbaa', or 'pdbnt'. 
                      Default: 'nt' for nucleotide sequences; 'nr' for amino acid sequences.
                      More info on BLAST databases: https://ncbi.github.io/blast-cloud/blastdb/available-blastdbs.html
     - ncbi_gi        True/False whether to return NCBI GI identifiers. Default False.
     - descriptions   int or None. Limit number of descriptions to return. Default 500.
     - alignments     int or None. Limit number of alignments to return. Default 500.
     - hitlist_size   int or None. Limit number of hits to return. Default 50.
     - expect         float or None. An expect value cutoff. Default 10.0.
     - low_comp_filt  True/False whether to apply low complexity filter. Default False.
     - megablast      True/False whether to use the MegaBLAST algorithm (blastn only). Default True.
     - verbose        True/False whether to print progress information. Default True.

    NCBI server rule: 
    Run scripts weekends or between 9 pm and 5 am Eastern time
    on weekdays if more than 50 searches will be submitted.

    This function does not check the validity of the arguments
    and passes the values to the server as is. 
    """
    # Server rules:
    # 1. Do not contact the server more often than once every 10 seconds.
    # 2. Do not poll for any single RID more often than once a minute.
    # 3. Use the URL parameter email and tool, so that the NCBI
    #    can contact you if there is a problem.
    # 4. Run scripts weekends or between 9 pm and 5 am Eastern time
    #    on weekdays if more than 50 searches will be submitted.
    # Reference: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=DeveloperInfo
    
    # Please note that NCBI uses the new Common URL API for BLAST searches
    # on the internet (http://ncbi.github.io/blast-cloud/dev/api.html). Thus,
    # some of the arguments used by this function are not (or are no longer)
    # officially supported by NCBI. Although they are still functioning, this
    # may change in the future.

    # Define server URL and client
    url = BLAST_URL
    client = BLAST_CLIENT

    ## Clean up arguments
    # If the path to a fasta file was provided instead of a nucleotide sequence,
    # read the file and extract the first sequence
    if ".fa" in sequence:
        # Read the fasta file and append its entries to 'fasta_values'
        from Bio.SeqIO import FastaIO
        fasta_values = []
        with open(sequence) as handle:
            for value in FastaIO.SimpleFastaParser(handle):
                fasta_values.append(value)
        # Set the first sequence from the fasta file as 'sequence'
        sequence = fasta_values[0][1]
        
        if len(fasta_values) > 2:
            logging.warning(
                "More than one sequence was passed in the .fa file. "
                "Only the first sequence will be submitted to BLAST."
            )
    
    ## Set program and database
    
    # Convert program and database to lower case
    program = program.lower()
    database = database.lower()
    # Valid program and database options
    programs = ["blastn", "blastp", "blastx", "tblastn", "tblastx"]
    dbs = ["nt", "nr", "refseq_rna", "refseq_protein", "swissprot", "pdbaa", "pdbnt"]
    
    # If user does not specify the program, 
    # check if a nulceotide or amino acid sequence was passed
    if program == "default":
        # Set of all possible nucleotides and amino acids
        nucleotides = set("ATGC")
        amino_acids = set("ARNDCQEGHILKMFPSTWYVBZ") 
        
        # If sequence is a nucleotide sequence, set program to blastn
        if set(sequence) <= nucleotides:
            program = "blastn"
            
            # Set database to nt (unless user specified another database)
            if database == "default":
                database = "nt"
                if verbose == True:
                    logging.warning("Sequence recognized as nucleotide sequence. "
                                    "BLAST will use program 'blastn' with database 'nt'.")
            else:
                # Check if the user specified database is valid
                if database not in dbs:
                    raise ValueError(
                        f"Database specified is {database}. Expected one of {', '.join(dbs)}"
                    )
                else:
                    if verbose == True:
                        logging.warning("Sequence recognized as nucleotide sequence. "
                                        "BLAST will use program 'blastn' with user-specified database.")
        # If sequence is an amino acid sequence, set program to blastp        
        elif set(sequence) <= amino_acids:
            program = "blastp"
            
            # Set database to nr (unless user specified another database)
            if database == "default":
                database = "nr"
                if verbose == True:
                    logging.warning("Sequence recognized as amino acid sequence. "
                                    "BLAST will use program 'blastp' with database 'nr'.")
            else:
                # Check if the user specified database is valid
                if database not in dbs:
                    raise ValueError(
                        f"Database specified is {database}. Expected one of {', '.join(dbs)}"
                    )
                else:
                    if verbose == True:
                        logging.warning("Sequence recognized as amino acid sequence. "
                                        "BLAST will use program 'blastp' with user-specified database.")
        else:
            raise ValueError(f"""
                Sequence not automatically recognized as a nucleotide or amino acid sequence.
                Please specify 'program' and 'database'.
                Program options: {', '.join(programs)} 
                Database options:  {', '.join(dbs)} 
                """)
    # Check if the user specified program is valid
    else:
        if program not in programs:
            raise ValueError(
                f"Program specified is {program}. Expected one of {', '.join(programs)}"
            )

    ## Translate filter and ncbi_gi arguments
    if low_comp_filt == False:
        low_comp_filt = None
    else:
        low_comp_filt = "T"

    if ncbi_gi == False:
        ncbi_gi = None
    else:
        ncbi_gi = "T"

    if megablast == False:
        megablast = None
    else:
        megablast = "on"

    ## Submit search
    #  The following code was partly adapted from the Biopython BLAST NCBIWWW project written
    #  by Jeffrey Chang (Copyright 1999), Brad Chapman, and Chris Wroe distributed under the
    #  Biopython License Agreement and BSD 3-Clause License
    #  https://github.com/biopython/biopython/blob/171697883aca6894f8367f8f20f1463ce7784d0c/LICENSE.rst
                             
    # Args for the PUT command
    put_args = [
        ("PROGRAM", program),
        ("DATABASE", database),
        ("QUERY", sequence),
        ("NCBI_GI", ncbi_gi),
        ("DESCRIPTIONS", descriptions),
        ("ALIGNMENTS", alignments),
        ("HITLIST_SIZE", hitlist_size),
        ("EXPECT", expect),
        ("FILTER", low_comp_filt),
        ("MEGABLAST", megablast),
        ("CMD", "Put"),
    ]

    # Define query
    put_query = [x for x in put_args if x[1] is not None]
    put_message = urlencode(put_query).encode()

    # Submit search to server
    request = Request(url, put_message, {"User-Agent": client})
    handle = urlopen(request)

    ## Fetch Request ID (RID) and estimated time to completion (RTOE)
    RID, RTOE = parse_blast_ref_page(handle)
        
    # Wait for search to complete
    # (At least 11 seconds to comply with server rule 1)
    if RTOE < 11:
        # Communicate RTOE
        if verbose == True:
            logging.warning(f"BLAST initiated. Estimated time to completion: 11 seconds.")
        time.sleep(11)
    else:
        # Communicate RTOE
        if verbose == True:
            logging.warning(f"BLAST initiated. Estimated time to completion: {RTOE} seconds.")  
        time.sleep(int(RTOE))

    ## Poll server for status and fetch search results
    # Args for the GET command
    get_args = [
        ("RID", RID),
        ("ALIGNMENTS", alignments),
        ("DESCRIPTIONS", descriptions),
        ("HITLIST_SIZE", hitlist_size),
        ("FORMAT_TYPE", "HTML"),
        ("NCBI_GI", ncbi_gi),
        ("CMD", "Get"),
    ]
    get_query = [x for x in get_args if x[1] is not None]
    get_message = urlencode(get_query).encode()

    ## Poll NCBI until the results are ready
    searching = True
    i = 0
    while searching:
        if i > 0:
            # Sleep for 61 seconds if first fetch was not succesful
            # to comply with server rules
            time.sleep(61)

        # Query for search status
        request = Request(url, get_message, {"User-Agent": client})
        handle = urlopen(request)
        results = handle.read().decode()
        
        # Fetch search status
        i = results.index("Status=")
        j = results.index("\n", i)
        status = results[i + len("Status=") : j].strip()

        if status == "WAITING":
            if verbose == True:
                logging.warning("BLASTING...")
            i += 1
            continue

        if status == "FAILED":
            raise ValueError(f"Search {RID} failed; please report to blast-help@ncbi.nlm.nih.gov.")

        if status == "UNKNOWN":
            raise ValueError(f"NCBI status {status}. Search {RID} expired.")

        if status == "READY":
            if verbose == True:
                logging.warning("Retrieving results...")
            # Stop search
            searching = False

            ## Return results
#             if verbose == True:
#                 logging.warning("BLAST complete.")
                
            # Parse HTML results
            soup = BeautifulSoup(results, "html.parser")
            # Get the descriptions table
            dsc_table = soup.find(lambda tag: tag.name=='table' and tag.has_attr('id') and tag['id']=="dscTable") 
            results_df = pd.read_html(str(dsc_table))[0]
            # Drop the first column
            results_df = results_df.iloc[: , 1:]
            
            return results_df

        else:
            raise ValueError(f"""
            Something unexpected happened. \n
            Search {RID} possibly failed; please report to blast-help@ncbi.nlm.nih.gov\n
            or post an issue on Github: https://github.com/lauraluebbert/gget\n
            """
            )
