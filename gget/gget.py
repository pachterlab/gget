# Packages for gget search
import pandas as pd
import mysql.connector as sql
import time

# Packages for gget ref/info
from bs4 import BeautifulSoup
import requests
import re
import numpy as np
import json

# Package to write standard error output
import sys

# Packages for use from terminal
# from . import __version__
import argparse
import os
from tabulate import tabulate

version = "0.0.16"

def rest_query(server, query, content_type):
    """
    Function to query a 

    Parameters:
    - server
    Serve to query.
    - Query
    Query that is passed to server.
    - content_type
    Contect type requested from server.

    Returns server output.
    """

    r = requests.get(
        server + query, 
        headers={ "Content-Type" : content_type}
    )

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    if content_type == 'application/json':
        return r.json()
    else:
        return r.text
        
def gget_species_options(release=105):
    """
    Function to find all available species core databases for gget.

    Parameters:
    - release
    Ensembl release for which the databases are fetched.

    Returns list of available core databases.
    """
    # Find all available databases
    url = f"http://ftp.ensembl.org/pub/release-{release}/mysql/"
    html = requests.get(url)
    
    # Raise error if status code not "OK" Response
    if html.status_code != 200:
        raise RuntimeError(f"HTTP response status code {html.status_code}. Please try again.")
        
    soup = BeautifulSoup(html.text, "html.parser")

    # Return list of all available databases
    databases = []
    for subsoup in soup.body.findAll('a'):
        if "core" in subsoup["href"]:
            databases.append(subsoup["href"].split("/")[0])
    
    return databases
    
    
def ref_species_options(which, release=None):
    """
    Function to find all available species for gget ref.

    Parameters:
    - release
    Ensembl release for which available species should be fetched.
    - which
    Which type of FTP. Possible entries: 'dna', 'cdna', 'gtf'.

    Returns list of available species.
    """
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

    # Find all available species for this release and FTP type
    if which == "gtf":
        url = f"http://ftp.ensembl.org/pub/release-{ENS_rel}/gtf/"
    if which == "dna" or which == "cdna":
        url = f"http://ftp.ensembl.org/pub/release-{ENS_rel}/fasta/"
    html = requests.get(url)
    
    # Raise error if status code not "OK" Response
    if html.status_code != 200:
        raise RuntimeError(f"HTTP response status code {html.status_code}. Please try again.")
    
    soup = BeautifulSoup(html.text, "html.parser")

    sps = []
    for subsoup in soup.body.findAll('a'):
        sps.append(subsoup["href"].split("/")[0])

    species_list = sps[1:]
    
    # Return list of all available species
    return species_list

## gget info
def info(ens_ids, expand=False, homology=False, xref=False, save=False):
    """
    Looks up information about Ensembl IDs.

    Parameters:
    - ens_ids
    One or more Ensembl IDs to look up (passed as string or list of strings).
    -expand
    Expand returned information (default: False). For genes: add isoform information. For transcripts: add translation and exon information.
    - homology
    If True, returns homology information of ID (default: False).
    - xref
    If True, returns information from external references (default: False).
    -save
    If True, saves json with query results in current working directory.

    Returns a dictionary/json file containing the requested information about the Ensembl IDs.
    """
    # Define Ensembl REST API server
    server = "http://rest.ensembl.org/"
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
        # Raise error if ID not found
        except:
            sys.stderr.write(f"Ensembl ID {ensembl_ID} not found. Please double-check spelling.\n")
            
        ## Delete superfluous entries
        # Delete superfluous entries in general info
        keys_to_delete = ["version", "source", "db_type", "logic_name", "id"]
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
            # Submit query
            df_temp = rest_query(server, query, content_type)
                
            # Add results to main dict
            try:
                results_dict[ensembl_ID].update({"homology":df_temp["data"][0]["homologies"]})
            except:
                sys.stderr.write(f"No homology information found for {ensembl_ID}.\n")

        ## xrefs/id/ query: Retrieves external reference information by Ensembl gene id
        if xref == True:
            # Define the REST query
            query = "xrefs/id/" + ensembl_ID + "?"
            # Submit query
            df_temp = rest_query(server, query, content_type)

            # Add results to main dict
            try:
                results_dict[ensembl_ID].update({"xrefs":df_temp})
            except:
                sys.stderr.write(f"No external reference information found for {ensembl_ID}.\n")
    
        # Add results to master dict
        master_dict.update(results_dict)

    # Save
    if save == True:
        with open('info_results.json', 'w', encoding='utf-8') as f:
            json.dump(master_dict, f, ensure_ascii=False, indent=4)

    # Return dictionary containing results
    return master_dict   


## gget search
def search(searchwords, species, d_type="gene", andor="or", limit=None, save=False):
    """
    Function to query Ensembl for genes based on species and free form search terms. 
    
    Parameters:
    - searchwords
    The parameter "searchwords" is a list of one or more strings containing free form search terms 
    (e.g.searchwords = ["GABA", "gamma-aminobutyric"]).
    All results that contain at least one of the search terms are returned.
    The search is not case-sensitive.
    - species
    Species or database. 
    Species can be passed in the format 'genus_species', e.g. 'homo_sapiens'.
    To pass a specific database (e.g. specific mouse strain),
    enter the name of the core database without "/", e.g. 'mus_musculus_dba2j_core_105_1'. 
    All availabale species databases can be found here: http://ftp.ensembl.org/pub/release-105/mysql/
    -d_type
    "gene" (default) or "transcript". 
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
        
    sys.stderr.write(f"Fetching results from database: {db}\n")

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
        if d_type == "gene":
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
        
        if d_type == "transcript":
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
    df = df.rename(columns={"stable_id": "Ensembl_ID", 
                            "display_label":"Gene_name",
                            "biotype": "Biotype"})
    # Changing description columns name by column index since they were returned with the same name ("description")
    df.columns.values[2] = "Ensembl_description"
    df.columns.values[3] = "Ext_ref_description"

    # Remove any duplicate search results from the master data frame and reset the index
    df = df.drop_duplicates().reset_index(drop=True)

    # Add URL to gene summary on Ensembl
    df["URL"] = "https://uswest.ensembl.org/" + "_".join(db.split("_")[:2]) + "/Gene/Summary?g=" + df["Ensembl_ID"]
    
    # Print query time and number of genes fetched
    sys.stderr.write(f"Query time: {round(time.time() - start_time, 2)} seconds\n")
    sys.stderr.write(f"Genes fetched: {len(df)}\n")
    
    # Save
    if save == True:
        df.to_csv("gget_results.csv", index=False)
    
    # Return data frame
    return df

## gget ref
def ref(species, which="all", release=None, ftp=False, save=False):
    """
    Function to fetch GTF and FASTA (cDNA and DNA) URLs from the Ensemble FTP site.
    
    Parameters:
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
    
    ## Return results
    # If single which passed as string, convert to list
    if type(which) == str:
        which = [which]

    # Raise error if several values are passed and 'all' is included
    if len(which) > 1 and "all" in which:
        raise ValueError("Parameter 'which' must be 'all', or any one or a combination of the following: 'gtf', 'cdna', 'dna'.")

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
            else:
                raise ValueError("Parameter 'which' must be 'all', or any one or a combination of the following: 'gtf', 'cdna', 'dna'.")

        if save == True:
            file = open("ref_results.txt", "w")
            for element in results:
                file.write(element + "\n")
            file.close()

        return results
    
# gget seq
def seq(ens_ids, isoforms=False, save=False):
    """
    Fetch DNA sequences from gene or transcript Ensembl IDs. 

    Parameters:
    - ens_ids
    One or more Ensembl IDs (passed as string or list of strings).
    - isoforms
    If true: If a gene Ensembl ID is passed, this returns sequences of all known transcript isoforms.
    - save
    If True: Save output FASTA to current directory.
    
    Returns a FASTA file containing sequences of the Ensembl IDs.
    """
    # Define Ensembl REST API server
    server = "http://rest.ensembl.org/"
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
        
        ## SEQUENCE
        # sequence/id/ query: Request sequence by stable identifier
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
                pass
            
#         ## OVERLAP
#         # overlap/id/ query: Retrieves features (e.g. genes, transcripts, variants and more) that overlap a region defined by the given identifier.
#         query = "overlap/id/" + ensembl_ID + "?" + "feature=gene"
#         # Submit query
#         df_temp = rest_query(server, query, content_type)
#         print(df_temp)


#         ## MAPPING
#         # map/cdna/:id/:region: Convert from cDNA coordinates to genomic coordinates. Output reflects forward orientation coordinates as returned from the Ensembl API.
#         query = "map/cdna/" + ensembl_ID + "/100..300" + "?"
#         # Submit query
#         df_temp = rest_query(server, query, content_type)
#         print(df_temp)      
        
        
#         # REGULATORY
#         regulatory/species/:species/id/:id query: Returns a RegulatoryFeature given its stable ID
        
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
                
    
    # Save
    if save == True:
        file = open("seq_results.fa", "w")
        for element in fasta:
            file.write(element + "\n")
        file.close()

    # Return dictionary containing results
    return fasta   


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
    
    
def main():
    """
    Function containing argparse parsers and arguments to allow the use of gget from the terminal.
    """
    # Define parent parser 
    parent_parser = argparse.ArgumentParser(description=f"gget v{version}", add_help=False)
    # Initiate subparsers
    parent_subparsers = parent_parser.add_subparsers(dest="command")
    # Define parent (not sure why I need both parent parser and parent, but otherwise it does not work)
    parent = argparse.ArgumentParser(add_help=False)
    
    # Add custom help argument to parent parser
    parent_parser.add_argument(
            "-h","--help",
            action="store_true",
            help="Print manual. Recommendation: pipe into less by running 'gget -h | less'"
    )
    # Add custom version argument to parent parser
    parent_parser.add_argument(
            "-v","--version",
            action="store_true",
            help="Print version."
    )
    
    ## gget ref subparser
    parser_ref = parent_subparsers.add_parser(
        "ref",
        parents=[parent],
        description="Fetch FTP links for a specific species from Ensemble.",
        help="Fetch FTP links for a specific species from Ensemble.",
        add_help=True
        )
    # ref parser arguments
    parser_ref.add_argument(
        "-s", "--species", 
        default=None,
        type=str,
        help="Species for which the FTPs will be fetched, e.g. homo_sapiens."
    )
    # ref parser arguments
    parser_ref.add_argument(
        "-l", "--list", 
        default=None, 
        action="store_true",
        required=False,
        help="List out all available species."
    )
    parser_ref.add_argument(
        "-w", "--which", 
        default="all", 
        type=str,
        nargs='+',
        required=False,
        help=("Defines which results to return.\n" 
              "Possible entries are:\n"
              "'all' - Returns GTF, cDNA, and DNA links and associated info (default).\n" 
              "Or one or a combination of the following:\n"  
              "'gtf' - Returns the GTF FTP link and associated info.\n" 
              "'cdna' - Returns the cDNA FTP link and associated info.\n"
              "'dna' - Returns the DNA FTP link and associated info."
             )
        )
    parser_ref.add_argument(
        "-r", "--release",
        default=None,  
        type=int, 
        required=False,
        help="Ensemble release the FTPs will be fetched from, e.g. 104 (default: latest Ensembl release).")
    parser_ref.add_argument(
        "-ftp", "--ftp",  
        default=False, 
        action="store_true",
        required=False,
        help="Return only the FTP link instead of a json.")
    parser_ref.add_argument(
        "-d", "--download",  
        default=False, 
        action="store_true",
        required=False,
        help="Download FTPs to the current directory using wget.")
    parser_ref.add_argument(
        "-o", "--out",
        type=str,
        required=False,
        help=(
            "Path to the json file the results will be saved in, e.g. path/to/directory/results.json.\n" 
            "Default: None (just prints results)."
        )
    )

    ## gget search subparser
    parser_gget = parent_subparsers.add_parser(
        "search",
         parents=[parent],
         description="Query Ensembl for genes based on species and free form search terms.", 
         help="Query Ensembl for genes based on species and free form search terms.",
         add_help=True
         )
    # Search parser arguments
    parser_gget.add_argument(
        "-sw", "--searchwords", 
        type=str, 
        nargs="+",
        required=True, 
        help="One or more free form searchwords for the query, e.g. gaba, nmda."
    )
    parser_gget.add_argument(
        "-s", "--species",
        type=str,  
        required=True, 
        help="Species to be queried, e.g. homo_sapiens."
    )
    parser_gget.add_argument(
        "-t", "--d_type",
        choices=["gene", "transcript"],
        default="gene",
        type=str,  
        required=False, 
        help=(
            "'gene': Returns genes that match the searchwords. (default).\n"
            "'transcript': Returns transcripts that match the searchwords. \n"
        )
    )
    parser_gget.add_argument(
        "-ao", "--andor",
        choices=["and", "or"],
        default="or",
        type=str,  
        required=False, 
        help=(
            "'or': Gene descriptions must include at least one of the searchwords (default).\n"
            "'and': Only return genes whose descriptions include all searchwords.\n"
        )
    )
    parser_gget.add_argument(
        "-l", "--limit",
        type=int, 
        default=None,
        required=False,
        help="Limits the number of results, e.g. 10 (default: None)."
    )
    parser_gget.add_argument(
        "-o", "--out",
        type=str,
        required=False,
        help=(
            "Path to the json file the results will be saved in, e.g. path/to/directory/results.json.\n" 
            "Default: None (just prints results)."
        )
    )
    
    ## gget info subparser
    parser_info = parent_subparsers.add_parser(
        "info",
        parents=[parent],
        description="Look up information about Ensembl IDs.", 
        help="Look up information about Ensembl IDs.",
        add_help=True
        )
    # info parser arguments
    parser_info.add_argument(
        "-id", "--ens_ids", 
        type=str,
        nargs="+",
        required=True, 
        help="One or more Ensembl IDs."
    )
    parser_info.add_argument(
        "-e", "--expand", 
        default=False, 
        action="store_true",
        required=False, 
        help="Expand returned information (default: False). For genes: add isoform information. For transcripts: add translation and exon information."
    )
    parser_info.add_argument(
        "-H", "--homology", 
        default=False, 
        action="store_true",
        required=False, 
        help="Returns homology information of ID (default: False)."
    )
    parser_info.add_argument(
        "-x", "--xref", 
        default=False, 
        action="store_true",
        required=False, 
        help="Returns information from external references (default: False)."
    )
    parser_info.add_argument(
        "-o", "--out",
        type=str,
        required=False,
        help=(
            "Path to the json file the results will be saved in, e.g. path/to/directory/results.json.\n" 
            "Default: None (just prints results)."
        )
    )
    
    ## gget seq subparser
    parser_seq = parent_subparsers.add_parser(
        "seq",
        parents=[parent],
        description="Look up DNA sequences from Ensembl IDs.", 
        help="Look up DNA sequences from Ensembl IDs.",
        add_help=True
        )
    # info parser arguments
    parser_seq.add_argument(
        "-id", "--ens_ids", 
        type=str,
        nargs="+",
        required=True, 
        help="One or more Ensembl IDs."
    )
    parser_seq.add_argument(
        "-i", "--isoforms", 
        default=False, 
        action="store_true",
        required=False, 
        help="If searching a gene ID, returns sequences of all known transcripts (default: False)."
    )
    parser_seq.add_argument(
        "-o", "--out",
        type=str,
        required=False,
        help=(
            "Path to the FASTA file the results will be saved in, e.g. path/to/directory/results.fa.\n" 
            "Default: None (just prints results)."
        )
    )
    
    
    ## Show help when no arguments are given
    if len(sys.argv) == 1:
        parent_parser.print_help(sys.stderr)
        sys.exit(1)

    args = parent_parser.parse_args()

    ### Define return values
    ## Help return
    if args.help:
        help_()
        
    ## Version return
    if args.version:        
        print(f"gget version: {version}")
        
    ## ref return
    if args.command == "ref":
        # If list flag but no release passed, return all available species for latest release
        if args.list and args.release is None:
                # Find all available species for GTFs for this Ensembl release
                species_list_gtf = ref_species_options('gtf')
                # Find all available species for FASTAs for this Ensembl release
                species_list_dna = ref_species_options('dna') 

                # Find intersection of the two lists 
                # (Only species which have GTF and FASTAs available can continue)
                species_list = list(set(species_list_gtf) & set(species_list_dna))
                
                # Print available species list
                print(species_list)
                
        # If list flag and release passed, return all available species for this release
        if args.list and args.release:
                # Find all available species for GTFs for this Ensembl release
                species_list_gtf = ref_species_options('gtf', release=args.release)
                # Find all available species for FASTAs for this Ensembl release
                species_list_dna = ref_species_options('dna', release=args.release) 

                # Find intersection of the two lists 
                # (Only species which have GTF and FASTAs available can continue)
                species_list = list(set(species_list_gtf) & set(species_list_dna))
                
                # Print available species list
                print(species_list)
        
        # Raise error if neither species nor list flag passed
        if args.species is None and args.list is None:
            parser_ref.error("\n\nThe following arguments are required to fetch FTPs: -s/--species, e.g. '-s homo_sapiens'\n\n"
                             "gget ref --list -> lists out all available species. " 
                             "Combine with [-r] to define specific Ensembl release (default: latest release).")
        
        ## Clean up 'which' entry if passed
        if type(args.which) != str:
            which_clean = []
            # Split by comma (spaces are automatically split by nargs:"+")
            for which in args.which:
                which_clean.append(which.split(","))
            # Flatten which_clean
            which_clean_final = [item for sublist in which_clean for item in sublist]   
            # Remove empty strings resulting from split
            while("" in which_clean_final):
                which_clean_final.remove("")   
        else:
            which_clean_final = args.which

        if args.species:
            
            # Query Ensembl for requested FTPs using function ref
            ref_results = ref(args.species, which_clean_final, args.release, args.ftp)

            # Print or save list of URLs (ftp=True)
            if args.ftp == True:
                # Save in specified directory if -o specified
                if args.out:
                    directory = "/".join(args.out.split("/")[:-1])
                    if directory != "":
                        os.makedirs(directory, exist_ok=True)
                    file = open(args.out, "w")
                    for element in ref_results:
                        file.write(element + "\n")
                    file.close()
                    sys.stderr.write(
                        f"\nResults saved as {args.out}.\n"
                    )
                    
                    if args.download == True:
                        # Download list of URLs
                        for link in ref_results:
                            command = "wget " + link
                            os.system(command)
                    else:
                        sys.stderr.write(
                            "To download the FTPs to the current directory, add flag [-d].\n"
                        )
                
                # Print results if no directory specified
                else:
                    # Print results
                    results = " ".join(ref_results)
                    print(results)
                    sys.stderr.write(
                        "\nTo save these results, use flag '-o' in the format: '-o path/to/directory/results.txt'.\n"
                    )
                    
                    if args.download == True:
                        # Download list of URLs
                        for link in ref_results:
                            command = "wget " + link
                            os.system(command)
                    else:
                        sys.stderr.write(
                            "To download the FTPs to the current directory, add flag [-d].\n"
                        )
                    
            # Print or save json file (ftp=False)
            else:
                # Save in specified directory if -o specified
                if args.out:
                    directory = "/".join(args.out.split("/")[:-1])
                    if directory != "":
                        os.makedirs(directory, exist_ok=True)
                    with open(args.out, 'w', encoding='utf-8') as f:
                        json.dump(ref_results, f, ensure_ascii=False, indent=4)
                    sys.stderr.write(
                        f"\nResults saved as {args.out}.\n"
                    )
                    
                    if args.download == True:
                        # Download the URLs from the dictionary
                        for link in ref_results:
                            for sp in ref_results:
                                for ftp_type in ref_results[sp]:
                                    link = ref_results[sp][ftp_type]['ftp']
                                    command = "wget " + link
                                    os.system(command)    
                    else:
                        sys.stderr.write(
                            "To download the FTPs to the current directory, add flag [-d].\n"
                        )
                    
                # Print results if no directory specified
                else:
                    print(json.dumps(ref_results, ensure_ascii=False, indent=4))
                    sys.stderr.write(
                        "\nTo save these results, use flag '-o' in the format: '-o path/to/directory/results.json'.\n"
                    )
                    
                    if args.download == True:
                        # Download the URLs from the dictionary
                        for link in ref_results:
                            for sp in ref_results:
                                for ftp_type in ref_results[sp]:
                                    link = ref_results[sp][ftp_type]['ftp']
                                    command = "wget " + link
                                    os.system(command)
                    else:
                        sys.stderr.write(
                            "To download the FTPs to the current directory, add flag [-d].\n"
                        )
        
    ## search return
    if args.command == "search":
        
        ## Clean up args.searchwords
        sw_clean = []
        # Split by comma (spaces are automatically split by nargs:"+")
        for sw in args.searchwords:
            sw_clean.append(sw.split(","))
        # Flatten which_clean
        sw_clean_final = [item for sublist in sw_clean for item in sublist]   
        # Remove empty strings resulting from split
        while("" in sw_clean_final) :
            sw_clean_final.remove("")  
        
        # Query Ensembl for genes based on species and searchwords using function search
        gget_results = search(sw_clean_final, 
                              args.species,
                              d_type=args.d_type,
                              andor=args.andor, 
                              limit=args.limit)
        
        # Save in specified directory if -o specified
        if args.out:
            directory = "/".join(args.out.split("/")[:-1])
            if directory != "":
                os.makedirs(directory, exist_ok=True)
            gget_results.to_csv(args.out, index=False)
            sys.stderr.write(f"\nResults saved as {args.out}.\n")
        
        # Print results if no directory specified
        else:
            print(tabulate(gget_results, headers = 'keys', tablefmt = 'plain'))
            sys.stderr.write("\nTo save these results, use flag '-o' in the format: '-o path/to/directory/results.csv'.\n")
            
    ## info return
    if args.command == "info":

        ## Clean up args.ens_ids
        ids_clean = []
        # Split by comma (spaces are automatically split by nargs:"+")
        for id_ in args.ens_ids:
            ids_clean.append(id_.split(","))
        # Flatten which_clean
        ids_clean_final = [item for sublist in ids_clean for item in sublist]   
        # Remove empty strings resulting from split
        while("" in ids_clean_final) :
            ids_clean_final.remove("")  

        # Look up requested Ensembl IDs
        info_results = info(ids_clean_final, expand=args.expand, homology=args.homology, xref=args.xref)

        # Print or save json file
        # Save in specified directory if -o specified
        if args.out:
            directory = "/".join(args.out.split("/")[:-1])
            if directory != "":
                os.makedirs(directory, exist_ok=True)
            with open(args.out, 'w', encoding='utf-8') as f:
                json.dump(info_results, f, ensure_ascii=False, indent=4)
            sys.stderr.write(f"\nResults saved as {args.out}.\n")
        # Print results if no directory specified
        else:
            print(json.dumps(info_results, ensure_ascii=False, indent=4))
            sys.stderr.write("\nTo save these results, use flag '-o' in the format: '-o path/to/directory/results.json'.\n")
            
    ## seq return
    if args.command == "seq":

        ## Clean up args.ens_ids
        ids_clean = []
        # Split by comma (spaces are automatically split by nargs:"+")
        for id_ in args.ens_ids:
            ids_clean.append(id_.split(","))
        # Flatten which_clean
        ids_clean_final = [item for sublist in ids_clean for item in sublist]   
        # Remove empty strings resulting from split
        while("" in ids_clean_final) :
            ids_clean_final.remove("")  

        # Look up requested Ensembl IDs
        seq_results = seq(ids_clean_final, isoforms=args.isoforms)

        # Save in specified directory if -o specified
        if args.out:
            directory = "/".join(args.out.split("/")[:-1])
            if directory != "":
                os.makedirs(directory, exist_ok=True)
            file = open(args.out, "w")
            for element in seq_results:
                file.write(element + "\n")
            file.close()
            sys.stderr.write(
                f"\nResults saved as {args.out}.\n"
            )
            
        # Print results if no directory specified
        else:
            print(seq_results)
            sys.stderr.write(
                "\nTo save these results, use flag '-o' in the format: '-o path/to/directory/results.fa'.\n"
            )
    
# Python interpreter to run main()
if __name__ == '__main__':
    main()
