## File description:
# gget.py contains all of the primary functions for the gget package (supporting functions are stored in utils.py).

## Import packages
# Import gget main and utils
import main
import utils

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
            df_temp = utils.rest_query(server, query, content_type)
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
            df_temp = utils.rest_query(server, query, content_type)
                
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
            df_temp = utils.rest_query(server, query, content_type)

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
    databases = utils.gget_species_options(release=105)
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
    species_list_gtf = utils.ref_species_options('gtf', release=ENS_rel)
    # Find all available species for FASTAs for this Ensembl release
    species_list_dna = utils.ref_species_options('dna', release=ENS_rel) 

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
        df_temp = utils.rest_query(server, query, content_type)
        
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
                    df_temp = utils.rest_query(server, query, content_type)
                    
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
                        df_temp = utils.rest_query(server, query, content_type)
                        
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
#         df_temp = utils.rest_query(server, query, content_type)
#         print(df_temp)


#         ## MAPPING
#         # map/cdna/:id/:region: Convert from cDNA coordinates to genomic coordinates. Output reflects forward orientation coordinates as returned from the Ensembl API.
#         query = "map/cdna/" + ensembl_ID + "/100..300" + "?"
#         # Submit query
#         df_temp = utils.rest_query(server, query, content_type)
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
    
# Python interpreter to run main()
if __name__ == '__main__':
    main.main()
