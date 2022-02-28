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
def info(ens_ids, homology=False, xref=False, save=False):
    """
    Looks up information about Ensembl IDs.

    Parameters:
    - ens_ids
    One or more Ensembl IDs to look up (passed as string or list of strings).
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
        query = "lookup/id/" + ensembl_ID + "?"
        # Submit query
        try:
            df_temp = utils.rest_query(server, query, content_type)
        # Raise error if ID not found
        except:
            sys.stderr.write(f"Ensembl ID {ensembl_ID} not found. Please double-check spelling.\n")
        # Delete superfluous entries
        try:
            del df_temp["version"], df_temp["source"], df_temp["db_type"], df_temp["logic_name"], df_temp["id"]
        except:
            continue

        # Add results to main dict
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
def search(searchwords, species, andor="or", limit=None, save=False):
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
        # If limit is specified, fetch only the first {limit} genes for which the searchword appears in the description
        if limit != None:
            query = f"""
            SELECT gene.stable_id, gene.description, xref.description, gene.biotype
            FROM gene
            LEFT JOIN xref ON gene.display_xref_id = xref.xref_id
            WHERE (gene.description LIKE '%{searchword}%' OR xref.description LIKE '%{searchword}%')
            LIMIT {limit}
            """
        # Else, fetch all genes for which the searchword appears in the description
        else:
            query = f"""
            SELECT gene.stable_id, gene.description, xref.description, gene.biotype
            FROM gene
            LEFT JOIN xref ON gene.display_xref_id = xref.xref_id
            WHERE (gene.description LIKE '%{searchword}%' OR xref.description LIKE '%{searchword}%')
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
                            "biotype": "Biotype"})
    # Changing description columns name by column index since they were returned with the same name ("description")
    df.columns.values[1] = "Ensembl_description"
    df.columns.values[2] = "Ext_ref_description"

    # Remove any duplicate search results from the master data frame and reset the index
    df = df.drop_duplicates().reset_index(drop=True)

    # Find name of gene using info function and add to df
    gene_names = []
    for ens_id in df["Ensembl_ID"].values:
        try:
            gene_names.append(info(ens_id)[ens_id]["display_name"])
        # If no gene name is found, add "None" instead
        except KeyError:
            gene_names.append(None)
    # Add gene names to df
    df["Gene_name"] = gene_names

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

    ## Find latest Ensembl release
    url = "http://ftp.ensembl.org/pub/"
    html = requests.get(url)
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
        if f"{ENS_rel}.gtf.gz" in string.text:
            gtf_str = string
            
    gtf_url = f"http://ftp.ensembl.org/pub/release-{ENS_rel}/gtf/{species}/{gtf_str['href']}"
            
    # Get release date and time of this GTF link by searching the <None> elements
    for i, string in enumerate(nones):
        if f"{ENS_rel}.gtf.gz" in string.text:
            gtf_date_size = nones[i+1]
            
    gtf_date = gtf_date_size.strip().split("  ")[0]
    gtf_size = gtf_date_size.strip().split("  ")[-1]

    ## Get cDNA FASTA link for this species and release
    url = f"http://ftp.ensembl.org/pub/release-{ENS_rel}/fasta/{species}/cdna"
    html = requests.get(url)
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
        if "cdna.all.fa" in string.text:
            cdna_str = string
            
    cdna_url = f"http://ftp.ensembl.org/pub/release-{ENS_rel}/fasta/{species}/cdna/{cdna_str['href']}"
    
    # Get release date and time of this url by searching the <None> elements
    for i, string in enumerate(nones):
        if "cdna.all.fa" in string.text:
            cdna_date_size = nones[i+1]
            
    cdna_date = cdna_date_size.strip().split("  ")[0]
    cdna_size = cdna_date_size.strip().split("  ")[-1]
    
    ## Get DNA FASTA link for this species and release
    url = f"http://ftp.ensembl.org/pub/release-{ENS_rel}/fasta/{species}/dna"
    html = requests.get(url)
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
    for string in a_elements:
        if ".dna.primary_assembly.fa" in string.text:
            dna_str = string
            dna_search = ".dna.primary_assembly.fa"
            
    # Find the <a> element containing the url        
    if dna_str == None:
        for string in a_elements:
            if ".dna.toplevel.fa" in string.text:
                dna_str = string
                dna_search = ".dna.toplevel.fa"
        
    dna_url = f"http://ftp.ensembl.org/pub/release-{ENS_rel}/fasta/{species}/dna/{dna_str['href']}"
    
    # Get release date and time of this url by searching the <None> elements          
    for i, string in enumerate(nones):
        if dna_search in string.text:
            dna_date_size = nones[i+1]  
      
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

    
# Python interpreter to run main()
if __name__ == '__main__':
    main.main()
