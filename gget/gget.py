# Packages for gget search
import pandas as pd
import mysql.connector as sql
import time

# Packages for ref
from bs4 import BeautifulSoup
import requests
import re
import numpy as np
import json

# Packages for use from terminal
import argparse
import sys
import os

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

def spy(ens_ids, seq=False, homology=False, xref=False, save=False):
    """
    Looks up information about Ensembl IDs.

    Parameters:
    - ens_ids
    One or more Ensembl IDs to look up (passed as string or list of strings).
    - seq
    If True, returns bp sequence of gene (or parent gene if transcript ID passed) (default: False).
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

    master_dict = {}

    # If single Ensembl ID passed as string, convert to list
    if type(ens_ids) == str:
        ens_ids = [ens_ids]

    # Query REST APIs from https://rest.ensembl.org/
    for ensembl_ID in ens_ids:
        # Create dict to save query results
        results_dict = {ensembl_ID:{}}

        ## lookup/id/ query: Find the species and database for a single identifier 
        # Define the REST query
        query = "lookup/id/" + ensembl_ID + "?"
        # Submit query
        try:
            df_temp = rest_query(server, query, content_type)
        # Raise error if ID not found
        except:
            raise ValueError(f"Ensembl ID {ensembl_ID} not found. Please double-check spelling.")
        # Delete superfluous entries
        try:
            del df_temp["version"], df_temp["source"], df_temp["db_type"], df_temp["logic_name"], df_temp["id"]
        except:
            continue

        # Add results to main dict
        results_dict[ensembl_ID].update(df_temp)

        ## sequence/id/ query: Request sequence by stable identifier
        if seq == True:
            # Define the REST query
            query = "sequence/id/" + ensembl_ID + "?"
            # Submit query
            df_temp = rest_query(server, query, content_type)

            # Add results to main dict
            results_dict[ensembl_ID].update({"seq":df_temp["seq"]})

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
                print("No homology information found for this ID.")

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
                print("No external reference information found for this ID.")
    
        # Add results to master dict
        master_dict.update(results_dict)

    # Save
    if save == True:
        with open('spy_results.json', 'w', encoding='utf-8') as f:
            json.dump(master_dict, f, ensure_ascii=False, indent=4)

    # Return dictionary containing results
    return master_dict   

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
    soup = BeautifulSoup(html.text, "html.parser")

    # Return list of all available databases
    databases = []
    for subsoup in soup.body.findAll('a'):
        if "core" in subsoup["href"]:
            databases.append(subsoup["href"].split("/")[0])
    
    return databases

def search(searchwords, species, limit=None, save=False):
    """
    Function to query Ensembl for genes based on species and free form search terms. 
    
    Parameters:
    - searchwords
    The parameter "searchwords" is a list of one or more strings containing free form search terms 
    (e.g.searchwords = ["GABA", "gamma-aminobutyric acid"]).
    All results that contain at least one of the search terms are returned.
    The search is not case-sensitive.
    - species
    Species or database. 
    Species can be passed in the format 'genus_species', e.g. 'homo_sapiens'.
    To pass a specific database (e.g. specific mouse strain),
    enter the name of the core database without "/", e.g. 'mus_musculus_dba2j_core_105_1'. 
    All availabale species databases can be found here: http://ftp.ensembl.org/pub/release-105/mysql/
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
            "Species matches more than one database. "
            "Please double-check spelling or pass specific database. " 
            "All available databases can be found here: "
            "http://ftp.ensembl.org/pub/release-105/mysql/"
            )
    # Raise error if no matching database was found 
    elif len(db) == 0:
        raise ValueError(
            "Species not found in database. "
            "Please double-check spelling or pass specific database. " 
            "All available databases can be found here: "
            "http://ftp.ensembl.org/pub/release-105/mysql/"
            )

    else:
        db = db[0]
        
    print(f"Fetching results from database: {db}")

    ## Connect to data base
    db_connection = sql.connect(host='ensembldb.ensembl.org', 
                                database=db, 
                                user='anonymous', 
                                password='')

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

        # Fetch the search results form the host using the specified query
        df_temp = pd.read_sql(query, con=db_connection)
        # Order by ENSEMBL ID (I am using pandas for this instead of SQL to increase speed)
        df_temp = df_temp.sort_values("stable_id").reset_index(drop=True)

        # In the first iteration, make the search results equal to the master data frame
        if i == 0:
            df = df_temp.copy()
        # Add new search results to master data frame
        else:
            df = pd.concat([df, df_temp])

    # Rename columns
    df = df.rename(columns={"stable_id": "Ensembl_ID", 
                            "biotype": "Biotype"})
    # Changing description columns name by column index since they were returned with the same name ("description")
    df.columns.values[1] = "Ensembl_description"
    df.columns.values[2] = "Ext_ref_description"

    # Remove any duplicate search results from the master data frame and reset the index
    df = df.drop_duplicates().reset_index(drop=True)

    # Find name of gene using spy function and add to df
    gene_names = []
    for ens_id in df["Ensembl_ID"].values:
        try:
            gene_names.append(spy(ens_id)[ens_id]["display_name"])
        # If no gene name is found, add "None" instead
        except KeyError:
            gene_names.append(None)
    # Add gene names to df
    df["Gene_name"] = gene_names

    # Add URL to gene summary on Ensembl
    df["URL"] = "https://uswest.ensembl.org/" + "_".join(db.split("_")[:2]) + "/Gene/Summary?g=" + df["Ensembl_ID"]
    
    # Print query time and number of genes fetched
    print(f"Query time: {round(time.time() - start_time, 2)} seconds")
    print(f"Genes fetched: {len(df)}")
    
    # Save
    if save == True:
        df.to_csv("gget_results.csv", index=False)
    
    # Return data frame
    return df


def ref_species_options(release, which):
    """
    Function to find all available species for gget ref.

    Parameters:
    - release
    Ensembl release for which available species should be fetched.
    - which
    Which type of FTP. Possible entries: 'dna', 'cdna', 'gtf'.

    Returns list of available species.
    """

    # Find all available species for this release and FTP type
    if which == "gtf":
        url = f"http://ftp.ensembl.org/pub/release-{release}/gtf/"
    if which == "dna" or which == "cdna":
        url = f"http://ftp.ensembl.org/pub/release-{release}/fasta/"
    html = requests.get(url)
    soup = BeautifulSoup(html.text, "html.parser")

    sps = []
    for subsoup in soup.body.findAll('a'):
        sps.append(subsoup["href"].split("/")[0])

    species_list = sps[1:]
    
    # Return list of all available species
    return species_list

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
    species_list_gtf = ref_species_options(ENS_rel, 'gtf')
    # Find all available species for FASTAs for this Ensembl release
    species_list_dna = ref_species_options(ENS_rel, 'dna') 

    # Find intersection of the two lists 
    # (Only species which have GTF and FASTAs available can continue)
    species_list = list(set(species_list_gtf) & set(species_list_dna))

    if species not in species_list:
        raise ValueError(
            f"Species does not match any available species for Ensembl release {ENS_rel}. "
            f"All available species are: {species_list} "
            "Please double-check spelling. "
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

        print(f"Fetching from Ensembl release: {ENS_rel}")
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
 
def main():
    """
    Function containing argparse parsers and arguments to allow the use of gget from the terminal.
    """
    # Define parent parser 
    parent_parser = argparse.ArgumentParser(description="gget parent parser")
    # Initiate subparsers
    parent_subparsers = parent_parser.add_subparsers(dest="command")
    # Define parent (not sure why I need both parent parser and parent, but otherwise it does not work)
    parent = argparse.ArgumentParser(add_help=False)
    # Add debug argument to parent parser
    parent.add_argument(
            '--debug',
            action='store_true',
            help='Print debug info'
        )
    
    ## gget ref subparser
    parser_ref = parent_subparsers.add_parser(
        "ref",
        parents=[parent],
        description="Fetch FTP links for a specific species from Ensemble.",
        add_help=False
        )
    # ref parser arguments
    parser_ref.add_argument(
        "-s", "--species", 
        type=str,
        metavar="",
        required=True, 
        help="Species for which the FTPs will be fetched, e.g. homo_sapiens."
    )
    parser_ref.add_argument(
        "-w", "--which", 
        default="all", 
        type=str,
        nargs='+',
        required=False,
        help=("Defines which results to return." 
              "Possible entries are:"
              "'all' - Returns GTF, cDNA, and DNA links and associated info (default)." 
              "Or one or a combination of the following:"  
              "'gtf' - Returns the GTF FTP link and associated info." 
              "'cdna' - Returns the cDNA FTP link and associated info."
              "'dna' - Returns the DNA FTP link and associated info."
             )
        )
    parser_ref.add_argument(
        "-r", "--release",
        default=None,  
        type=int, 
        required=False,
        metavar="",
        help="Ensemble release the FTPs will be fetched from, e.g. 104 (default: latest Ensembl release).")
    parser_ref.add_argument(
        "-ftp", "--ftp",  
        default=False, 
        action='store_true',
        required=False,
        help="If True: return only the FTP link instead of a json.")
    parser_ref.add_argument(
        "-o", "--out",
        type=str,
        required=False,
        metavar="",
        help=(
            "Path to the json file the results will be saved in, e.g. path/to/directory/results.json." 
            "Default: None (just prints results)."
        )
    )

    ## gget search subparser
    parser_gget = parent_subparsers.add_parser(
        "search",
         parents=[parent],
         description="Query Ensembl for genes based on species and free form search terms.", 
         add_help=False
         )
    # Search parser arguments
    parser_gget.add_argument(
        "-sw", "--searchwords", 
        type=str, 
        nargs="+",
        required=True, 
        metavar="",  
        help="One or more free form searchwords for the query, e.g. gaba, nmda."
    )
    parser_gget.add_argument(
        "-s", "--species",
        type=str,  
        required=True, 
        metavar="",
        help="Species to be queried, e.g. homo_sapiens."
    )
    parser_gget.add_argument(
        "-l", "--limit", 
        type=int, 
        required=False,
        metavar="",
        help="Limits the number of results, e.g. 10 (default: None)."
    )
    parser_gget.add_argument(
        "-o", "--out",
        type=str,
        required=False,
        help=(
            "Path to the csv file the results will be saved in, e.g. path/to/directory/results.csv." 
            "Default: None (just prints results)."
        )
    )
    
    ## gget spy subparser
    parser_spy = parent_subparsers.add_parser(
        "spy",
        parents=[parent],
        description="Look up information about Ensembl IDs.", 
        add_help=False
        )
    # spy parser arguments
    parser_spy.add_argument(
        "-id", "--ens_ids", 
        type=str,
        nargs="+",
        required=True, 
        metavar="",    # Cleans up help message
        help="One or more Ensembl IDs."
    )
    parser_spy.add_argument(
        "-seq", "--seq",
        default=False, 
        action='store_true',
        required=False, 
        help="Returns bp sequence of gene (or parent gene if transcript ID passed) (default: False)."
    )
    parser_spy.add_argument(
        "-H", "--homology",
        default=False, 
        action='store_true',
        required=False, 
        help="Returns homology information of ID (default: False)."
    )
    parser_spy.add_argument(
        "-x", "--xref",
        default=False, 
        action='store_true',
        required=False, 
        help="Returns information from external references (default: False)."
    )
    parser_spy.add_argument(
        "-o", "--out",
        type=str,
        required=False,
        help=(
            "Path to the json file the results will be saved in, e.g. path/to/directory/results.json." 
            "Default: None (just prints results)."
        )
    )
    
    ## Show help when no arguments are given
    if len(sys.argv) == 1:
        parent_parser.print_help(sys.stderr)
        sys.exit(1)

    args = parent_parser.parse_args()

    ### Define return values
    ## Debug return
    if args.debug:
        print("debug: " + str(args))
        
    ## ref return
    if args.command == "ref":
        
        ## Clean up 'which' entry if passed
        if type(args.which) != str:
            which_clean = []
            # Split by comma (spaces are automatically split by nargs:"+")
            for which in args.which:
                which_clean.append(which.split(","))
            # Flatten which_clean
            which_clean_final = [item for sublist in which_clean for item in sublist]   
            # Remove empty strings resulting from split
            while("" in which_clean_final) :
                which_clean_final.remove("")   
        else:
            which_clean_final = args.which

        # Query Ensembl for requested FTPs using function ref
        ref_results = ref(args.species, which_clean_final, args.release, args.ftp)

        # Print or save list of URLs
        if args.ftp==True:
            if args.out:
                os.makedirs("/".join(args.out.split("/")[:-1]), exist_ok=True)
                file = open(args.out, "w")
                for element in ref_results:
                    file.write(element + "\n")
                file.close()
                print(f"Results saved as {args.out}.")

            else:
                print(" ".join(ref_results))
        
        # Print or save json file
        else:
            # Save in specified directory if -o specified
            if args.out:
                os.makedirs("/".join(args.out.split("/")[:-1]), exist_ok=True)
                with open(args.out, 'w', encoding='utf-8') as f:
                    json.dump(ref_results, f, ensure_ascii=False, indent=4)
                print(f"Results saved as {args.out}.")
            # Print results if no directory specified
            else:
                print(json.dumps(ref_results, ensure_ascii=False, indent=4))
                print("To save these results, use flag '-o' in the format: '-o path/to/directory/results.json'.")
        
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
        gget_results = search(sw_clean_final, args.species, limit=args.limit)
        
        # Save in specified directory if -o specified
        if args.out:
            os.makedirs("/".join(args.out.split("/")[:-1]), exist_ok=True)
            gget_results.to_csv(args.out, index=False)
            print(f"Results saved as {args.out}.")
        
        # Print results if no directory specified
        else:
            print(gget_results)
            print("To save these results, use flag '-o' in the format: '-o path/to/directory/results.csv'.")
            
    ## spy return
    if args.command == "spy":

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
        spy_results = spy(ids_clean_final, args.seq, args.homology, args.xref)

        # Print or save json file
        # Save in specified directory if -o specified
        if args.out:
            os.makedirs("/".join(args.out.split("/")[:-1]), exist_ok=True)
            with open(args.out, 'w', encoding='utf-8') as f:
                json.dump(spy_results, f, ensure_ascii=False, indent=4)
            print(f"Results saved as {args.out}.")
        # Print results if no directory specified
        else:
            print(json.dumps(spy_results, ensure_ascii=False, indent=4))
            print("To save these results, use flag '-o' in the format: '-o path/to/directory/results.json'.")

if __name__ == '__main__':
    main()
