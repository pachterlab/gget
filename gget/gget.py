# Packages for gget
import pandas as pd
import mysql.connector as sql
import time

# Packages for ftpget
from bs4 import BeautifulSoup
import requests
import re
import numpy as np

# Packages for use from terminal
import argparse
import sys
import os

def gget(searchwords, species, limit=None, save=False):
    """
    Function to query Ensembl for genes based on species and free form search terms. 
    
    Parameters:
    - searchwords
    The parameter "searchwords" is a list of one or more strings containing free form search terms 
    (e.g.searchwords = ["GABA", "gamma-aminobutyric acid"]).
    All results that contain at least one of the search terms are returned.
    The search is not case-sensitive.
    
    - species
    Possible entries for species are:
    "homo_sapiens" (or "human")
    "mus_musculus" (or "mouse")
    "taeniopygia_guttata" (or "zebra finch")
    "caenorhabditis_elegans" (or "roundworm")
    Note: If you would like to access results from a different species, or restrain your results to a certain mouse strain, 
    you can instead enter the core database as the "species" variable (e.g. species = "rattus_norvegicus_core_105_72").
    You can find all availabale species databases here: http://ftp.ensembl.org/pub/release-105/mysql/
    
    - limit
    "Limit" limits the number of search results to the top {limit} genes found.
    
    - save
    If "save=True", the data frame is saved as a csv in the current directory.
    
    Returns a data frame with the query results.
    """
    start_time = time.time()
    
    if species == "caenorhabditis_elegans" or species == "roundworm":
        db = "caenorhabditis_elegans_core_105_269"   
    elif species == "homo_sapiens" or species == "human":
        db = "homo_sapiens_core_105_38"
    elif species == "mus_musculus" or species == "mouse":
        db = "mus_musculus_core_105_39"
    elif species == "taeniopygia_guttata" or species == "zebra_finch":
        db = "taeniopygia_guttata_core_105_12"
    else: 
        db = species
        
    print(f"Results fetched from database: {db}")

    db_connection = sql.connect(host='ensembldb.ensembl.org', 
                                database=db, 
                                user='anonymous', 
                                password='')
    
    # If single searchword passed as string, convert to list
    if type(searchwords) == str:
        searchwords = [searchwords]
        
    # For human and mouse, the gene name is saved in gene_attrib.value where gene_attrib.attrib_type_id = 4
    if db == "homo_sapiens_core_105_38" or db == "mus_musculus_core_105_39":
        for i, searchword in enumerate(searchwords):
            # If limit is specified, fetch only the first {limit} genes for which the searchword appears in the description
            if limit != None:
                query = f"""
                SELECT gene.stable_id, gene_attrib.value, gene.description, xref.description, gene.biotype
                FROM gene
                LEFT JOIN xref ON gene.display_xref_id = xref.xref_id
                LEFT JOIN gene_attrib ON gene.gene_id = gene_attrib.gene_id
                WHERE (gene_attrib.attrib_type_id = 4) 
                AND (gene_attrib.value LIKE '%{searchword}%' OR 
                gene.description LIKE '%{searchword}%' OR 
                xref.description LIKE '%{searchword}%')
                LIMIT {limit}
                """
            # Else, fetch all genes for which the searchword appears in the description
            else:
                query = f"""
                SELECT gene.stable_id, gene_attrib.value, gene.description, xref.description, gene.biotype
                FROM gene
                LEFT JOIN xref ON gene.display_xref_id = xref.xref_id
                LEFT JOIN gene_attrib ON gene.gene_id = gene_attrib.gene_id
                WHERE (gene_attrib.attrib_type_id = 4) 
                AND (gene_attrib.value LIKE '%{searchword}%' OR 
                gene.description LIKE '%{searchword}%' OR 
                xref.description LIKE '%{searchword}%')
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
                                "value": "Gene_name",
                                "biotype": "Biotype"})
        # Changing description columns name by column index since they were returned with the same name ("description")
        df.columns.values[2] = "Ensembl_description"
        df.columns.values[3] = "Ext_ref_description"
        
    # For other species, the gene name will not be fetched      
    else: 
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
    
    # Add URL to gene summary on Ensembl
    df["URL"] = "https://uswest.ensembl.org/" + "_".join(db.split("_")[:2]) + "/Gene/Summary?g=" + df["Ensembl_ID"]
    
    # Remove any duplicate search results from the master data frame and reset the index
    df = df.drop_duplicates().reset_index(drop=True)
    
    print(f"Query time: {round(time.time() - start_time, 2)} seconds")
    print(f"Genes fetched: {len(df)}")
    
    if save == True:
        df.to_csv("gget_results.csv", index=False)
    
    return df

def ref(species, which="all", release=None, save=False, FTP=False):
    """
    Function to fetch GTF and FASTA (cDNA and DNA) files from the Ensemble FTP site.
    
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
    
    - save
    If "save=True", the json containing all results is saved in the current directory. Only works if "returnval='json'".
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
    
    # If single which passed as string, convert to list
    if type(which) == str:
        which = [which]
    
    ## Return results
    if FTP == False:
        # Single entry for which
        if len(which) == 1:
            if which == ["all"]:
                ref_dict = {
                    species: {
                        "transcriptome": {
                            "ftp":cdna_url,
                            "ensembl_release": int(ENS_rel),
                            "release_date": cdna_date.split(" ")[0],
                            "release_time": cdna_date.split(" ")[1],
                            "bytes": cdna_size
                        },
                        "genome": {
                            "ftp":dna_url,
                            "ensembl_release": int(ENS_rel),
                            "release_date": dna_date.split(" ")[0],
                            "release_time": dna_date.split(" ")[1],
                            "bytes": dna_size
                        },
                        "annotation": {
                            "ftp":gtf_url,
                            "ensembl_release": int(ENS_rel),
                            "release_date": gtf_date.split(" ")[0],
                            "release_time": gtf_date.split(" ")[1],
                            "bytes": gtf_size
                        }
                    }
                }

            elif which == ["gtf"]:
                ref_dict = {
                    species: {
                        "annotation": {
                        "ftp":gtf_url,
                        "ensembl_release": int(ENS_rel),
                        "release_date": gtf_date.split(" ")[0],
                        "release_time": gtf_date.split(" ")[1],
                        "bytes": gtf_size
                        }
                    }
                }

            elif which == ["cdna"]:
                ref_dict = {
                    species: {
                        "transcriptome": {
                            "ftp":cdna_url,
                            "ensembl_release": int(ENS_rel),
                            "release_date": cdna_date.split(" ")[0],
                            "release_time": cdna_date.split(" ")[1],
                            "bytes": cdna_size
                        },
                    }
                }

            elif which == ["dna"]:
                ref_dict = {
                    species: {
                        "genome": {
                            "ftp":dna_url,
                            "ensembl_release": int(ENS_rel),
                            "release_date": dna_date.split(" ")[0],
                            "release_time": dna_date.split(" ")[1],
                            "bytes": dna_size
                        }
                    }
                }

            else:
                raise ValueError("Parameter 'which' must be 'json', or any one or a combination of the following: 'gtf', 'cdna', 'dna'.")

            if save == True:
                import json
                with open('ref.json', 'w', encoding='utf-8') as f:
                    json.dump(ref_dict, f, ensure_ascii=False, indent=4)

            return ref_dict

        # Return multiple results
        else:
            results = []
            for return_val in which:
                if return_val == "all":
                    raise ValueError("Parameter 'which' must be 'all', or one or a combination of the following: 'gtf', 'cdna', 'dna'.")
                elif return_val == "gtf":

                    species: {
                        "annotation": {
                        "ftp":gtf_url,
                        "ensembl_release": int(ENS_rel),
                        "release_date": gtf_date.split(" ")[0],
                        "release_time": gtf_date.split(" ")[1],
                        "bytes": gtf_size
                        }
                    }
                }
                    results.append(gtf_url)
                elif return_val == "cdna":
                    results.append(cdna_url)
                elif return_val == "dna":
                    results.append(dna_url)
                else:
                    raise ValueError("Parameter 'which' must be 'all', or any one or a combination of the following: 'gtf', 'cdna', 'dna'.")

            print(f"Fetching from Ensembl release: {ENS_rel}")
            return results
        
    # Return only the URLs as a list (instead of a json) if FTP==True    
    if FTP == True:
        # Single entry for which
        if len(which) == 1:
            if which == ["all"]:
                ref_dict = {
                    species: {
                        "transcriptome": {
                            "ftp":cdna_url,
                            "ensembl_release": int(ENS_rel),
                            "release_date": cdna_date.split(" ")[0],
                            "release_time": cdna_date.split(" ")[1],
                            "bytes": cdna_size
                        },
                        "genome": {
                            "ftp":dna_url,
                            "ensembl_release": int(ENS_rel),
                            "release_date": dna_date.split(" ")[0],
                            "release_time": dna_date.split(" ")[1],
                            "bytes": dna_size
                        },
                        "annotation": {
                            "ftp":gtf_url,
                            "ensembl_release": int(ENS_rel),
                            "release_date": gtf_date.split(" ")[0],
                            "release_time": gtf_date.split(" ")[1],
                            "bytes": gtf_size
                        }
                    }
                }

            elif which == ["gtf"]:
                ref_dict = {
                    species: {
                        "annotation": {
                        "ftp":gtf_url,
                        "ensembl_release": int(ENS_rel),
                        "release_date": gtf_date.split(" ")[0],
                        "release_time": gtf_date.split(" ")[1],
                        "bytes": gtf_size
                        }
                    }
                }

            elif which == ["cdna"]:
                ref_dict = {
                    species: {
                        "transcriptome": {
                            "ftp":cdna_url,
                            "ensembl_release": int(ENS_rel),
                            "release_date": cdna_date.split(" ")[0],
                            "release_time": cdna_date.split(" ")[1],
                            "bytes": cdna_size
                        },
                    }
                }

            elif which == ["dna"]:
                ref_dict = {
                    species: {
                        "genome": {
                            "ftp":dna_url,
                            "ensembl_release": int(ENS_rel),
                            "release_date": dna_date.split(" ")[0],
                            "release_time": dna_date.split(" ")[1],
                            "bytes": dna_size
                        }
                    }
                }

            else:
                raise ValueError("Parameter 'which' must be 'json', or any one or a combination of the following: 'gtf', 'cdna', 'dna'.")

            if save == True:
                import json
                with open('ref.json', 'w', encoding='utf-8') as f:
                    json.dump(ref_dict, f, ensure_ascii=False, indent=4)

            return ref_dict

        # Return multiple results
        else:
            results = []
            for return_val in which:
                if return_val == "all":
                    raise ValueError("Parameter 'which' must be 'all', or one or a combination of the following: 'gtf', 'cdna', 'dna'.")
                elif return_val == "gtf":

                    species: {
                        "annotation": {
                        "ftp":gtf_url,
                        "ensembl_release": int(ENS_rel),
                        "release_date": gtf_date.split(" ")[0],
                        "release_time": gtf_date.split(" ")[1],
                        "bytes": gtf_size
                        }
                    }
                }
                    results.append(gtf_url)
                elif return_val == "cdna":
                    results.append(cdna_url)
                elif return_val == "dna":
                    results.append(dna_url)
                else:
                    raise ValueError("Parameter 'which' must be 'all', or any one or a combination of the following: 'gtf', 'cdna', 'dna'.")

            print(f"Fetching from Ensembl release: {ENS_rel}")
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

    # gget search subparser
    parser_gget = parent_subparsers.add_parser("search",
                                               parents=[parent],
                                               description="Query Ensembl for genes based on species and free form search terms.", 
                                               add_help=False)
    # Search arguments
    parser_gget.add_argument(
        "-sw", "--searchwords", 
        nargs="*",     # 0 or more values expected => creates a list
        type=str, 
        required=True, 
        metavar="",    # Cleans up help message
        help="One or more free form searchwords for the query (if more than one: use space between searchwords), e.g. gaba nmda."
    )
    parser_gget.add_argument(
        "-s", "--species",  
        required=True, 
        metavar="",
        help="Species to be queried, e.g. homo_sapiens or human."
    )
    parser_gget.add_argument(
        "-l", "--limit", 
        type=int, 
        metavar="",
        help="Limits the number of results, e.g. 10 (default: None)."
    )
    parser_gget.add_argument(
        "-o", "--out",
        type=str,
        help=(
            "Path to the csv file the results will be saved in, e.g. '-o path/to/directory/results.csv'." 
            "Default: None (just prints results)."
        )
    )


    # gget ref subparser
    parser_ref = parent_subparsers.add_parser("ref",
                                              parents=[parent],
                                              description="Fetch GTF and/or FASTA (cDNA and/or DNA) files for a specific species from the Ensemble FTP site.",
                                              add_help=False)
    # ref arguments
    parser_ref.add_argument(
        "-s", "--species", 
        required=True,
        type=str,
        metavar="", 
        help="Species for which the FTPs will be fetched, e.g. homo_sapiens."
    )
    parser_ref.add_argument(
        "-w", "--which", 
        default="all", 
        nargs="*", 
        type=str,
        metavar="",
        help=("Defines which results to return." 
              "Possible entries are: 
              "'all' - Returns GTF, cDNA, and DNA links and associated info (default)." 
              "Or one or a combination of the following:"  
              "'gtf' - Returns the GTF FTP link and associated info." 
              "'cdna' - Returns the cDNA FTP link and associated info."
              "'dna' - Returns the DNA FTP link and associated info."
             )
    parser_ref.add_argument(
        "-r", "--release",  
        type=int, 
        metavar="",
        help="Ensemble release the FTPs will be fetched from, e.g. 104 (default: latest Ensembl release).")
    parser_ref.add_argument(
        "-ftp", "--ftp",  
        default=False, 
        action='store_true'
        help="If True: return only the FTP link instead of a json.")
    parser_ref.add_argument(
        "-o", "--out",
        type=str,
        help=(
            "Path to the json file the results will be saved in, e.g. '-o path/to/directory/results.json'." 
            "Default: None (just prints results)."
        )
    )
    
    # Show help when no arguments are given
    if len(sys.argv) == 1:
        parent_parser.print_help(sys.stderr)
        sys.exit(1)

    args = parent_parser.parse_args()

    ## Define return values
    # debug return
    if args.debug:
        print("debug: " + str(args))
        
    # search return
    if args.command == "search":
        gget_results = gget(args.searchwords, args.species, args.limit)
        
        # Save in specified directory if -o specified
        if args.out:
            os.makedirs(args.out, exist_ok=True)
            gget_results.to_csv(args.out, index=False)
            print(f"Results saved as {args.out}.")
        
        # Print results if no directory specified
        else:
            print(gget_results)
            print("To save these results, use flag '-o' in the format: '-o path/to/directory/results.csv'.")
            
    # ref return
    if args.command == "ref":
        ref_results = ref(args.species, args.which, args.release, args.ftp)
        
        if args.ftp:
        
        
        
        # Print or save json file
        else:
            import json
            # Print if '-o' flag not True
            if args.out == False:
                print(json.dumps(ref_results, ensure_ascii=False, indent=4))
                print("To save these results in the current working directory, add the flag '-o'.")
            # Save in current working directory if args.out=True
            if args.out == True:
                with open('ref.json', 'w', encoding='utf-8') as f:
                    json.dump(ref_results, f, ensure_ascii=False, indent=4)
                print("Results saved in current working directory.")
        # If not json, return space sparated list of requested urls
        else:
            if args.out == False:
                print(" ".join(ref_results))
#             if args.out == True:
                
if __name__ == '__main__':
    main()
