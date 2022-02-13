import pandas as pd
# pure-Python MySQL client library for ENSEMBL webscraping
import pymysql

def gget(searchwords, species, limit=None):
    """
    Function to query Ensembl for genes based on species and free form search terms. 
    
    The variable "searchwords" is a list strings containing the free form search terms 
    (e.g.searchwords = ["GABA", "gamma-aminobutyric acid"]).
    The search is not case-sensitive.
    
    "Limit" limits the number of search results to the top {limit} genes found.
    
    Possible entries for species are:
    "homo_sapiens" (or "human")
    "mus_musculus" (or "mouse")
    "taeniopygia_guttata" (or "zebra finch")
    "caenorhabditis_elegans" (or "roundworm")
    
    If you would like to access results from a different species or restrain your results to a certain mouse strain, 
    you can instead enter the core database as the "species" variable (e.g. species = "rattus_norvegicus_core_105_72").
    You can find all availabale species databases here: http://ftp.ensembl.org/pub/release-105/mysql/
    
    Returns a data frame which contains the cleaned up, combined results for all search words and
    a dictionary containing the raw results for each searchword (the searchwords are the dictionary keys and the values 
    for each key are the search results).
    """
    
    if species == "caenorhabditis_elegans" or species == "roundworm":
        db = "caenorhabditis_elegans_core_105_269"   
    elif species == "homo_sapiens" or species == "human":
        db = "homo_sapiens_core_105_38"
    elif species == "mus_musculus" or species == "mouse":
        db = "mus_musculus_core_105_39"
    elif species == "taeniopygia_guttata" or species == "zebra finch":
        db = "taeniopygia_guttata_core_105_12"
    else: 
        db = species
        
    print(f"Results fetched from database: {db}")
    
    database = pymysql.connect(
        host="ensembldb.ensembl.org",
        user="anonymous",
        password="",
    )

    results = {}

    for searchword in searchwords:
        cursor = database.cursor()
        # If limit is specified, fetch only the first {limit} genes for which the searchword appears in the description
        if limit != None:
            cursor.execute(
                f"select * from gene where description like '%{searchword}%' limit {limit};"
            )
        # Else, fetch all genes for which the searchword appears in the description
        else:
            cursor.execute(f"select * from gene where description like '%{searchword}%';")
        # Fetch the search results and save as "res"
        res = cursor.fetchall()

        # Add search results for this searchword to the dictionary
        results[searchword] = res

    # Add the results for all search words to a data frame
    for i, res in enumerate(results):
        # Create dataframe from search results
        df_temp = pd.DataFrame(
            results[res],
            columns=[
                "0", "Biotype", "2", "3", "4", "5", "6", "7", "8", "Description", "10", "11", "Gene ID", "13", "14", "15",
            ],
        )
        # Drop columns containing irrelevant data
        df_temp = df_temp.drop(
            ["0", "2", "3", "4", "5", "6", "7", "8", "10", "11", "13", "14", "15"], axis=1
        )

        # In the first iteration, make the search results equal to the master data frame
        if i == 0:
            df = df_temp.copy()
        # Add new search results to master data frame
        else:
            df = pd.concat([df, df_temp])

    # Remove any duplicate search results from the master data frame and reset the index
    df = df.drop_duplicates()
    df = df.reset_index(drop=True)

    # Close cursor and database
    cursor.close()
    database.close()
    
    return df, results