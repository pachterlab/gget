import numpy as np
import json
import logging
# Add and format time stamp in logging messages
logging.basicConfig(
    format="%(asctime)s %(levelname)s %(message)s", 
    level=logging.INFO,
    datefmt="%d %b %Y %H:%M:%S",
)
# Custom functions
from .utils import rest_query
# Constants
from .constants import ENSEMBL_REST_API

## gget info
def info(
    ens_ids, 
    expand=False, 
    homology=False, 
    xref=False, 
    save=False,
    verbose=True
):
    """
    Fetch gene and transcript metadata using Ensembl IDs.

    Args:
    - ens_ids   One or more Ensembl IDs to look up (string or list of strings).
    - expand    Expand returned information (only for genes and transcripts) (default: False). 
                For genes, this adds isoform information. 
                For transcripts, this adds translation and exon information.
    - homology  If True, also returns homology information of ID (default: False).
    - xref      If True, also returns information from external references (default: False).
    - save      If True, saves json with query results in current working directory (default: False).
    - verbose   If True, prints progress information (default: True).

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
        except RuntimeError:
            # Try submitting query without expand (expand does not work for exons and translation IDs)
            try:
                query = "lookup/id/" + ensembl_ID + "?"
                df_temp = rest_query(server, query, content_type)
            # Raise error if this also did not work
            except RuntimeError:
                if verbose == True:
                    logging.error(
                        f"Ensembl ID {ensembl_ID} not found. "
                        "Please double-check spelling/arguments and try again."
                    )
                return
            
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
                if verbose == True:
                    logging.warning(f"No homology information found for {ensembl_ID}.")

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
                if verbose == True:
                    logging.warning(f"No external reference information found for {ensembl_ID}.")
    
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
