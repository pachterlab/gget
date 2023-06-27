import requests
import pandas as pd
import numpy as np
import logging
from bs4 import BeautifulSoup
from io import StringIO

# Add and format time stamp in logging messages
logging.basicConfig(
    format="%(asctime)s %(levelname)s %(message)s",
    level=logging.INFO,
    datefmt="%c",
)
# Mute numexpr threads info
logging.getLogger("numexpr").setLevel(logging.WARNING)

# Call elm api to get elm id, start, stop and boolean values
# Returns tab separated values
def get_response_api(seq):
    url = "http://elm.eu.org/start_search/"
    # Build URL
    html = requests.get(url + seq)

    # Raise error if status code not "OK" Response
    if html.status_code != 200:
        raise RuntimeError(
            f"The ELM server returned error status code {html.status_code}. Please try again."
        )
    soup = BeautifulSoup(html.text, "html.parser")
    soup_string = str(soup)
    return soup_string

# Scrapes webpage for information about functional site class, description, pattern probability
# Return html tags in text format
def get_html(elm_id):
    url = "http://elm.eu.org/elms/"
    # Build URL
    resp = requests.get(url + elm_id)
    html = resp.text

    # Raise error if status code not "OK" Response
    if resp.status_code != 200:
        raise RuntimeError(
            f"The ELM server returned error status code {resp.status_code}. Please try again."
        )

    return html

# Searches through separated tab values soup for tags corresponding to field param
# Return text string inside html tags 
def get_additional_info(field, soup):
    if (field=="Interaction Domain:"):
        info = soup.find(id="interaction_domain").findNext('td').text
    else:
        info = soup.find(text=field).findNext('td').text
    return info

# Convert tab separated values to dataframe format
def tsv_to_df(tab_separated_values):
    df = pd.read_csv(StringIO(tab_separated_values), sep='\t', names=['elm_identifier', 'start', 'stop', 'is_annotated', 'is_phiblastmatch', 'is_filtered', 'phiblast', 'topodomfiler', 'taxonfilter', 'structure'])
    return df

def elm(
    sequence,
    uniprot=False,
):
    """
    Searches the Eukaryotic Linear Motif resource for Functional Sites in Proteins.
    Args:
     - sequence       amino acid sequence or Uniprot ID
                      (If more than one sequence in FASTA file, only the first will be submitted to BLAST.)

     - uniprot        If True, searches using Uniprot ID instead of amino acid sequence. Default: False

    Returns a data frame with the ELM results.

    
    """
    # Server rules:

    # Searches via Uniprot ID are limited to a maximum of 1 every 3 minutes.
    # If you exceed this limit, you will recieve a "429 Too many requests" error.
    # 
    # Searches of a single sequence is to a maximum of 1 per minute.
    # If you exceed this limit, you will recieve a "429 Too many requests" error.
    # 
    # Note: this does not always work for sequences longer than 2000 amino acids:
    # URLs may be truncated beyond this length.

    
    if uniprot:
        sequence = sequence + ".tsv"

    tab_separated_values = get_response_api(sequence)
    df = tsv_to_df(tab_separated_values)

    # for amino acid sequences, more information can be scraped from the webpage using elm ids returned
    if not uniprot:
        amino_acids = set("ARNDCQEGHILKMFPSTWYVBZXBJZ")

        # If sequence is not a valid amino sequence, raise error
        if not set(sequence) <= amino_acids:
            raise ValueError(
            f"'Input sequence is not a valid amino acid sequence. Please try again with either amino acid sequence or Uniprot ID"
            )

            
        column_names = ['Accession:',
        'Functional site class:',
        'Functional site description:',
        'ELM Description:',
        'ELMs with same func. site:',
        'Pattern:',
        'Pattern Probability:',
        'Present in taxons:',
        'Interaction Domain:', 
        ]

        # Creates new dataframe to store information from scraping
        df_2 = pd.DataFrame()
        # dropping the first row containing just the column name
        df = df.drop(index=0)
        # Grab elm identifiers column from dataframe
        elm_ids = df['elm_identifier'].values 

        # Add column of elm identifiers to new dataframe
        df_2['elm_identifier'] = elm_ids

        # Index dataframe using ELM id
        df_2 = df_2.set_index('elm_identifier')
        
        elm_id_index = 0
        # Loop through each elm identifier, get and parse html content 
        for elm_id in elm_ids:
            html = get_html(elm_id)
            soup = BeautifulSoup(html, "html.parser")
        
            for column in column_names:
                column_ignored_colon = column[:-1]
                try:
                    value = get_additional_info(column, soup)
                except AttributeError:
                    try:
                        column = "Present in taxon:" 
                        value = get_additional_info(column, soup)
                    except AttributeError:
                        logging.debug(f"No values for ELM ID: {elm_id} for {column_ignored_colon}")
                        continue
                # Clean up results and add to corresponding position in new dataframe
                value = value.strip().replace("\n", " ").replace("\t", " ")
                df_2.loc[elm_id, column_ignored_colon] = value

            # calculate position of motifs associated with each elm id compared to original sequence
            start_str = df.iloc[elm_id_index].loc['start']
            stop_str = df.iloc[elm_id_index]['stop']
            if (start_str != ""):
                start = int (start_str)
            if (stop_str != ""):
                stop = int (stop_str)
            df_2.loc[elm_id, 'Motif position in original sequence'] = sequence[start-1: stop-1]
            elm_id_index += 1

        # if a motif is found more than once in a UniProt ID or sequence, df_2 will have duplicate rows
        df_2 = df_2.drop_duplicates()

        # Merge two dataframes and sort by pattern probability 
        df_merge = df_2.merge(df, on='elm_identifier', how="right")
        df_merge = df_merge.sort_values(by='Pattern Probability', ascending=False)
        return df_merge

    # if search using Uniprot ID, return only dataframe from api call
    return df

