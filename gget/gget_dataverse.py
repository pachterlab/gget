import os
import requests
from tqdm import tqdm
import pandas as pd
from .utils import print_sys
from .constants import DATAVERSE_GET_URL

def dataverse_downloader(url, path, file_name):
    """dataverse download helper with progress bar

    Args:
        url (str): the url of the dataset to download
        path (str): the path to save the dataset locally
        file_name (str): the name of the file to save locally
    """
    save_path = os.path.join(path, file_name)
    response = requests.get(url, stream=True)
    total_size_in_bytes = int(response.headers.get("content-length", 0))
    block_size = 1024
    progress_bar = tqdm(total=total_size_in_bytes, unit="iB", unit_scale=True)
    with open(save_path, "wb") as file:
        for data in response.iter_content(block_size):
            progress_bar.update(len(data))
            file.write(data)
    progress_bar.close()


def download_wrapper(entry, path, return_type=None):
    """wrapper for downloading a dataset given the name and path, for csv,pkl,tsv or similar files

    Args:
        entry (dict): the entry of the dataset to download. Must include 'id', 'name', 'type' keys
        path (str): the path to save the dataset locally
        return_type (str, optional): the return type. Defaults to None. Can be "url", "filename", or ["url", "filename"]

    Returns:
        str: the exact dataset query name
    """
    url = DATAVERSE_GET_URL + str(entry['id'])

    if not os.path.exists(path):
        os.mkdir(path)

    filename = f"{entry['name']}.{entry['type']}"

    if os.path.exists(os.path.join(path, filename)):
        print_sys(f"Found local copy for {entry['id']} datafile as {filename} ...")
        os.path.join(path, filename)
    else:
        print_sys(f"Downloading {entry['id']} datafile as {filename} ...")
        dataverse_downloader(url, path, filename)
    
    if return_type == "url":
        return url
    elif return_type == "filename":
        return filename
    elif return_type == ["url", "filename"]:
        return url, filename


def dataverse(df, path, sep=","):
    """download datasets from dataverse for a given dataframe
    Input dataframe must have 'name', 'id', 'type' columns.
    - 'name' is the dataset name for single file
    - 'id' is the unique identifier for the file
    - 'type' is the file type (e.g. csv, tsv, pkl)

    Args:
        df (pd.DataFrame or str): the dataframe or path to the csv/tsv file
        path (str): the path to save the dataset locally
    """
    if type(df) == str:
        if os.path.exists(df):
            df = pd.read_csv(df, sep=sep)
        else:
            raise FileNotFoundError(f"File {df} not found")
    elif type(df) == pd.DataFrame:
        pass
    else:
        raise ValueError("Input must be a pandas dataframe or a path to a csv / tsv file")
    
    print_sys(f"Searching for {len(df)} datafiles in dataverse ...")

    # run the download wrapper for each entry in the dataframe
    for _, entry in df.iterrows():
        download_wrapper(entry, path)
    
    print_sys(f"Download completed, saved to `{path}`.")
