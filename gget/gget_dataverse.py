import os
import requests
from tqdm import tqdm
import urllib.request
import json
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


def process_local_json(filename):
    """Process a local JSON file.

    Args:
        filename (str): The path to the local JSON file.
    
    Returns:
        dict: The local JSON file information as a dictionary.
    """
    
    f = open(filename, 'r')
    data = json.load(f)

    return data


def process_remote_json(url, save=False):
    """Process a remote JSON file.

    Args:
        url (str): The URL of the remote JSON file.
        save (bool, optional): Whether to save the JSON file locally. Defaults to False.

    Returns:
        dict: The remote JSON file information as a dictionary.
    """
    response = urllib.request.urlopen(url)
    data = json.loads(response.read())

    # Save JSON
    if save:
        with open(save, "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False)

    return data


def dataverse(data, path=None, run_download=False, save_json=None):
    """process a json file including the dataverse datasets information and download the datasets

    Args:
        data (str or dict): list of datasets to download in JSON format. URL, path to a local file, or a python dictionary.
        path (str, optional): the path to save the datasets. Defaults to None.
        run_download (bool, optional): whether to download the datasets. Defaults to True.
        save_json (str): path to save JSON file. Defaults to None.
    """
    if "https" in data or "http" in data:
        data = process_remote_json(data, save=save_json)
    elif type(data) == str and '.json' in data:
        data = process_local_json(data)
    elif type(data) == dict:
        pass

    if "datasets" not in data:
        # TODO: Add more error handling
        raise ValueError("The json file must include proper 'datasets' key")

    if not path and not run_download:
        pass
    elif not path and run_download:
        raise ValueError("Please provide a path to save the datasets and set run_download=True")
    elif run_download:
        for entry in data['datasets']:
            download_wrapper(entry, path)
