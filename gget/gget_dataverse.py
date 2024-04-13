import os
import requests
from tqdm import tqdm
from utils import print_sys


def dataverse_download(url, path, name, types):
    """dataverse download helper with progress bar

    Args:
        url (str): the url of the dataset
        path (str): the path to save the dataset
        name (str): the dataset name
        types (dict): a dictionary mapping from the dataset name to the file format
    """
    save_path = os.path.join(path, f"{name}.{types[name]}")
    response = requests.get(url, stream=True)
    total_size_in_bytes = int(response.headers.get("content-length", 0))
    block_size = 1024
    progress_bar = tqdm(total=total_size_in_bytes, unit="iB", unit_scale=True)
    with open(save_path, "wb") as file:
        for data in response.iter_content(block_size):
            progress_bar.update(len(data))
            file.write(data)
    progress_bar.close()


def download_wrapper(name, path, return_type=None):
    """wrapper for downloading a dataset given the name and path, for csv,pkl,tsv or similar files

    Args:
        name (str): the rough dataset query name
        path (str): the path to save the dataset
        return_type (str, optional): the return type. Defaults to None. Can be "url", "name", or ["url", "name"]

    Returns:
        str: the exact dataset query name
    """
    server_path = "https://dataverse.harvard.edu/api/access/datafile/"

    url = server_path + str(name2id[name])

    if not os.path.exists(path):
        os.mkdir(path)

    file_name = f"{name}.{name2type[name]}"

    if os.path.exists(os.path.join(path, file_name)):
        print_sys("Found local copy...")
        os.path.join(path, file_name)
    else:
        print_sys("Downloading...")
        dataverse_download(url, path, name, name2type)
    
    if return_type == "url":
        return url
    elif return_type == "name":
        return file_name
    elif return_type == ["url", "name"]:
        return url, file_name
