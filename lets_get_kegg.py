import os
import time

import requests
from requests import Session
from requests.adapters import HTTPAdapter
from urllib3.util import Retry


def file_list(mypath=None):
    """
    This function generates a list of all files in a specified directory.
    If no directory is specified, it defaults to the current working directory.

    Parameters:
    mypath (str, optional): The path to the directory. Defaults to None, which means the current working directory.

    Returns:
    list: A list of all files in the specified directory.
    """
    mypath = mypath or os.getcwd()  # If no path is provided, use the current working directory
    return [f for f in os.listdir(mypath) if
            os.path.isfile(os.path.join(mypath, f))]  # Return a list of all files in the directory


def format_mol_id(id_number, prefix="D"):
    """
    This function formats a molecule ID by adding a prefix and ensuring the ID number is five digits long.

    Parameters:
    id_number (int): The ID number of the molecule.
    prefix (str, optional): The prefix to be added to the ID number. Defaults to "D".

    Returns:
    str: The formatted molecule ID.
    """
    return prefix + '{:05}'.format(id_number)


def get_data(id, save_dir, kegg_website="https://rest.kegg.jp/get/", request_sleep=0.6):
    # Get the data type from the id
    data_type = id[0]
    # Get the full file path
    if data_type == "R":
        full_file = os.path.join(save_dir, f"{id}.data")
    else:
        full_file = os.path.join(save_dir, f"{id}.mol")

    # Make the session
    s = Session()
    # Add retries
    retries = Retry(
        total=5,
        backoff_factor=0.1,
        status_forcelist=[502, 503, 504],
        allowed_methods={'POST'},
    )
    # Mount the session
    s.mount('https://', HTTPAdapter(max_retries=retries))

    # Check if the file exists
    if not os.path.exists(full_file):
        # Limit the number of requests
        time.sleep(request_sleep)
        # Get the data
        try:
            # Get the response
            if data_type == "R":
                response = s.get(f"{kegg_website}{id}", timeout=10.0)
            else:
                response = s.get(f"{kegg_website}{id}/mol", timeout=10.0)
        except requests.exceptions.RequestException as e:
            # Some error in the connection
            print(f"Error in ID {id}, connection exception {e}")
            return False
        # Check if the response is ok
        if response.ok:
            # Strip the last section of the file
            res = response.text.split("> <ENTRY>")[0]
            # Save the file
            with open(full_file, "w") as f:
                f.write(res)
        else:
            # Some error in the response
            print(f"Error in ID {id}, response {response.status_code}")
            return False
    else:
        # Skip the download as the file already exists
        print(f"File {full_file} already exists")
    return True


def get_kegg(target_dir, prefix="D", max_idx=12897):
    """
    This function downloads molecule files from the KEGG database and saves them to a specified directory.
    It creates a subdirectory for each molecule and downloads the molecule file into that subdirectory.
    The download is done in a loop, starting from molecule ID 1 and going up to a maximum index.
    If a molecule file cannot be downloaded, the function prints a message and continues with the next molecule.
    If the maximum index is reached, the function prints a message and breaks the loop.

    Parameters:
    target_dir (str): The directory where the molecule files are to be saved.
    prefix (str, optional): The prefix to be added to the molecule ID. Defaults to "D".
    max_idx (int, optional): The maximum index of the molecule ID to be downloaded. Defaults to 12897.

    Returns:
    None
    """
    # Check path for saving
    os.makedirs(target_dir, exist_ok=True)
    # Start the loop
    i = 0
    while True:
        i += 1
        # Get the formatted id
        id = format_mol_id(i, prefix=prefix)
        # Get the full path
        full_path = os.path.join(target_dir, id)
        # Make subdirectory
        os.makedirs(full_path, exist_ok=True)
        print(f"Downloading {id}")
        # Check if the drug is downloaded
        if not get_data(id, full_path):
            print(f"Download failed for {id}")
        # Check if the maximum index is reached
        if max_idx is not None and i >= max_idx:
            print(f"Maximum index reached {max_idx}")
            break
    return None


def list_empty_folders(target_dir):
    """
    List all empty folders within a specified directory.

    Parameters:
    target_dir (str): The path to the directory to search in.

    Returns:
    list: A list of paths to all empty folders within the specified directory.
    """
    empty_folders = []
    for entry in os.listdir(target_dir):
        entry_path = os.path.join(target_dir, entry)
        if os.path.isdir(entry_path) and not os.listdir(entry_path):
            empty_folders.append(entry_path)
    return empty_folders

def clean_empty_folders(target_dir,size=False):
    """
    This function removes all empty folders from a specified directory.
    It first generates a list of all empty folders in the directory,
    then removes each folder in the list.
    It prints the number of folders removed and returns this number.

    Parameters:
    target_dir (str): The directory from which empty folders are to be removed.

    Returns:
    int: The number of folders removed.
    """
    # Make a list of empty folders using the size of the folder
    e1 = [folder for folder in os.listdir(target_dir) if
                    os.path.getsize(os.path.join(target_dir, folder)) == 0]
    # Make a list of empty folders using the os.listdir function
    e2 = list_empty_folders(target_dir)
    # Combine the two lists
    if size:
        empty_folder = set(e1 + e2)
    else:
        empty_folder = set(e2)
    # Remove the empty folder
    n_rm = 0
    for folder in empty_folder:
        try:
            os.rmdir(os.path.join(target_dir, folder))
            n_rm += 1
        except OSError as e:
            print(f"Error removing folder {folder}: {e}")
    print(f"Removed {n_rm} empty folders")
    return n_rm


def get_total_n(database="reaction", kegg_website=r"https://rest.kegg.jp/list/", request_sleep=0.6):
    # Make the session
    s = Session()
    # Add retries
    retries = Retry(
        total=5,
        backoff_factor=0.1,
        status_forcelist=[502, 503, 504],
        allowed_methods={'POST'},
    )
    # Mount the session
    s.mount('https://', HTTPAdapter(max_retries=retries))

    # Limit the number of requests
    time.sleep(request_sleep)
    # Get the data
    try:
        # Get the response
        response = s.get(f"{kegg_website}{database}", timeout=10.0)
    except requests.exceptions.RequestException as e:
        # Some error in the connection
        print(f"Error connection exception {e}")
        exit()
    # Check if the response is ok
    if response.ok:
        # Strip the last section of the file
        res = response.text.split("> <ENTRY>")[0]
        # Split the response into lines
        lines = res.split("\n")
        # Count the number of lines
        n = len(lines)
        # find the index of the last line
        idx = int(lines[-2].split()[0][1:])
        return n, idx
    else:
        # Some error in the response
        print(f"Error response {response.status_code}")
        exit()


def get_kegg_all(target_dir="kegg_data", target="C"):
    if target == "D":
        _, max_idx = get_total_n(database="drug")
    elif target == "C":
        _, max_idx = get_total_n(database="compound")
    elif target == "R":
        _, max_idx = get_total_n(database="reaction")
    else:
        raise ValueError(f"Unknown target {target}")

    # Get the data
    get_kegg(target_dir + f"_{target}", prefix=target, max_idx=max_idx)


if __name__ == "__main__":
    print("Program started", flush=True)
    target = "R"
    target_dir = r"C:\Users\louie\skunkworks\data\kegg_data"

    # Get the data
    #get_kegg_all(target_dir=target_dir, target=target)
    # Clean the data
    clean_empty_folders(f"{target_dir}_{target}")
    print("Program finished", flush=True)
