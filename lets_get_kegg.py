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


def get_mol(id_drug, save_dir, kegg_website="https://rest.kegg.jp/get/", request_sleep=0.6):
    """
    This function downloads a molecule file from the KEGG database and saves it to a specified directory.
    If the file already exists in the directory, the download is skipped.

    Parameters:
    id_drug (str): The ID of the drug for which the molecule file is to be downloaded.
    save_dir (str): The directory where the molecule file is to be saved.
    kegg_website (str, optional): The base URL of the KEGG database. Defaults to "https://rest.kegg.jp/get/".
    request_sleep (float, optional): The time to sleep between requests in seconds. Defaults to 0.6.

    Returns:
    bool: True if the download was successful or the file already exists, False otherwise.
    """
    # Get the full file path
    full_file = os.path.join(save_dir, f"{id_drug}.mol")

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
            response = s.get(f"{kegg_website}{id_drug}/mol", timeout=10.0)
        except requests.exceptions.RequestException as e:
            # Some error in the connection
            print(f"Error in drug ID {id_drug}, connection exception {e}")
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
            print(f"Error in drug ID {id_drug}, response {response.status_code}")
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
        # Get the formatted drug id
        id_drug = format_mol_id(i, prefix=prefix)
        # Get the full path
        full_path = os.path.join(target_dir, id_drug)
        # Make subdirectory
        os.makedirs(full_path, exist_ok=True)
        print(f"Downloading drug {id_drug}")
        # Check if the drug is downloaded
        if not get_mol(id_drug, full_path):
            print(f"Download failed for drug {id_drug}")
        # Check if the maximum index is reached
        if max_idx is not None and i >= max_idx:
            print(f"Maximum index reached {max_idx}")
            break
    return None


def clean_empty_folders(target_dir):
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
    # Make a list of empty folders
    empty_folder = [folder for folder in os.listdir(target_dir) if
                    os.path.getsize(os.path.join(target_dir, folder)) == 0]
    # Remove the empty folder
    for folder in empty_folder:
        os.rmdir(os.path.join(target_dir, folder))
    n_rm = len(empty_folder)
    print(f"Removed {n_rm} empty folders")
    return n_rm


def get_total_n_reactions(database="reaction", kegg_website=r"https://rest.kegg.jp/list/", request_sleep=0.6):
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
    """
    This function downloads all molecule files from the KEGG database for a specified target and saves them to a specified directory.
    The target can be either "D" for drugs or "C" for compounds.
    The maximum index for the molecule ID to be downloaded is determined based on the target.
    The function calls the get_kegg function to perform the download.

    Parameters:
    target_dir (str, optional): The directory where the molecule files are to be saved. Defaults to "kegg_data".
    target (str, optional): The target for which the molecule files are to be downloaded. Can be either "D" for drugs or "C" for compounds. Defaults to "C".

    Raises:
    ValueError: If an unknown target is specified.

    Returns:
    None
    """

    if target == "D":
        _, max_idx = get_total_n_reactions(database="drug")
    elif target == "C":
        _, max_idx = get_total_n_reactions(database="compound")
    else:
        raise ValueError(f"Unknown target {target}")

    # Get the data
    get_kegg(target_dir + f"_{target}", prefix=target, max_idx=max_idx)


if __name__ == "__main__":
    # Get the data
    get_kegg_all(target="C")

    # Clean the data
    clean_empty_folders("kegg_data_C")
