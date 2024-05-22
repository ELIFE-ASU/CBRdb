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
    return prefix + '{:05}'.format(id_number)


def get_mol(id_drug, save_dir, kegg_website="https://rest.kegg.jp/get/", request_sleep=0.6):
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


def clean_empty_files(target_dir):
    # Make a list of empty files
    empty_files = [file for file in os.listdir(target_dir) if os.path.getsize(os.path.join(target_dir, file)) == 0]
    # Remove the empty files
    for file in empty_files:
        os.remove(os.path.join(target_dir, file))
    n_rm = len(empty_files)
    print(f"Removed {n_rm} empty files")
    return n_rm


def get_kegg_all(target_dir="kegg_data", target="C"):
    max_drugs = 12897
    max_compounds = 22919

    if target == "D":
        max_idx = max_drugs
    elif target == "C":
        max_idx = max_compounds
    else:
        raise ValueError(f"Unknown target {target}")

    # Get the data
    get_kegg(target_dir + f"_{target}", prefix=target, max_idx=max_idx)


if __name__ == "__main__":
    # Get the data
    get_kegg_all(target="C")

    # Clean the data
    clean_empty_files("kegg_data_C")
