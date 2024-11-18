import os
import time

import requests
from requests import Session
from requests.adapters import HTTPAdapter
from urllib3.util import Retry


def load_bad_entries(bad_file, target_str="molless"):
    with open(bad_file, 'r') as file:
        return [line.split(',')[0].strip() for line in file if target_str in line]


def prepare_session():
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
    return s


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


def clean_empty_folders(target_dir, size=False):
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
            print(f"Error removing folder {folder}: {e}", flush=True)
    print(f"Removed {n_rm} empty folders", flush=True)
    return n_rm


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


def get_total_n(session, database="reaction", kegg_website=r"https://rest.kegg.jp/list/", request_sleep=0.2):
    # Get the data
    try:
        # Get the response
        response = session.get(f"{kegg_website}{database}", timeout=10.0)
        # Limit the number of requests
        time.sleep(request_sleep)
    except requests.exceptions.RequestException as e:
        # Some error in the connection
        print(f"Error connection exception {e}", flush=True)
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


def check_kegg_valid_id(id, session, kegg_website="https://rest.kegg.jp/get/", request_sleep=0.2):
    # Get the data
    try:
        response = session.get(f"{kegg_website}{id}", timeout=10.0)
        # Limit the number of requests
        time.sleep(request_sleep)
    except requests.exceptions.RequestException as e:
        # Some error in the connection
        print(f"Error in ID {id}, connection exception {e}", flush=True)
        return False
    # Check if the response is ok
    if response.ok:
        return True
    else:
        # Some error in the response
        print(f"Error in ID {id}, not a valid ID", flush=True)
        return False


def get_data(id, save_dir, session, kegg_website="https://rest.kegg.jp/get/", request_sleep=0.4, full=False):
    # Get the data type from the id
    data_type = id[0]
    # Get the full file path
    if data_type == "R" in save_dir or full:
        full_file = os.path.join(save_dir, f"{id}.data")
    else:
        full_file = os.path.join(save_dir, f"{id}.mol")

    # Check if the file exists
    if not os.path.exists(full_file):
        # Get the data
        try:
            # Get the response
            if data_type == "R" in save_dir or full:
                response = session.get(f"{kegg_website}{id}", timeout=10.0)
            else:
                response = session.get(f"{kegg_website}{id}/mol", timeout=10.0)
            # Limit the number of requests
            time.sleep(request_sleep)
        except requests.exceptions.RequestException as e:
            # Some error in the connection
            print(f"Error in ID {id}, connection exception {e}", flush=True)
            return False
        # Check if the response is ok
        if response.ok:
            # Strip the last section of the file
            res = response.text.split("> <ENTRY>")[0]
            # Make subdirectory
            os.makedirs(save_dir, exist_ok=True)
            # Save the file
            with open(full_file, "w") as f:
                f.write(res)
        else:
            # Some error in the response
            print(f"Error in ID {id}, response {response.status_code}", flush=True)
            return False
    else:
        # Skip the download as the file already exists
        print(f"{id} already exists", flush=True)
    return True


def get_kegg(target_dir, session,
             prefix="D",
             max_idx=12897):
    bad_file = os.path.abspath(f"data/{prefix.replace("_full", "")}_IDs_bad.dat")
    # Check if the prefix is to download the full data
    if "_full" in prefix:
        full = True
    else:
        full = False
    print(f"Downloading {prefix} data to {target_dir}", flush=True)
    # Clean up the prefix
    prefix = prefix.strip("_full").upper()
    # Check path for saving
    os.makedirs(target_dir, exist_ok=True)

    if os.path.exists(bad_file):
        molless_ids = load_bad_entries(bad_file, target_str="molless")
        print(f"Loaded {len(molless_ids)} Mol-less IDs", flush=True)
        invalid_ids = load_bad_entries(bad_file, target_str="invalid")
        print(f"Loaded {len(invalid_ids)} invalid IDs", flush=True)
    else:
        molless_ids = []
        invalid_ids = []

    # Start the loop
    i = 0
    while True:
        i += 1
        # Get the formatted id
        id = format_mol_id(i, prefix=prefix)
        # Get the full path
        full_path = os.path.join(target_dir, id)
        print(f"Downloading {id}", flush=True)
        # Check if the ID is in the invalid list
        if id in invalid_ids:
            print(f"Skipping {id} as it is invalid", flush=True)
            continue
        # Check if the ID is in the Mol-less list
        if id in molless_ids:
            print(f"Skipping {id} as it is Mol-less", flush=True)
            continue

        # Check if the id is downloaded
        if not get_data(id, full_path, session, full=full):
            # check if the id is valid
            if check_kegg_valid_id(id, session):
                # The ID is valid but the data is not downloaded, write the Mol-less id to file
                print(f"{id} is Mol-less, written to file", flush=True)
                molless_ids.append(id)
            else:
                # The ID is invalid, write the invalid id to file
                print(f"{id} is invalid, written to file", flush=True)
                invalid_ids.append(id)

        # Check if the maximum index is reached and break
        if max_idx is not None and i >= max_idx:
            print(f"Maximum index reached {max_idx}", flush=True)
            break

    # Write a log file
    with open(bad_file, "w") as f:
        f.write("# Bad IDs, reason\n")
        for id in molless_ids:
            f.write(f"{id}, molless\n")
        for id in invalid_ids:
            f.write(f"{id}, invalid\n")

    return None


def get_kegg_all(target_dir="kegg_data", target="C"):
    # make the session
    session = prepare_session()
    # Check if the target is a valid target and get the maximum index
    if target == "D" or target == "D_full":
        _, max_idx = get_total_n(session, database="drug")
    elif target == "C" or target == "C_full":
        _, max_idx = get_total_n(session, database="compound")
    elif target == "R":
        _, max_idx = get_total_n(session, database="reaction")
    else:
        raise ValueError(f"Unknown target {target}")

    print(f"Total number of entries {max_idx}", flush=True)

    # Get the data
    get_kegg(os.path.abspath(target_dir + f"_{target}"), session, prefix=target, max_idx=max_idx)


def main(target="R", target_dir=r"..\data\kegg_data"):
    target_dir = os.path.abspath(target_dir)
    # Clean the data
    clean_empty_folders(f"{target_dir}_{target}")
    # Get the data
    get_kegg_all(target_dir=target_dir, target=target)
    # Clean the data
    clean_empty_folders(f"{target_dir}_{target}")


if __name__ == "__main__":
    print("Program started", flush=True)
    main(target="C")
    main(target="R")
    main(target="C_full")
    print("Program finished", flush=True)
