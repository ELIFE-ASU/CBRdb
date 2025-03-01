import os
import time
import pandas as pd

import requests

from .tools_files import clean_empty_folders, make_custom_id, file_list_all
from .tools_requests import prepare_session


def load_bad_entries(bad_file, target_str="molless"):
    """
    Loads entries from a file that contain a specific target string.

    Parameters:
    bad_file (str): The path to the file containing the entries.
    target_str (str, optional): The target string to search for in each line of the file. Defaults to "molless".

    Returns:
    list: A list of entries from the file that contain the target string.
    """
    with open(bad_file, 'r') as file:
        return [line.split(',')[0].strip() for line in file if target_str in line]


def get_total_n(session,
                database="reaction",
                kegg_website="https://rest.kegg.jp/list/",
                request_sleep=0.2,
                timeout=10.0):
    """
    Retrieves the total number of entries and the index of the last entry from the KEGG database.

    Parameters:
    session (requests.Session): The session object to use for making requests.
    database (str): The name of the KEGG database to query. Defaults to "reaction".
    kegg_website (str): The base URL of the KEGG REST API. Defaults to "https://rest.kegg.jp/list/".
    request_sleep (float): The time to sleep between requests to avoid overloading the server. Defaults to 0.2 seconds.
    timeout (float): The timeout for the request in seconds. Defaults to 10.0 seconds.

    Returns:
    tuple: A tuple containing the total number of entries (n) and the index of the last entry (idx).

    Raises:
    ConnectionError: If there is an issue with the connection or the response is not OK.
    """
    try:
        response = session.get(f"{kegg_website}{database}", timeout=timeout)
        time.sleep(request_sleep)
    except requests.exceptions.RequestException as e:
        print(f"Error connection exception {e}", flush=True)
        raise ConnectionError

    if response.ok:
        lines = response.text.split("\n")
        n = len(lines)
        idx = int(lines[-2].split()[0][1:])
        valid_ids = [line.split()[0] for line in lines[:-1]]
        return n, idx, valid_ids
    else:
        print(f"Error response {response.status_code}")
        raise ConnectionError


def check_kegg_valid_id(id,
                        session,
                        kegg_website="https://rest.kegg.jp/get/",
                        request_sleep=0.2,
                        timeout=10.0):
    """
    Checks if a given KEGG ID is valid by making a request to the KEGG REST API.

    Parameters:
    id (str): The KEGG ID to be checked.
    session (requests.Session): The session object to use for making requests.
    kegg_website (str, optional): The base URL of the KEGG REST API. Defaults to "https://rest.kegg.jp/get/".
    request_sleep (float, optional): The time to sleep between requests to avoid overloading the server. Defaults to 0.2 seconds.
    timeout (float, optional): The timeout for the request in seconds. Defaults to 10.0 seconds.

    Returns:
    bool: True if the KEGG ID is valid, False otherwise.
    """
    try:
        response = session.get(f"{kegg_website}{id}", timeout=timeout)
        time.sleep(request_sleep)
    except requests.exceptions.RequestException as e:
        print(f"Error in ID {id}, connection exception {e}", flush=True)
        return False

    if response.ok:
        return True
    else:
        print(f"Error in ID {id}, not a valid ID", flush=True)
        return False


def get_data(id,
             save_dir,
             session,
             kegg_website="https://rest.kegg.jp/get/",
             request_sleep=0.4,
             timeout=10.0,
             full=False):
    """
    Retrieves data for a given KEGG ID and saves it to a specified directory.

    Parameters:
    id (str): The KEGG ID to retrieve data for.
    save_dir (str): The directory to save the retrieved data.
    session (requests.Session): The session object to use for making requests.
    kegg_website (str, optional): The base URL of the KEGG REST API. Defaults to "https://rest.kegg.jp/get/".
    request_sleep (float, optional): The time to sleep between requests to avoid overloading the server. Defaults to 0.4 seconds.
    timeout (float, optional): The timeout for the request in seconds. Defaults to 10.0 seconds.
    full (bool, optional): Whether to retrieve the full data or not. Defaults to False.

    Returns:
    bool: True if the data was successfully retrieved and saved, False otherwise.
    """
    data_type = id[0]
    full_file = os.path.join(save_dir, f"{id}.data" if data_type == "R" in save_dir or full else f"{id}.mol")

    if not os.path.exists(full_file):
        try:
            response = session.get(
                f"{kegg_website}{id}" if data_type == "R" in save_dir or full else f"{kegg_website}{id}/mol",
                timeout=timeout)
            time.sleep(request_sleep)
        except requests.exceptions.RequestException as e:
            print(f"Error in ID {id}, connection exception {e}", flush=True)
            return False

        if response.ok:
            os.makedirs(save_dir, exist_ok=True)
            with open(full_file, "w") as f:
                f.write(response.text.split("> <ENTRY>")[0])
        else:
            print(f"Error in ID {id}, response {response.status_code}", flush=True)
            return False
    else:
        print(f"{id} already exists", flush=True)
    return True


def get_kegg(target_dir,
             session,
             prefix="D",
             max_idx=12897,
             valid_ids=False):
    """
    Downloads data from the KEGG database and saves it to the specified directory.

    Parameters:
    target_dir (str): The directory to save the downloaded data.
    session (requests.Session): The session object to use for making requests.
    prefix (str, optional): The prefix for the KEGG IDs to download. Defaults to "D".
    max_idx (int, optional): The maximum index to download. Defaults to 12897.
    valid_ids (bool or list, optional): If False, generate IDs based on the prefix. If a list, use the provided valid IDs. Defaults to False.

    Returns:
    None
    """
    bad_file = f"../data/{prefix.replace('_full', '')}_IDs_bad.dat"
    if os.path.exists('data/'):
        bad_file = bad_file[3:]
    bad_file = os.path.abspath(bad_file)
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
        os.makedirs(os.path.dirname(bad_file), exist_ok=True)
        molless_ids = []
        invalid_ids = []

    # Start the loop
    i = 0
    while True:
        # Get the formatted id
        if not valid_ids:
            i += 1
            id = make_custom_id(i, prefix=prefix)
        else:
            id = valid_ids[i]
            i += 1
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
        if (max_idx is not None and i >= max_idx) or (valid_ids and i >= len(valid_ids)):
            print(f"Maximum index reached {max_idx}", flush=True)
            break

    # Write a log file
    with open(bad_file, "w") as f:
        f.write("id,reason\n")
        for id in molless_ids:
            f.write(f"{id},molless\n")
        for id in invalid_ids:
            f.write(f"{id},invalid\n")

    return None


def get_kegg_all(target_dir="kegg_data",
                 target="C",
                 skip=None):
    """
    Retrieves all data from the KEGG database for a specified target and saves it to the specified directory.

    Parameters:
    target_dir (str, optional): The directory to save the downloaded data. Defaults to "kegg_data".
    target (str, optional): The target type to download from the KEGG database. Can be "D", "D_full", "C", "C_full", or "R". Defaults to "C".
    skip (list, optional): A list of IDs to force skip during the download process. Defaults to None.

    Returns:
    None

    Raises:
    ValueError: If the target is unknown.
    """
    # make the session
    session = prepare_session()
    # Check if the target is a valid target and get the maximum index
    if target == "D" or target == "D_full":
        _, max_idx, valid_ids = get_total_n(session, database="drug")
    elif target == "C" or target == "C_full":
        _, max_idx, valid_ids = get_total_n(session, database="compound")
    elif target == "R":
        _, max_idx, valid_ids = get_total_n(session, database="reaction")
    else:
        raise ValueError(f"Unknown target {target}")

    print(f"Total number of current entries {len(valid_ids)}", flush=True)

    if skip is not None:
        valid_ids = sorted([id for id in valid_ids if id not in skip])
        print(f"Checking data for {len(valid_ids)} entries", flush=True)

    # Get the data
    get_kegg(os.path.abspath(target_dir + f"_{target}"),
             session,
             prefix=target,
             max_idx=max_idx,
             valid_ids=valid_ids)


def download_data(target="R",
                  target_dir=r"../../data/kegg_data",
                  skip=None):
    """
    Downloads and cleans KEGG data for a specified target.

    Parameters:
    target (str, optional): The target type to download from the KEGG database. Defaults to "R".
    target_dir (str, optional): The directory to save the downloaded data. Defaults to "../../data/kegg_data".
    skip (list, optional): A list of IDs to force skip during the download process. Defaults to None.

    Returns:
    pd.DataFrame: A DataFrame containing the entries found in the downloaded data.
    """
    target_dir = os.path.abspath(target_dir)
    if not os.path.exists(f"{target_dir}_{target}"):
        os.makedirs(f"{target_dir}_{target}")
    # Clean the data
    clean_empty_folders(f"{target_dir}_{target}")
    # Get the data
    get_kegg_all(target_dir=target_dir, target=target, skip=skip)
    # Clean the data
    clean_empty_folders(f"{target_dir}_{target}")
    # Note attributes found
    entries = [f for f in file_list_all(f"{target_dir}_{target}") if not os.path.basename(f).startswith('.')]
    if target in ['R', 'C_full']:
        entries = log_attributes_found(f"{target_dir}_{target}")
    elif target == 'C':
        entries = pd.DataFrame({f.split('/')[-2]: {'mol_file': f} for f in entries if f.endswith('.mol')}).T
    else:
        entries = pd.DataFrame({f.split('/')[-2]: {'file': f} for f in entries}).T
    return entries


def log_attributes_found(target_dir=r'../../data/kegg_data_C_full'):
    """
    Checks .data files for available attributes. 
    For compounds, ATOM and BOND fields ID molless files before calling API. BRACKET helps ID n/m/x/y/z/w coefficients.

    Parameters:
    target_dir (str): The directory containing the .data files.

    Returns:
    pd.DataFrame: a boolean DataFrame keyed on database entry IDs.
    """
    today = time.strftime('%Y%m%d')
    out_file = os.path.abspath(f'{target_dir}_attributes_{today}.csv')
    paths = {f.split('/')[-2]: f for f in file_list_all(target_dir) if f.endswith('.data')}

    compound_fields = ['ATOM', 'BOND', 'BRACKET', 'BRITE', 'COMMENT', 'DBLINKS', 'ENTRY',
                       'ENZYME', 'EXACT_MASS', 'FORMULA', 'GENE', 'MODULE', 'MOL_WEIGHT',
                       'NAME', 'NETWORK', 'ORGANISM', 'ORIGINAL', 'PATHWAY', 'REACTION',
                       'REMARK', 'REPEAT', 'SEQUENCE', '  TYPE']
    reaction_fields = ['DEFINITION', 'ENTRY', 'EQUATION', 'ENZYME', 'BRITE', 'RCLASS', 'NAME',
                       'PATHWAY', 'ORTHOLOGY', 'DBLINKS', 'COMMENT', 'MODULE', 'AUTHORS',
                       'REFERENCE', 'JOURNAL', 'TITLE', 'REMARK']

    print('Logging attributes found in data files...')
    has_field = dict()
    for id, path in paths.items():
        txt = open(path, 'r').read()
        has_field[id] = {s.strip(): txt.__contains__(f'\n{s}') for s in compound_fields + reaction_fields}

    df = pd.DataFrame(has_field).T
    df = df[df[df].dropna(axis=1, how='all').count().sort_values(ascending=False).index].sort_index()

    df.astype(int).to_csv(out_file)
    return df.assign(path=paths)
