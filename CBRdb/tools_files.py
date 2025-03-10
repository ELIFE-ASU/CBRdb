import os
import shutil

import pandas as pd


def file_list(mypath=None):
    """
    Returns a list of all files in the specified directory.

    If no directory path is provided, the current working directory is used.

    Parameters:
    mypath (str): The path to the directory. Defaults to None.

    Returns:
    list: A list of file names in the specified directory.
    """
    mypath = mypath or os.getcwd()  # If no path is provided, use the current working directory
    return [f for f in os.listdir(mypath) if
            os.path.isfile(os.path.join(mypath, f))]  # Return a list of all files in the directory


def file_list_all(mypath=None):
    """
    Returns a list of all files in the specified directory and its subdirectories.

    If no directory path is provided, the current working directory is used.

    Parameters:
    mypath (str): The path to the directory. Defaults to None.

    Returns:
    list: A list of file paths in the specified directory and its subdirectories.
    """
    mypath = mypath or os.getcwd()  # If no path is provided, use the current working directory
    files = []
    # os.walk generates the file names in a directory tree by walking the tree either top-down or bottom-up
    for dirpath, dirnames, filenames in os.walk(mypath):
        for filename in filenames:
            # os.path.join joins one or more path components intelligently
            files.append(os.path.join(dirpath, filename))
    return files


def list_empty_folders(target_dir):
    """
    Returns a list of empty folders in the specified directory.

    Parameters:
    target_dir (str): The path to the target directory.

    Returns:
    list: A list of paths to empty folders in the specified directory.
    """
    empty_folders = []
    for entry in os.listdir(target_dir):
        entry_path = os.path.join(target_dir, entry)
        if os.path.isdir(entry_path) and not os.listdir(entry_path):
            empty_folders.append(entry_path)
    return empty_folders


def clean_empty_folders(target_dir, size=False):
    """
    Removes empty folders in the specified directory.

    Parameters:
    target_dir (str): The path to the target directory.
    size (bool): If True, also considers folders with size 0 as empty. Defaults to False.

    Returns:
    int: The number of empty folders removed.
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


def remove_filepath(path):
    """
    Removes a file or directory at the specified path.

    Parameters:
    path (str): The path to the file or directory to be removed.

    Raises:
    ValueError: If the specified path is neither a file nor a directory.
    Exception: If there is an error during the removal process.
    """
    try:
        if os.path.isfile(path) or os.path.islink(path):
            os.remove(path)  # remove the file
        elif os.path.isdir(path):
            shutil.rmtree(path)  # remove dir and all contains
        else:
            raise ValueError("file {} is not a file or dir.".format(path))
    except Exception as e:
        print("Failed to delete", path, e, flush=True)


def delete_files_substring(target_dir, substring):
    """
    Deletes files in the specified directory and its subdirectories that contain the given substring in their names.

    Parameters:
    target_dir (str): The path to the target directory.
    substring (str): The substring to search for in file names.

    Returns:
    int: The number of files deleted.
    """
    # Get a list of all files in the directory and its subdirectories
    files = file_list_all(target_dir)
    # Iterate through the files
    count = 0
    for file in files:
        # Check if the substring is in the file name
        if substring in file:
            # Delete the file
            remove_filepath(file)
            count += 1
    print(f"Deleted {count} files", flush=True)
    return count


def make_custom_id(idx, prefix="C", digits=5):
    """
    Generates a custom ID with a specified prefix and number of digits.

    Parameters:
    idx (int): The index to be included in the ID.
    prefix (str): The prefix for the ID. Default is "C".
    digits (int): The number of digits for the ID. Default is 5.

    Returns:
    str: The generated custom ID.
    """
    return f"{prefix}{int(idx):0{digits}d}"


def add_suffix_to_file(fname, suffix):
    """
    Adds a suffix to the base name of a file before the file extension.

    Parameters:
    fname (str): The full path to the file.
    suffix (str): The suffix to add to the file name.

    Returns:
    str: The new file name with the suffix added.
    """
    dn = os.path.dirname(fname)  # Get the directory name from the file path
    bn = os.path.basename(fname)  # Get the base name of the file
    fn, ext = bn[:bn.find(".")], bn[bn.find("."):]  # Split the base name into the file name and extension
    return f"{dn}/{fn}_{suffix}{ext}"  # Return the new file name with the suffix added


def reaction_csv(df_R: pd.DataFrame, file_address: str):
    """
    Prints a reaction dataframe to CSV file, in standardized format.

    Parameters:
    df_R (pd.DataFrame): The reaction DataFrame to be printed as a CSV.
    file_address (str): The file address i.e. where to save the CSV.

    Returns:
    str : The file address where the CSV was saved.
    """
    df = df_R.copy(deep=True)
    file_address = os.path.abspath(file_address)
    params = {'encoding': 'utf-8', 'index': False, 'float_format': '%.3f'}
    listlike_cols = df.columns.intersection(['ec', 'rclass', 'pathway', 'orthology', 'rhea', 'module'])
    col_order = ['id', 'reaction'] + sorted(list(listlike_cols)) + sorted(
        df.columns.difference(listlike_cols.union(['id', 'reaction'])))
    for col in listlike_cols:
        first_entry = df[col].dropna().iloc[0]
        if type(first_entry) is str:
            df.loc[:, col] = df.loc[:, col].str.replace('  ', '__').str.split(' ')
        if hasattr(first_entry, '__iter__'):
            df.loc[:, col] = df.loc[:, col].map(lambda x: ' '.join(sorted(list(x))), na_action='ignore')
    df = df.sort_values(by='id').reset_index(drop=True).loc[:, col_order]
    df.to_csv(file_address, **params)
    return file_address


def compound_csv(df_C: pd.DataFrame, file_address: str):
    """
    Prints a reaction dataframe to CSV file, in standardized format.

    Parameters:
    df_C (pd.DataFrame): The compound DataFrame to be printed as a CSV.
    file_address (str): The file address i.e. where to save the CSV.

    Returns:
    str : The file address where the CSV was saved.
    """
    df = df_C.copy(deep=True)
    file_address = os.path.abspath(file_address)
    params = {'encoding': 'utf-8', 'index': False, 'float_format': '%.3f'}
    first_cols = ['compound_id', 'smiles', 'formula', 'molecular_weight', 'n_heavy_atoms', 'n_chiral_centers',
                  'smiles_capped', 'inchi_capped', 'name', 'comment', 'glycan_ids', 'drug_ids']
    col_order = [i for i in first_cols if i in df.columns] + sorted(list(df.columns.difference(first_cols)),
                                                                    reverse=True)
    df = df.sort_values(by='compound_id').reset_index(drop=True).loc[:, col_order]
    df.to_csv(file_address, **params)
    return file_address


def count_df_entries(dict_of_dbs, pipeline_stage):
    """
    Counts the number of unique IDs in a dictionary containing reaction and/or compound DataFrames.

    Parameters:
    dict_of_dbs (pd.DataFrame): dictionary containing DataFrames of reactions and/or compounds.
    pipeline_stage (str): stage of the pipeline where the count is being done. 

    Returns:
    pd.DataFrame: DataFrame containing the number of unique IDs in each DataFrame within the dict.
    """
    id_cols = ['id', 'compound_id']
    dbs = {i: j for i, j in dict_of_dbs.items() if type(j) == pd.DataFrame and len(j.columns.intersection(id_cols)) > 0}
    dbs = pd.Series(
        {i: len(j[j.columns.intersection(id_cols)].drop_duplicates().index) for i, j in dbs.items()}).to_frame(
        name=pipeline_stage)
    return dbs
