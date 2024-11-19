import os

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
