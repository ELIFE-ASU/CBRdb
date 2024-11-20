import os
import shutil


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


def file_list_all(mypath=None):
    """
    This function generates a list of all files in a specified directory and its subdirectories.
    If no directory is specified, it defaults to the current working directory.

    Parameters:
    mypath (str, optional): The path to the directory. Defaults to None, which means the current working directory.

    Returns:
    list: A list of all files in the specified directory and its subdirectories.
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


def remove_filepath(path):
    """ param <path> could either be relative or absolute. """
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
    This function deletes all files in a directory and its subdirectories that contain a specified substring.

    Parameters:
    target_dir (str): The path to the directory to be searched.
    substring (str): The substring to be searched for in the file names.

    Returns:
    None
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
