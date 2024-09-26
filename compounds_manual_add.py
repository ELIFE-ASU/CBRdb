import os

import pandas as pd


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


def main():
    molless_path = 'Data/C_IDs_molless.dat'
    data = pd.read_csv(molless_path, sep='\t').values.flatten()
    print("Data loaded", flush=True)
    print("Data shape", data.shape, flush=True)
    print("Data head", data[:4], flush=True)

    # Target directory
    target_dir = os.path.abspath('../data/kegg_data_C_full')

    # Remove empty folders
    clean_empty_folders(target_dir)
    # List all files in the target directory
    files = file_list(target_dir)
    print(f"Files in the target directory: {files}", flush=True)

    # Good files to follow-up on
    good_file = "Data/C_IDs_good.dat"
    r_list = []
    good_list = []
    # loop over the data
    for i, id in enumerate(data):
        # print(f"ID: {id}", flush=True)
        # load the data
        f_path = os.path.abspath(f'{target_dir}/{id}/{id}.data')
        with open(f_path, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                # Check that the cid has an associated reaction
                if line.startswith('REACTION'):
                    r_list.append(id)

    # Let us refine the list by removing the generic compounds
    for i, id in enumerate(r_list):
        good_flag = True
        # load the data
        f_path = os.path.abspath(f'{target_dir}/{id}/{id}.data')
        with open(f_path, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip().lower()
                # Check if there is a comment line
                if line.startswith('comment'):
                    # if "Generic compound" in line:
                    #     print(f"ID: {id} is a generic compound", flush=True)
                    #     good_flag = False
                    if "protein" in line:
                        print(f"ID: {id} is a protein", flush=True)
                        good_flag = False
                if line.startswith('name'):
                    # search if substring is in the line
                    if "glycan" in line:
                        print(f"ID: {id} is a glycan", flush=True)
                        good_flag = False
                    if "protein" in line:
                        print(f"ID: {id} is a protein", flush=True)
                        good_flag = False
                    if "peptide" in line:
                        print(f"ID: {id} is a peptide", flush=True)
                        good_flag = False
                    if "rna" in line:
                        print(f"ID: {id} is a RNA", flush=True)
                        good_flag = False
                    if "dna" in line:
                        print(f"ID: {id} is a DNA", flush=True)
                        good_flag = False
                    if "lase" in line:
                        print(f"ID: {id} is a enzyme", flush=True)
                        good_flag = False
                    if "steroid" in line:
                        print(f"ID: {id} is a steroid", flush=True)
                        good_flag = False
                    if "lipid" in line:
                        print(f"ID: {id} is a lipid", flush=True)
                        good_flag = False
                    if "lignin" in line:
                        print(f"ID: {id} is a lignin", flush=True)
                        good_flag = False

                if line.startswith("sequence"):
                    print(f"ID: {id} is a sequence", flush=True)
                    good_flag = False
            if good_flag:
                # append the good list
                good_list.append(id)

    # Now we have a list of all the reactions
    good_list = list(set(good_list))
    good_list.sort()
    print(f"Good list: {good_list}", flush=True)
    print(f"Good list length: {len(good_list)}", flush=True)
    # Save the good list
    with open(good_file, "w") as f:
        f.write("# Good compound IDs to follow up on\n")
        for id in good_list:
            f.write(f"{id}\n")


if __name__ == "__main__":
    print("Program started", flush=True)
    main()
    print("Program finished", flush=True)
