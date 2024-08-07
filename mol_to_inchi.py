import os
import shutil

import numpy as np
import pandas as pd
from rdkit.Chem import AllChem as Achem
from rdkit import Chem as Chem

from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Draw
from rdkit import RDLogger

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)


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


def check_for_R_group(target_file, re_target=["R# ", "R ", "* "]):
    with open(target_file, "r") as f:
        lines = [line for line in f if "M  " not in line]
        return any(target in line for line in lines for target in re_target)


def replace_r_group(target_file, new_file, re_atom="C", re_target=["R# ", "R ", "* "]):
    with open(target_file, "r") as f, open(new_file, "w") as nf:
        lines = [line for line in f if "M  " not in line]
        for target in re_target:
            lines = [line.replace(target, re_atom + "  " if len(target) == 2 else re_atom + " ") for line in lines]
        for line in lines[:-1]:
            nf.write(line)
        nf.write("M  END\n\n")


def check_for_problem_group(target_file, re_target=["OH"]):
    with open(target_file, "r") as f:
        lines = [line for line in f if "M  " not in line]
        return any(target in line for line in lines for target in re_target)


def replace_problem_group(target_file, new_file, re_target=["OH"]):
    with open(target_file, "r") as f, open(new_file, "w") as nf:
        lines = [line for line in f if "M  " not in line]
        for target in re_target:
            if target in "OH":
                lines = [line.replace(target, "O  ") for line in lines]
        for line in lines[:-1]:
            nf.write(line)
        nf.write("M  END\n\n")


def check_for_X_group(target_file):
    with open(target_file, "r") as f:
        for line in f:
            if "X" in line:
                return True


def embed_mol(mol):
    # Sanitise the molecule
    Achem.SanitizeMol(mol, catchErrors=True)
    # Update the properties
    mol.UpdatePropertyCache()
    # kekulize the molecule
    Achem.Kekulize(mol)
    # Add hydrogens
    mol = Chem.AddHs(mol)
    return mol


def standardize_mol(mol):
    # Standardize the molecule
    nrm = rdMolStandardize.Normalizer()
    nrm.normalizeInPlace(mol)
    # Embed the molecule
    # mol = embed_mol(mol)
    return mol


def remove(path):
    """ param <path> could either be relative or absolute. """
    if os.path.isfile(path) or os.path.islink(path):
        os.remove(path)  # remove the file
    elif os.path.isdir(path):
        shutil.rmtree(path)  # remove dir and all contains
    else:
        raise ValueError("file {} is not a file or dir.".format(path))


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
            remove(file)
            count += 1
    print(f"Deleted {count} files", flush=True)
    return count


def convert_mol_to_inchi(target_dir, bad_list, man_dict, outfile="kegg_inchi_C.csv"):
    # Get a list of all files in the directory
    files = file_list_all(target_dir)
    # Clean up the files
    delete_files_substring(target_dir, "mod")
    # Filter the files to only include .mol files
    files = [f for f in files if "mod" not in f]
    # Create arrays to store the InChI and CID
    n_data = len(files)
    arr_inchi = []
    arr_cid = []
    # Loop over the files
    for i, file in enumerate(files):
        # Get the CID
        cid = os.path.basename(file).split(".")[0]
        print(f"Processing file {i + 1}/{len(files)}: {cid}", flush=True)
        # Check if the CID is in the bad list
        if cid in bad_list:
            continue
        # Check if the CID is in the manual fix dictionary
        if cid in man_dict:
            mol = Chem.MolFromInchi(man_dict[cid])
        else:
            # check for R groups
            flag_r = check_for_R_group(file)
            if flag_r:
                f_load_r = file.split(".")[0] + "_r.mol"
                replace_r_group(file, f_load_r)
                file = f_load_r
            # check for problem groups
            flag_p = check_for_problem_group(file)
            if flag_p:
                f_load_p = file.split(".")[0] + "_p.mol"
                replace_problem_group(file, f_load_p)
                file = f_load_p
            # get the molecule
            mol = Chem.MolFromMolFile(file, removeHs=True, sanitize=True)
        # remove the temporary files
        if flag_r:
            remove(f_load_r)
        if flag_p:
            remove(f_load_p)
        # Standardize and embed the molecule
        mol = standardize_mol(mol)
        # Get the InChI
        arr_inchi.append(Chem.MolToInchi(mol))
        arr_cid.append(cid)

    # Create a dataframe
    df = pd.DataFrame(data={"CID": arr_cid, "InChI": arr_inchi})
    # Remove any entries with zero
    df = df[df["InChI"] != "0"]
    # Save the dataframe
    df.to_csv(outfile, index=False)


if __name__ == "__main__":
    print("Program started", flush=True)
    target_dir = os.path.abspath(r"C:/Users/louie/skunkworks/data/kegg_data_C")
    # mostly valence errors
    bad_list = ["C01322", "C02202", "C00210", "C00462", "C19040", "C21014", "C22680", "C01365", "C01706", "C01812",
                "C01813"]
    man_dict = {
        "C02202": "InChI=1S/C27H26N4O4/c28-29-18-25(32)23(16-20-10-4-1-5-11-20)30-26(33)24(17-21-12-6-2-7-13-21)31-27(34)35-19-22-14-8-3-9-15-22/h1-15,18,23-24H,16-17,19H2,(H,30,33)(H,31,34)/t23-,24-/m0/s1",
        "C13681": "InChI=1S/C8H10N3.BF4/c1-11(2)8-5-3-7(10-9)4-6-8;2-1(3,4)5/h3-6H,1-2H3;/q+1;-1",
        "C13932": "InChI=1S/6ClH.14H3N.2O.3Ru/h6*1H;14*1H3;;;;;/q;;;;;;;;;;;;;;;;;;;;;;3*+2/p-6",
        "C18368": "InChI=1S/ClO2/c2-1-3",
        "C19040": "InChI=1S/F6Si.2Na/c1-7(2,3,4,5)6;;/q-2;2*+1"}
    convert_mol_to_inchi(target_dir, bad_list, man_dict)
    print("Program finished", flush=True)
