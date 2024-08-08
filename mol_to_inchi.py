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
            remove(file)
            count += 1
    print(f"Deleted {count} files", flush=True)
    return count


def convert_mol_to_inchi(target_dir, bad_list, man_dict, outfile="kegg_inchi_C.csv"):
    # Get a list of all files in the directory
    files = file_list_all(target_dir)
    # Clean up the files
    delete_files_substring(target_dir, "_r")
    delete_files_substring(target_dir, "_p")
    # Filter the files to only include .mol files
    files = [f for f in files if "_r" not in f]
    files = [f for f in files if "_p" not in f]
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
            mol = Chem.MolFromSmiles(man_dict[cid])
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
        try:
            mol = standardize_mol(mol)
        except:
            print(f"Failed to standardize {cid}", flush=True)
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
    bad_list = ["C00210", "C02202","C13681","C13932","C18368","C19040", "C21014", "C22680"]
    man_dict = {
        "C02202": r"[C@H](Cc1ccccc1)(C(=O)N[C@H](Cc1ccccc1)C[N+]#[NH-])NCc1ccccc1",
        "C00210": r"[Co+](N1C2[C@@]3(N=C([C@H]([C@@]3(CC(=O)N)C)CCC(=O)N)C(=C3N=C([C@H]([C@@]3(CC(=O)N)C)CCC(=O)N)C=C3N=C(C(=C1[C@@]([C@H]2CC(=O)N)(CCC(=O)NC[C@H](OP(=O)(O[C@H]1[C@H]([C@H](O[C@@H]1CO)[*])O)[O-])C)C)C)[C@H](C3(C)C)CCC(=O)N)C)C)[*]",
        "C00462": r"[*H]",
        "C19040": r"[Si-2](F)(F)(F)(F)(F)F.[Na+].[Na+]",
        "C01365": r"c12[nH]cc(c1cccc2)C[C@@H](C(=O)O)N[*]",
        "C01706": r"C(=C\[*])/[*]",
        "C01812": r"[*]CC(=O)O",
        "C01813": r"C(O)([*])[*]",
        "C13681": r"c1c(ccc(c1)N(C)C)[N+]#N.[B+3]([F-])([F-])([F-])[F-]",
        "C13932": r"N.N.N.N.[Ru+10]([O-2][Ru])[O-2][Ru].N.N.N.N.N.N.N.N.N.N.[Cl-].[Cl-].[Cl-].[Cl-].[Cl-].[Cl-]",
        "C18368": r"O=[Cl]=O",
        "C19040": r"[Si-2](F)(F)(F)(F)(F)F.[Na+].[Na+]"}
    convert_mol_to_inchi(target_dir, bad_list, man_dict)
    print("Program finished", flush=True)
