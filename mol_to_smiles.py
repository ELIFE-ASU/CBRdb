import os
import shutil

import pandas as pd
from rdkit import Chem as Chem
from rdkit import RDLogger
from rdkit.Chem import AllChem as Achem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

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


def check_for_r_group(target_file, re_target=["R# ", "R ", "* "]):
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


def check_for_x_group(target_file):
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
    mol = embed_mol(mol)
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


def convert_mol_to_smiles(target_dir, bad_list, man_dict, outfile="kegg_data_C.csv.zip"):
    # Get a list of all files in the directory
    files = file_list_all(target_dir)
    # Clean up the files
    delete_files_substring(target_dir, "_r")
    delete_files_substring(target_dir, "_p")
    # Filter the files to only include .mol files
    files = [f for f in files if "_r" not in f]
    files = [f for f in files if "_p" not in f]
    # Create lists to store the outputs
    arr_smiles = []
    arr_cid = []
    arr_formula = []
    arr_mw = []
    arr_n_heavy = []
    # print(bad_list, flush=True)
    # Loop over the files
    for i, file in enumerate(files):
        # Get the CID
        cid = os.path.basename(file).split(".")[0]
        print(f"Processing file {i + 1}/{len(files)}: {cid}", flush=True)
        # Init flags
        flag_r = False
        flag_p = False
        f_load_r = None
        f_load_p = None
        # Check if the CID is in the bad list
        if cid in bad_list:
            print(f"Skipping {cid} due to bad list", flush=True)
            continue
        if check_for_x_group(file):
            print(f"Skipping {cid} due to X group", flush=True)
            continue
        # Check if the CID is in the manual fix dictionary
        if cid in man_dict:
            mol = Chem.MolFromSmiles(man_dict[cid])
        else:
            # Check for R groups
            flag_r = check_for_r_group(file)
            if flag_r:
                f_load_r = file.split(".")[0] + "_r.mol"
                replace_r_group(file, f_load_r)
                file = f_load_r
            # Check for problem groups
            flag_p = check_for_problem_group(file)
            if flag_p:
                f_load_p = file.split(".")[0] + "_p.mol"
                replace_problem_group(file, f_load_p)
                file = f_load_p
            # Get the molecule
            mol = Chem.MolFromMolFile(file, removeHs=True, sanitize=True)
        # Remove the temporary files
        if flag_r:
            remove(f_load_r)
        if flag_p:
            remove(f_load_p)
        # Standardize and embed the molecule
        mol = standardize_mol(mol)
        # Get the smiles
        smi = Chem.MolToSmiles(mol)
        try:
            # Ensure the mol can be converted back to mol
            standardize_mol(Chem.MolFromSmiles(smi))
            # Add the smiles to the array
            arr_smiles.append(smi)
            # Get the formula
            arr_formula.append(rdMolDescriptors.CalcMolFormula(mol))
            # Get the molecular weight
            arr_mw.append(rdMolDescriptors.CalcExactMolWt(mol))
            # Get the number of heavy atoms
            arr_n_heavy.append(mol.GetNumHeavyAtoms())
            # Add the ID
            arr_cid.append(cid)
        except:
            print(f"Error in {cid}, could not pass to SIMLES", flush=True)

    # Create a dataframe
    df = pd.DataFrame(data={
        "CID": arr_cid,
        "smiles": arr_smiles,
        "Formula": arr_formula,
        "Molecular Weight": arr_mw,
        "N Heavy": arr_n_heavy})
    # Save the dataframe
    df.to_csv(outfile, compression='zip', encoding='utf-8')


if __name__ == "__main__":
    print("Program started", flush=True)
    target_dir = os.path.abspath(r"C:/Users/louie/skunkworks/data/kegg_data_C")
    out_file = os.path.abspath(r"Data/kegg_data_C.csv.zip")

    list_x = ["C00462",
              "C01365",
              "C01706",
              "C01812",
              "C01813",
              "C01322",
              "C01872",
              "C02103",
              "C03122",
              "C13373",
              "C15564"]
    # Mostly valence errors
    bad_list = ["C02202",
                "C13681",
                "C13932",
                "C18368",
                "C19040",
                "C21014",
                "C22680"] + list_x
    man_dict = {
        "C02202": r"[C@H](Cc1ccccc1)(C(=O)N[C@H](Cc1ccccc1)C[N+]#[NH-])NCc1ccccc1",
        "C00210": r"[Co+](N1C2[C@@]3(N=C([C@H]([C@@]3(CC(=O)N)C)CCC(=O)N)C(=C3N=C([C@H]([C@@]3(CC(=O)N)C)CCC(=O)N)C=C3N=C(C(=C1[C@@]([C@H]2CC(=O)N)(CCC(=O)NC[C@H](OP(=O)(O[C@H]1[C@H]([C@H](O[C@@H]1CO)[*])O)[O-])C)C)C)[C@H](C3(C)C)CCC(=O)N)C)C)[*]",
        "C19040": r"[Si-2](F)(F)(F)(F)(F)F.[Na+].[Na+]",
        "C13681": r"c1c(ccc(c1)N(C)C)[N+]#N.[B+3]([F-])([F-])([F-])[F-]",
        "C13932": r"N.N.N.N.[Ru+10]([O-2][Ru])[O-2][Ru].N.N.N.N.N.N.N.N.N.N.[Cl-].[Cl-].[Cl-].[Cl-].[Cl-].[Cl-]",
        "C18368": r"O=[Cl]=O",
        "C19040": r"[Si-2](F)(F)(F)(F)(F)F.[Na+].[Na+]"}
    convert_mol_to_smiles(target_dir, bad_list, man_dict, outfile=out_file)
    print("Program finished", flush=True)
