import os
import re

import chemparse
import pandas as pd
from rdkit import Chem as Chem
from rdkit import RDLogger
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

from .tools_files import remove_filepath


def sanitize_mol(mol):
    """
    Standardizes and normalizes the given molecule.

    Parameters:
    mol (rdkit.Chem.Mol): The molecule to be sanitized.

    Returns:
    rdkit.Chem.Mol: The sanitized molecule.
    """
    # Standardize the molecule
    mol.UpdatePropertyCache(strict=False)
    Chem.SetConjugation(mol)
    Chem.SetHybridization(mol)
    # Normalize the molecule
    Chem.SanitizeMol(mol, sanitizeOps=(Chem.SANITIZE_ALL ^ Chem.SANITIZE_CLEANUP ^ Chem.SANITIZE_PROPERTIES))
    # Update the properties
    mol.UpdatePropertyCache(strict=False)
    return mol


def standardize_mol(mol):
    """
    Standardizes the given molecule by sanitizing, normalizing, and kekulizing it.

    Parameters:
    mol (rdkit.Chem.Mol): The molecule to be standardized.

    Returns:
    rdkit.Chem.Mol: The standardized molecule.
    """
    # Sanitize the molecule
    mol = sanitize_mol(mol)
    # Normalize the molecule
    rdMolStandardize.NormalizeInPlace(mol)
    # Kekulize the molecule
    Chem.Kekulize(mol)
    # Update the properties
    mol.UpdatePropertyCache(strict=False)
    return mol


def fix_r_group(mol, target="[H:0]", re="*"):
    """
    Replaces all occurrences of a specified target atom or group in the molecule with another specified atom or group.

    Parameters:
    mol (rdkit.Chem.Mol): The molecule to be modified.
    target (str): The atom or group to be replaced (default is "[H:0]").
    re (str): The replacement atom or group (default is "*").

    Returns:
    rdkit.Chem.Mol: The standardized molecule after replacement.
    """
    # Convert the molecule to a SMILES string
    smi = Chem.MolToSmiles(mol)
    # Replace all occurrences of the target with the replacement
    smi = smi.replace(target, re)
    # Convert the SMILES string back to a molecule without sanitizing
    mol = Chem.MolFromSmiles(smi, sanitize=False)
    # Standardize the molecule
    return standardize_mol(mol)


def get_chirality(mol):
    """
    Determines the number of chiral centers in the given molecule.

    Parameters:
    mol (rdkit.Chem.Mol): The molecule to be analyzed.

    Returns:
    int: The number of chiral centers in the molecule.
    float: NaN if an error occurs during the analysis.
    """
    try:
        return len(
            Chem.FindMolChiralCenters(mol,
                                      useLegacyImplementation=False,
                                      includeUnassigned=True,
                                      includeCIP=False))
    except:
        return float('NaN')


def mol_replacer(smi, target="[H]"):
    """
    Replaces all occurrences of the '*' atom in the given SMILES string with the specified target atom or group.

    Parameters:
    smi (str): The SMILES string representing the molecule.
    target (str): The atom or group to replace '*' with (default is "[H]").

    Returns:
    str: The modified SMILES string with all '*' atoms replaced by the target.
    """
    # Convert the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smi)
    # Create a molecule for the '*' atom
    star = Chem.MolFromSmiles('*')
    # Create a molecule for the target atom or group
    tar = Chem.MolFromSmiles(target)
    # Replace all occurrences of '*' with the target in the molecule
    mol_re = Chem.ReplaceSubstructs(mol, star, tar, replaceAll=True)
    # Convert the modified molecule back to a SMILES string and return it
    return Chem.MolToSmiles(mol_re[0])


def get_mol_descriptors(mol):
    """
    Calculates various molecular descriptors for the given molecule.

    Parameters:
    mol (rdkit.Chem.Mol): The molecule to be analyzed.

    Returns:
    tuple: A tuple containing the molecular formula (str), exact molecular weight (float),
           number of heavy atoms (int), and number of chiral centers (int or float).
    """
    mol = sanitize_mol(mol)
    return (rdMolDescriptors.CalcMolFormula(mol),
            rdMolDescriptors.CalcExactMolWt(mol),
            rdMolDescriptors.CalcNumHeavyAtoms(mol),
            get_chirality(mol))


def check_for_r_group(target_file, re_target=None):
    """
    Checks if a target file contains any R groups.

    Parameters:
    target_file (str): The path to the target file.
    re_target (list, optional): A list of R group patterns to search for. Default is ["R# ", "R ", "* "].

    Returns:
    bool: True if any R group pattern is found in the file, False otherwise.
    """
    if re_target is None:
        re_target = ["R# ", "R ", "* "]
    with open(target_file, "r") as f:
        lines = [line for line in f if "M  " not in line]
        return any(target in line for line in lines for target in re_target)


def replace_r_group(target_file, new_file, re_atom="H", re_target=None):
    """
    Replaces R groups in a target file with a specified replacement atom and writes the result to a new file.

    Parameters:
    target_file (str): The path to the target file.
    new_file (str): The path to the new file where the result will be written.
    re_atom (str): The replacement atom to use for R groups. Default is "H".
    re_target (list, optional): A list of R group patterns to search for. Default is ["R# ", "R ", "* "].

    Returns:
    None
    """
    if re_target is None:
        re_target = ["R# ", "R ", "* "]
    with open(target_file, "r") as f, open(new_file, "w") as nf:
        lines = [line for line in f if "M  " not in line]
        for target in re_target:
            lines = [line.replace(target, re_atom + "  " if len(target) == 2 else re_atom + " ") for line in lines]
        for line in lines[:-1]:
            nf.write(line)
        nf.write("M  END\n\n")


def check_for_problem_group(target_file, re_target=None):
    """
    Checks if a target file contains any problem groups.

    Parameters:
    target_file (str): The path to the target file.
    re_target (list, optional): A list of problem group patterns to search for. Default is ["OH"].

    Returns:
    bool: True if any problem group pattern is found in the file, False otherwise.
    """
    if re_target is None:
        re_target = ["OH"]
    with open(target_file, "r") as f:
        lines = [line for line in f if "M  " not in line]
        return any(target in line for line in lines for target in re_target)


def replace_problem_group(target_file, new_file, re_target=None):
    """
    Replaces problem groups in a target file with a specified replacement atom and writes the result to a new file.

    Parameters:
    target_file (str): The path to the target file.
    new_file (str): The path to the new file where the result will be written.
    re_target (list, optional): A list of problem group patterns to search for. Default is ["OH"].

    Returns:
    None
    """
    if re_target is None:
        re_target = ["OH"]
    with open(target_file, "r") as f, open(new_file, "w") as nf:
        lines = [line for line in f if "M  " not in line]
        for target in re_target:
            if target in "OH":
                lines = [line.replace(target, "O  ") for line in lines]
        for line in lines[:-1]:
            nf.write(line)
        nf.write("M  END\n\n")


def contains_x_not_xe(s):
    """
    Checks if a string contains the character 'X' not followed by 'e'.

    Parameters:
    s (str): The input string to search.

    Returns:
    bool: True if 'X' not followed by 'e' is found, False otherwise.
    """
    # Use a regular expression to find "X" not followed by "e"
    return bool(re.search(r'X(?!e)', s))


def check_for_x_group(target_file):
    """
    Checks if a target file contains any 'X' groups.

    Parameters:
    target_file (str): The path to the target file.

    Returns:
    bool: True if any 'X' group is found in the file, False otherwise.
    """
    with open(target_file, "r") as f:
        for line in f:
            if contains_x_not_xe(line):
                return True
    return False


def compound_super_safe_load(file):
    """
    Safely loads a compound from a file, handling special cases such as R groups and problem groups.

    Parameters:
    file (str): The path to the file containing the compound.

    Returns:
    rdkit.Chem.rdchem.Mol or None: The loaded molecule, or None if the file contains an 'X' group.
    """
    # Init flags
    f_load_r = None
    f_load_p = None
    if check_for_x_group(file):
        print(f"X group found {file}", flush=True)
        return None

    # Check for R groups
    flag_r = check_for_r_group(file)
    if flag_r:
        print(f"R group found {file}", flush=True)
        f_load_r = file.split(".")[0] + "_r.mol"
        replace_r_group(file, f_load_r)
        file = f_load_r

    # Check for problem groups
    flag_p = check_for_problem_group(file)
    if flag_p:
        print(f"Problem group found {file}", flush=True)
        f_load_p = file.split(".")[0] + "_p.mol"
        replace_problem_group(file, f_load_p)
        file = f_load_p

    # Get the molecule
    mol = Chem.MolFromMolFile(file, sanitize=False, removeHs=False)
    # Remove the temporary files
    if flag_r:
        remove_filepath(f_load_r)
    if flag_p:
        remove_filepath(f_load_p)
    return mol


def get_properties(mol):
    """
    Standardizes a molecule, fixes R groups, converts it to SMILES, and calculates molecular descriptors.

    Parameters:
    mol (rdkit.Chem.rdchem.Mol): The molecule to process.

    Returns:
    list: A list containing the SMILES string, molecular formula, molecular weight, number of heavy atoms, and number of chiral centers.
    """
    # Standardize and embed the molecule
    mol = standardize_mol(mol)
    # Fix the fix r group so that it can be converted to smiles
    mol = fix_r_group(mol)
    # Convert the molecule to smiles
    smi = Chem.MolToSmiles(mol)
    # Calculate the molecular descriptors
    formula, mw, n_heavy, nc = get_mol_descriptors(mol)
    return [smi, formula, mw, n_heavy, nc]


def filter_compounds_without_star(dataframe):
    """
    Filters the compounds DataFrame to return only the rows where the 'smiles' column does not contain a '*'.

    Parameters:
    dataframe (pd.DataFrame): A DataFrame containing compound data with a 'smiles' column.

    Returns:
    pd.DataFrame: A filtered DataFrame with rows where the 'smiles' column does not contain a '*'.
    """
    return dataframe[~dataframe['smiles'].str.contains('*', regex=False)]


def get_sorted_compounds(c_path="../data/kegg_data_C.csv.zip", filter_star=True):
    """
    Retrieves and sorts compound data from a CSV file, optionally filtering out rows where the 'smiles' column contains a '*'.

    Parameters:
    c_path (str): The path to the CSV file containing compound data.
    filter_star (bool): Whether to filter out rows where the 'smiles' column contains a '*'. Default is True.

    Returns:
    pd.DataFrame: A DataFrame containing the sorted compound data.
    """
    # Get the data
    data_c = pd.read_csv(os.path.abspath(c_path))

    # Filter out the star compounds or has * in the smiles
    if filter_star:
        data_c = filter_compounds_without_star(data_c)

    # Sort by the smile size
    return data_c.sort_values(by="smiles", key=lambda x: x.str.len())


def get_small_compounds(c_path="../data/kegg_data_C.csv.zip", filter_star=True, n=1):
    """
    Retrieves and filters small compounds from a CSV file.

    Parameters:
    c_path (str): The path to the CSV file containing compound data. Default is "../data/kegg_data_C.csv.zip".
    filter_star (bool): Whether to filter out rows where the 'smiles' column contains a '*'. Default is True.
    n (int): The maximum number of heavy atoms for a compound to be considered small. Default is 1.

    Returns:
    pd.DataFrame: A DataFrame containing the filtered small compounds.
    """
    data_c = get_sorted_compounds(c_path=c_path, filter_star=filter_star)

    # Filter compounds by the number of heavy atoms
    if n == 1:
        return data_c[data_c["n_heavy_atoms"] <= n]
    else:
        return data_c[data_c["n_heavy_atoms"] == n]


def get_small_compounds_all(c_path="../data/kegg_data_C.csv.zip", filter_star=True, n=1):
    """
    Retrieves and filters small compounds from a CSV file.

    Parameters:
    c_path (str): The path to the CSV file containing compound data. Default is "../data/kegg_data_C.csv.zip".
    filter_star (bool): Whether to filter out rows where the 'smiles' column contains a '*'. Default is True.
    n (int): The maximum number of heavy atoms for a compound to be considered small. Default is 1.

    Returns:
    pd.DataFrame: A DataFrame containing the filtered small compounds.
    """
    data_c = get_sorted_compounds(c_path=c_path, filter_star=filter_star)

    # Filter compounds by the number of heavy atoms
    return data_c[data_c["n_heavy_atoms"] <= n]


def get_compounds_with_matching_elements(data_c_1, diff_ele_react, diff_ele_prod):
    """
    Filters compounds that contain all the elements present in the union of two dictionaries' keys.

    Parameters:
    data_c_1 (DataFrame): A pandas DataFrame containing compound data with 'compound_id' and 'formula' columns.
    diff_ele_react (dict): A dictionary representing elements in reactants.
    diff_ele_prod (dict): A dictionary representing elements in products.

    Returns:
    list: A list of compound IDs that contain all the element symbols from the union of the keys in diff_ele_react and diff_ele_prod.
    """
    # Get the set of keys in react_ele and prod_ele
    element_symbols = list(set(diff_ele_react.keys()).union(set(diff_ele_prod.keys())))

    def contains_all_elements(formula, elements):
        """
        Checks if a chemical formula contains all specified elements.

        Parameters:
        formula (str): The chemical formula to check.
        elements (list): A list of element symbols to check for in the formula.

        Returns:
        bool: True if all elements are present in the formula, False otherwise.
        """
        # Convert the formula to a dictionary of elements and their counts
        formula_dict = chemparse.parse_formula(formula)
        # Check if all elements are in the formula
        return all(element in formula_dict for element in elements)

    # Filter the DataFrame to get compounds that contain all the element symbols
    filtered_compounds = data_c_1[data_c_1['formula'].apply(lambda x: contains_all_elements(x, element_symbols))]

    # Sort by the number of heavy atoms
    filtered_compounds = filtered_compounds.sort_values(by='n_heavy_atoms')

    # Return the list of compound IDs
    return filtered_compounds['compound_id'].tolist()
