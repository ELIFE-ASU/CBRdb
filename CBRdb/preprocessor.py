import csv
import os
import re

import pandas as pd
from rdkit import Chem as Chem
from rdkit import RDLogger

from .tools_mols import standardize_mol, get_mol_descriptors, fix_r_group

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

from .tools_files import file_list_all, remove_filepath, delete_files_substring
from .tools_eq import standardise_eq


def load_csv_to_dict(file_path):
    """
    Loads a CSV file and converts it to a dictionary.

    Parameters:
    file_path (str): The path to the CSV file.

    Returns:
    dict: A dictionary where the keys are the first column values and the values are the second column values.
    """
    result_dict = {}
    with open(file_path, mode='r', encoding='utf-8') as file:
        reader = csv.reader(file)
        for row in reader:
            if len(row) >= 2:  # Ensure there are at least two columns
                key, value = row[0].strip(), row[1].strip()
                result_dict[key] = value
    return result_dict


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


def convert_mol_to_smiles(target_dir, man_dict, outfile="kegg_data_C.csv.zip", calc_info=True, n_print=100):
    # Get a list of all files in the directory
    files = file_list_all(target_dir)
    # Clean up the files
    delete_files_substring(target_dir, "_r")
    delete_files_substring(target_dir, "_p")
    # Filter the files to only include .mol files
    files = [f for f in files if "_r" not in f]
    files = [f for f in files if "_p" not in f]

    # X group files
    ids_x = []
    ids_similes_fail = []

    # Get the number of files
    n = len(files)
    # Create lists to store the outputs
    arr_smiles = []
    arr_cid = []
    arr_formula = []
    arr_mw = []
    arr_n_heavy = []
    arr_nc = []
    # Loop over the files
    for i, file in enumerate(files):
        # Get the CID
        cid = os.path.basename(file).split(".")[0]
        if i % n_print == 0:
            print(f"Processing file {i}/{n}: {cid}", flush=True)
        # Init flags
        f_load_r = None
        f_load_p = None
        if check_for_x_group(file):
            print(f"Skipping {cid} due to X group", flush=True)
            ids_x.append(cid)
            continue

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
        mol = Chem.MolFromMolFile(file, sanitize=False, removeHs=False)
        # Remove the temporary files
        if flag_r:
            remove_filepath(f_load_r)
        if flag_p:
            remove_filepath(f_load_p)

        try:
            # Standardize and embed the molecule
            mol = standardize_mol(mol)
            # Fix the fix_r_group
            mol = fix_r_group(mol)
            # Convert the molecule to smiles
            smi = Chem.MolToSmiles(mol)
            # Add the ID
            arr_cid.append(cid)
            # Add the smiles to the array
            arr_smiles.append(smi)
            if calc_info:
                formula, mw, n_heavy, nc = get_mol_descriptors(mol)
                # Get the formula
                arr_formula.append(formula)
                # Get the molecular weight
                arr_mw.append(mw)
                # Get the number of heavy atoms
                arr_n_heavy.append(n_heavy)
                # Get the chirality
                arr_nc.append(nc)
        except:
            print(f"Error in {cid}, could not pass to SIMILES", flush=True)
            ids_similes_fail.append(cid)

    # Loop over the manual fixes and add them to the list
    for cid, smiles in man_dict.items():
        # Check if the CID is not already in the list
        if cid not in arr_cid:
            print(f"Adding manual fix for {cid}", flush=True)
            # Add the ID
            arr_cid.append(cid)
            mol = Chem.MolFromSmiles(smiles)
            # Standardize and embed the molecule
            mol = standardize_mol(mol)
            # Add the smiles to the array
            arr_smiles.append(Chem.MolToSmiles(mol))
            if calc_info:
                formula, mw, n_heavy, nc = get_mol_descriptors(mol)
                # Get the formula
                arr_formula.append(formula)
                # Get the molecular weight
                arr_mw.append(mw)
                # Get the number of heavy atoms
                arr_n_heavy.append(n_heavy)
                # Get the chirality
                arr_nc.append(nc)

    # Create a dataframe
    df = pd.DataFrame(data={
        "compound_id": arr_cid,
        "smiles": arr_smiles,
        "formula": arr_formula,
        "molecular_weight": arr_mw,
        "n_heavy_atoms": arr_n_heavy,
        "n_chiral_centers": arr_nc})
    # Sort the dataframe by the compound ID
    df = df.sort_values(by="compound_id")
    # Save the dataframe
    df.to_csv(outfile, compression='zip', encoding='utf-8', index=False)

    # Save the problematic CIDs
    with open("../data/C_IDs_bad.dat", "a") as f:
        for cid in ids_x:
            f.write(f"{cid}, X group\n")
        for cid in ids_similes_fail:
            f.write(f"{cid}, SIMILES fail\n")


def preprocess_kegg_r(target_dir, outfile, n_print=100):
    # Get a list of all files in the directory
    paths = [m for n in [[f'{i}/{k}' for k in j] for i, _, j in list(os.walk(target_dir))[1:]] for m in n]

    # Import reaction data (takes < 20 seconds)
    df = pd.DataFrame({os.path.basename(path).split(".")[0]:  # for each reaction ID
                           pd.read_fwf(path, colspecs=[(0, 12), (12, -1)], header=None,
                                       names=['id', 'line'])  # read file
                      .ffill(axis=0).set_index('id')  # indented lines relate to last-appearing header
                           ['line'].str.strip().groupby(level=0).apply(' '.join)  # combine all lines for each header
                       for path in paths}).drop('///').T  # indexes are reaction IDs; cols are info types
    df = df.set_axis(df.columns.str.strip().str.lower(), axis=1).drop(  # remove columns not needed currently
        ['reference', 'authors', 'journal', 'title', 'brite', 'definition'], axis=1)

    # Remove reactions with glycan IDs mixed in. "remark" column tags their equivalent reactions.
    df = df.loc[df['equation'].str.count(r"(\bG\d{5}\b)") == 0]

    # Store observed KO definitions in a file; old versions of this are used to annotate JGI (meta)genomes.
    ko_defs = df['orthology'].dropna().drop_duplicates()
    ko_defs = pd.Series(dict(zip(ko_defs.str.findall(r"(\bK\d{5}\b)").explode(),
                                 ko_defs.str.split(r"\bK\d{5}\b").apply(lambda x: x[1:]).explode())))
    ko_defs.to_csv(outfile.replace('.csv.zip', '_kodefs.csv.zip'), compression='zip', encoding='utf-8')
    del ko_defs

    # Extract reaction attributes and linkages
    df['reaction'] = df['equation'].apply(standardise_eq)  # standardize reaction formatting
    df['ec'] = df['enzyme'].fillna(' ').str.split().map(' '.join)  # combine all ECs, including partials

    patterns = {'orthology': r"(\bK\d{5}\b)", 'pathway': r"(\brn\d{5}\b)", 'module': r"(\bM\d{5}\b)",
                'rclass': r"(\bRC\d{5}\b  \bC\d{5}_C\d{5})", 'dblinks': r"( \d{5})", 'entry': 'Overall'}
    [df.update(df[k].str.findall(v).map(' '.join, na_action='ignore')) for k, v in patterns.items()]

    # Rename columns where appropriate
    df.rename(columns={'dblinks': 'rhea', 'entry': 'category'}, inplace=True)
    df['category'] = df['category'].replace('', float('nan'))
    df = df.loc[:, df.count().sort_values(ascending=False).index].drop(columns='enzyme')
    df.reset_index().to_csv(outfile, compression='zip', encoding='utf-8', index=False)
    return None


def preprocess(target="R",
               target_dir=r"../../data/kegg_data",
               out_file=r"../data/kegg_data",
               cid_manual_file=r"../data/C_IDs_manual.dat"):
    """
    Preprocesses KEGG data based on the specified target type.

    Parameters:
    target (str, optional): The target type to preprocess ('C' for compounds, 'R' for reactions). Defaults to 'R'.
    target_dir (str, optional): The directory containing the KEGG data. Defaults to "../../data/kegg_data".
    out_file (str, optional): The output file path for the preprocessed data. Defaults to "../data/kegg_data".
    cid_manual_file (str, optional): The file path for manual CID fixes. Defaults to "../data/C_IDs_manual.dat".

    Returns:
    None
    """
    # Set absolute paths
    target_dir = os.path.abspath(target_dir + f"_{target}")
    out_file = os.path.abspath(f"{out_file}_{target}.csv.zip")
    cid_manual_file = os.path.abspath(cid_manual_file)

    if target == "C":
        # Defines a dictionary of manual fixes
        man_dict = load_csv_to_dict(cid_manual_file)
        # Defines a list of bad CIDs to skip
        convert_mol_to_smiles(target_dir, man_dict, outfile=out_file)
        print("C preprocessing done", flush=True)
    elif target == "R":
        preprocess_kegg_r(target_dir, out_file)
        print("R preprocessing done", flush=True)
