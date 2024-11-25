import os

import pandas as pd
from rdkit import Chem as Chem
from rdkit import RDLogger

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

from .tools_mols import standardize_mol, get_mol_descriptors
from .tools_files import make_custom_id
from .preprocessor import check_for_x_group
from .tools_mp import tp_calc


def make_id_range(data, hal_exp):
    """
    Creates a range of IDs for combinations of data and halogen expansions.

    Parameters:
    data (list): A list of data entries.
    hal_exp (list): A list of halogen expansions.

    Returns:
    list: A list of ranges, each representing a combination of data and halogen expansions.
    """
    n_data = len(data)
    n_hal = len(hal_exp)
    n_comb = n_data * n_hal
    num_range = range(n_comb)
    # Reshape num_range to a 2D array
    return [num_range[i:i + n_hal] for i in range(0, n_comb, n_hal)]


def load_bad_entries(target_dir_c):
    """
    Loads and filters bad entries from a target directory.

    Parameters:
    target_dir_c (str): The path to the target directory containing the files.

    Returns:
    list: A list of filenames that are flagged as bad entries.
    """
    target_dir_c = os.path.abspath(target_dir_c)
    files = os.listdir(target_dir_c)
    files_full = [os.path.join(target_dir_c, f, f + ".mol") for f in files]
    flags = tp_calc(check_for_x_group, files_full)
    return [elem for elem, flag in zip(files, flags) if flag]


def get_reactions_with_substring(reactions_df, substring):
    """
    Filters reactions in a DataFrame that contain a specific substring.

    Parameters:
    reactions_df (DataFrame): A pandas DataFrame containing reaction data.
    substring (str): The substring to search for in the reaction data.

    Returns:
    DataFrame: A DataFrame containing reactions that match the substring.
    """
    return reactions_df[reactions_df['reaction'].str.contains(substring, case=False, na=False)]


def fix_halogen_compounds(
        target_dir_c=r"../../data/kegg_data_C",
        hal_exp=None,
        f_print=True):
    """
    Fixes halogen compounds by expanding halogen placeholders in molecular files and generating new compound IDs.

    Parameters:
    target_dir_c (str): The path to the target directory containing compound data files. Default is '../../data/kegg_data_C'.
    hal_exp (list): A list of halogens to expand. Default is ['F', 'Cl', 'Br', 'I'].
    f_print (bool): Flag to print debug information. Default is True.

    Returns:
    tuple: A tuple containing two dictionaries:
        - cids_dict (dict): A dictionary mapping original compound IDs to new compound IDs.
        - smis_dict (dict): A dictionary mapping original compound IDs to their SMILES representations.
    """
    # Prepare the full path of the files
    target_dir_c = os.path.abspath(target_dir_c)
    # Prepare the halogen list to expand over
    if hal_exp is None:
        hal_exp = ['F', 'Cl', 'Br', 'I']
    # Load the bad compound IDs
    data_bad_id = load_bad_entries(target_dir_c)
    if f_print:
        print(f"bad files ID with halogens: {data_bad_id}", flush=True)

    # Make the combinations of the data and the halogens
    num_range = make_id_range(data_bad_id, hal_exp)

    smis_dict = {}
    cids_dict = {}
    for i in range(len(data_bad_id)):
        # Get the full path of the file
        tmp_file = os.path.join(target_dir_c, data_bad_id[i], data_bad_id[i] + ".mol")
        with open(tmp_file, 'r') as f:
            file_data = f.read()
        # Initialize the list for the current data[i]
        cids_dict[data_bad_id[i]] = []
        smis_dict[data_bad_id[i]] = []
        # Replace the X with the halogen
        for j, hal in enumerate(hal_exp):
            idx = num_range[i][j]
            # Load the molecule
            mol = Chem.MolFromMolBlock(file_data.replace("X", hal))
            # Standardize the molecule
            mol = standardize_mol(mol)
            smi = Chem.MolToSmiles(mol, allHsExplicit=True)
            if f_print:
                print(f"Compound {data_bad_id[i]} with halogen {hal}, idx {idx} -> {smi}", flush=True)
            # Determine the compound id from the ones already given
            # If not found, generate a new one
            cid = {
                "C00462": {"F": "C16487", "Cl": "C01327", "Br": "C13645", "I": "C05590"},
                "C01322": {},  # full expansion required!
                "C01365": {},  # No reaction data!
                "C01706": {},  # No reaction data!
                "C01812": {"F": "C06108", "Cl": "C06755"},
                "C01813": {},  # No reaction data!
                "C01872": {},  # full expansion required?
                "C02103": {},  # full expansion required?
                "C03122": {},  # No reaction data!
                "C15564": {},  # full expansion required?

            }.get(data_bad_id[i], {}).get(hal)

            # Generate the compound id if not found
            if cid is None:
                cid = make_custom_id(99000 + idx, prefix="C")
            cids_dict[data_bad_id[i]].append(cid)
            smis_dict[data_bad_id[i]].append(smi)
    return cids_dict, smis_dict


def merge_halogen_compounds(cids_dict,
                            smis_dict,
                            c_id_file="../data/kegg_data_C.csv.zip",
                            int_file=None,
                            out_file=None):
    """
    Merges halogen compound data with existing compound data, removing duplicates and saving the result.

    Parameters:
    cids_dict (dict): A dictionary mapping original compound IDs to new compound IDs.
    smis_dict (dict): A dictionary mapping original compound IDs to their SMILES representations.
    c_id_file (str): The path to the file containing existing compound data. Default is '../data/kegg_data_C.csv.zip'.
    int_file (str, optional): The path to the intermediate file where merged data will be saved. Default is None.
    out_file (str, optional): The path to the output file where final merged data will be saved. Default is the same as c_id_file.

    Returns:
    None
    """
    if out_file is None:
        out_file = c_id_file
    # Prepare the full path of the files
    c_id_file = os.path.abspath(c_id_file)
    if int_file is not None:
        int_file = os.path.abspath(int_file)
    out_file = os.path.abspath(out_file)

    # Convert the cids_dict and smis_dict to a list
    cids_list = []
    smis_list = []
    for key in cids_dict.keys():
        cids_list.extend(cids_dict[key])
        smis_list.extend(smis_dict[key])

    # Loop over the cids_list
    arr_formula = []
    arr_mw = []
    arr_n_heavy = []
    arr_nc = []
    for i in range(len(cids_list)):
        # Load the molecule
        mol = Chem.MolFromSmiles(smis_list[i])
        # Standardize the molecule
        mol = standardize_mol(mol)
        # Get the formula, molecular weight, number of heavy atoms, and the number of chiral centers
        formula, mw, n_heavy, nc = get_mol_descriptors(mol)
        # Add the data to the arrays
        arr_formula.append(formula)
        arr_mw.append(mw)
        arr_n_heavy.append(n_heavy)
        arr_nc.append(nc)

    # Create a dataframe
    df = pd.DataFrame(data={
        "compound_id": cids_list,
        "smiles": smis_list,
        "formula": arr_formula,
        "molecular_weight": arr_mw,
        "n_heavy_atoms": arr_n_heavy,
        "n_chiral_centers": arr_nc})
    if int_file is not None:
        df.to_csv(int_file, compression='zip', encoding='utf-8', index=False)

    # Load the compounds data
    df_old = pd.read_csv(c_id_file, compression='zip')

    # Merge the dataframes
    df = pd.concat([df_old, df], ignore_index=True)
    # Drop the duplicates
    df = df.drop_duplicates(subset="compound_id")
    # Sort the dataframe by the compound ID
    df = df.sort_values(by="compound_id")
    # Save the dataframe
    df.to_csv(out_file, compression='zip', encoding='utf-8')
    return None


def fix_halogen_reactions(cids_dict,
                          r_id_file="../data/atlas_data_R.csv.zip",
                          int_file=None,
                          out_file=None,
                          f_print=False):
    """
    Fixes halogen reactions by replacing halogen placeholders in reaction equations and generating new reaction IDs.

    Parameters:
    cids_dict (dict): A dictionary mapping original compound IDs to new compound IDs.
    r_id_file (str): The path to the file containing existing reaction data. Default is '../data/atlas_data_R.csv.zip'.
    int_file (str, optional): The path to the intermediate file where processed reaction data will be saved. Default is None.
    out_file (str, optional): The path to the output file where final processed reaction data will be saved. Default is the same as r_id_file.
    f_print (bool): Flag to print debug information. Default is False.

    Returns:
    None
    """
    if out_file is None:
        out_file = r_id_file
    # Prepare the full path of the files
    r_id_file = os.path.abspath(r_id_file)
    if int_file is not None:
        int_file = os.path.abspath(int_file)
    out_file = os.path.abspath(out_file)
    idx_base = 99000
    prefix = "R"
    if "atlas" in r_id_file and "kegg" not in r_id_file:
        idx_base = 990000
        prefix = "A"

    # Prepare the full path of the files
    r_id_file = os.path.abspath(r_id_file)
    if f_print:
        print(f"cids_dict {cids_dict}", flush=True)

    # Load the bad compound IDs
    data_bad_id = list(cids_dict.keys())
    n_compounds = len(data_bad_id)
    if f_print:
        print(f"{n_compounds} bad compound files ID with halogens: {data_bad_id}", flush=True)

    # Load the reactions data
    reactions = pd.read_csv(r_id_file, compression='zip')
    if f_print:
        print(reactions.columns, flush=True)
        print(reactions.head(), flush=True)

    # Find the reactions with the halogens
    reactions_set = set()
    for i in range(len(data_bad_id)):
        # Get the equations
        equations = get_reactions_with_substring(reactions, data_bad_id[i])
        reactions_set.update(equations['id'].values)
    reactions_set = list(reactions_set)
    n_reactions = len(reactions_set)
    if f_print:
        print(f"{n_reactions} Reactions with the halogens {reactions_set}", flush=True)

    # Get the EC numbers of the reactions_set
    ec_numbers = reactions[reactions['id'].isin(reactions_set)]['ec'].values
    if f_print:
        print(f"EC numbers {ec_numbers}", flush=True)

    id_list = []
    reaction_list = []
    ec_list = []
    idx = 0
    # loop over the reactions
    for i in range(n_reactions):
        if f_print:
            print(f"{i}: Reaction {reactions_set[i]}", flush=True)
        # Get the reaction
        reaction = reactions[reactions['id'] == reactions_set[i]]
        eq = reaction['reaction'].values[0]
        eq_re = eq
        # for each key value in the dictionary of cid
        if f_print:
            print(f"input eq:       {eq}", flush=True)
        # Loop over the data_bad_id items
        for k in range(n_compounds):
            val = data_bad_id[k]
            if val in eq:
                if f_print:
                    print(f"{k} Replacing compound {val}", flush=True)
                for j in range(len(cids_dict[val])):
                    eq_re = eq_re.replace(val, cids_dict[val][j])
                    new_id = make_custom_id(idx_base + idx, prefix=prefix)
                    if f_print:
                        print(f"{j}: Replacing Halogen {val} with {cids_dict[val][j]}", flush=True)
                        print(f"output eq_re:   {eq_re}", flush=True)
                        print(f"New ID {new_id}", flush=True)
                    idx += 1
                    # # Create an entry for the reaction id
                    id_list.append(new_id)
                    reaction_list.append(eq_re)
                    ec_list.append(ec_numbers[i])
    # Create a dataframe
    df = pd.DataFrame(data={
        "id": id_list,
        "reaction": reaction_list,
        "ec": ec_list})
    if int_file is not None:
        df.to_csv(int_file, compression='zip', encoding='utf-8', index=False)
    # Load the reactions data
    df_old = pd.read_csv(r_id_file, compression='zip')
    # Merge the dataframes
    df = pd.concat([df_old, df], ignore_index=True)
    # Drop the duplicates
    df = df.drop_duplicates(subset="id")
    # Sort the dataframe by the ID
    df = df.sort_values(by="id")
    # Save the dataframe
    df.to_csv(out_file, compression='zip', encoding='utf-8')
    return None
