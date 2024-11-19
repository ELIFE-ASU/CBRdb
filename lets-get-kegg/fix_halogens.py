import os

import pandas as pd
from rdkit import Chem as Chem
from rdkit import RDLogger

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

from tools_mols import standardize_mol, get_mol_descriptors
from tools_eq import eq_to_dict


def load_bad_entries(bad_file, target_str="molless"):
    with open(bad_file, 'r') as file:
        return [line.split(',')[0].strip() for line in file if target_str in line]


def get_reactions_with_substring(reactions_df, substring):
    return reactions_df[reactions_df['reaction'].str.contains(substring, case=False, na=False)]


def fix_halogen_compounds(c_id_bad_file="../data/C_IDs_bad.dat",
                          target_dir_c=r"../../data/kegg_data_C",
                          hal_exp=None,
                          ):
    # Prepare the full path of the files
    c_id_bad_file = os.path.abspath(c_id_bad_file)
    target_dir_c = os.path.abspath(target_dir_c)
    # Prepare the halogen list to expand over
    if hal_exp is None:
        hal_exp = ['F', 'Cl', 'Br', 'I']
    # Load the bad compound IDs
    data_bad_id = load_bad_entries(c_id_bad_file, target_str="X group")
    print(f"bad files ID with halogens: {data_bad_id}", flush=True)

    # Make the combinations of the data and the halogens
    n_data = len(data_bad_id)
    n_hal = len(hal_exp)
    n_comb = n_data * n_hal
    num_range = range(n_comb)
    # Reshape num_range to a 2D array
    num_range = [num_range[i:i + n_hal] for i in range(0, n_comb, n_hal)]
    smis_dict = {}
    cids_dict = {}
    for i in range(n_data):
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
            print(f"Compound {data_bad_id[i]} with halogen {hal}, idx {idx} -> {smi}", flush=True)
            # Determine the compound id from the ones already given
            cid = {
                "C00462": {"F": "C16487", "Cl": "C01327", "Br": "C13645", "I": "C05590"},
                "C01322": {},  # full expansion and C added required!
                "C01365": {},  # No reaction data!
                "C01706": {},  # No reaction data!
                "C01812": {"F": "C06108", "Cl": "C06755"},
                "C01813": {},  # No reaction data!
                "C01872": {},  # full expansion and C added required!
                "C02103": {},  # full expansion and C added required!
                "C03122": {},  # No reaction data!
                "C15564": {},  # full expansion and C added required!

            }.get(data_bad_id[i], {}).get(hal)

            # Generate the compound id if not found
            if cid is None:
                cid = f"C{int(99000 + idx):05d}"
            cids_dict[data_bad_id[i]].append(cid)
            smis_dict[data_bad_id[i]].append(smi)
    return cids_dict, smis_dict


def merge_halogen_compounds(cids_dict, smis_dict, c_id_file="../data/kegg_data_C.csv.zip", out_file=None):
    if out_file is None:
        out_file = c_id_file
    # Prepare the full path of the files
    c_id_file = os.path.abspath(c_id_file)

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

    # Load the compounds data
    df_old = pd.read_csv(c_id_file, compression='zip', index_col=0)

    # Merge the dataframes
    df = pd.concat([df_old, df], ignore_index=True)
    # Drop the duplicates
    df = df.drop_duplicates(subset="compound_id")
    # Sort the dataframe by the compound ID
    df = df.sort_values(by="compound_id")
    # Save the dataframe
    df.to_csv(out_file, compression='zip', encoding='utf-8')
    return None


def fix_halogen_reactions(cids_dict, smis_dict, r_id_file = "../data/kegg_data_R.csv.zip"):
    # Prepare the full path of the files
    r_id_file = os.path.abspath(r_id_file)

    print(f"cids_dict {cids_dict}")
    print(f"smis_dict {smis_dict}")

    # Load the bad compound IDs
    data_bad_id = list(cids_dict.keys())

    # Load the reactions data
    reactions = pd.read_csv(r_id_file, compression='zip')

    # Find the reactions with the halogens
    re_set = set()
    for i in range(len(data_bad_id)):
        # Get the equations
        equations = get_reactions_with_substring(reactions, data_bad_id[i])
        re_set.update(equations['id'].values)
    re_set = list(re_set)
    print(f"Reactions with the halogens {re_set}", flush=True)

    # loop over the reactions
    for i in range(len(re_set)):
        print(f"Reaction {re_set[i]}", flush=True)
        # Get the reaction
        reaction = reactions[reactions['id'] == re_set[i]]
        eq = reaction['reaction'].values[0]
        eq_re = eq
        # for each key value in the dictionary of cid
        print(f"input eq:       {eq}", flush=True)
        # break the eq
        lhs, rhs = eq_to_dict(eq)
        # print(lhs, rhs)
        # Loop over the data_bad_id items
        for val in data_bad_id:
            print(f"Replacing {val}", flush=True)
            for j in range(len(cids_dict[val])):
                print(f"Replacing {val} with {cids_dict[val][j]}", flush=True)
                eq_re = eq_re.replace(val, cids_dict[val][j])
                print(f"output eq_re:   {eq_re}", flush=True)


if __name__ == "__main__":
    print("Program started", flush=True)
    cids_dict, smis_dict = fix_halogen_compounds()
    merge_halogen_compounds(cids_dict, smis_dict)
    # fix_halogen_reactions(cids_dict, smis_dict)
    print("Program finished", flush=True)
