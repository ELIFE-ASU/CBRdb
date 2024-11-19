import os

import pandas as pd
from rdkit import Chem as Chem
from rdkit import RDLogger

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

from tools_mols import standardize_mol
from tools_eq import eq_to_dict


def load_bad_entries(bad_file, target_str="molless"):
    with open(bad_file, 'r') as file:
        return [line.split(',')[0].strip() for line in file if target_str in line]


def get_reactions_with_substring(reactions_df, substring):
    return reactions_df[reactions_df['reaction'].str.contains(substring, case=False, na=False)]


def fix_halogen_compounds(c_id_bad_file="../data/C_IDs_bad.dat",
                          target_dir_c=r"../../data/kegg_data_C",
                          ):
    # Load the bad compound IDs
    data_bad_id = load_bad_entries(os.path.abspath(c_id_bad_file), target_str="X group")
    print(f"bad files ID with halogens: {data_bad_id}")

    # Prepare the halogen list to expand over
    hal_exp = ['F', 'Cl', 'Br', 'I']

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
        tmp_file = os.path.join(os.path.abspath(target_dir_c), data_bad_id[i], data_bad_id[i] + ".mol")
        with open(os.path.abspath(tmp_file), 'r') as f:
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
            print(f"Compound {data_bad_id[i]} with halogen {hal}, idx {idx} -> {cid}", flush=True)
            cids_dict[data_bad_id[i]].append(cid)
            smis_dict[data_bad_id[i]].append(smi)
    return cids_dict, smis_dict


def fix_halogen_reactions(cids_dict, smis_dict):
    c_id_file = "../data/kegg_data_C.csv.zip"
    r_id_file = "../data/kegg_data_R.csv.zip"
    r_id_file = "../data/atlas_data_kegg_R.csv.zip"

    print(f"cids_dict {cids_dict}")
    print(f"smis_dict {smis_dict}")

    # Load the bad compound IDs
    data_bad_id = list(cids_dict.keys())

    # load the compounds data
    compounds = pd.read_csv(c_id_file, compression='zip')
    # # Get the compounds with the halogens
    # compounds = compounds[compounds['compound_id'].isin(data)]
    # print(compounds['compound_id'].values)
    target_dir = r"../../data/kegg_data_C"
    target_dir = os.path.abspath(target_dir)

    # Load the reactions data
    reactions = pd.read_csv(r_id_file, compression='zip')

    # Find the reactions with the halogens
    re_set = set()
    for i in range(len(data_bad_id)):
        # Get the equations
        equations = get_reactions_with_substring(reactions, data_bad_id[i])
        re_set.update(equations['id'].values)
    re_set = list(re_set)
    print(f"Reactions with the halogens {re_set}")

    # loop over the reactions
    for i in range(len(re_set)):
        print(f"Reaction {re_set[i]}")
        # Get the reaction
        reaction = reactions[reactions['id'] == re_set[i]]
        eq = reaction['reaction'].values[0]
        eq_re = eq
        # for each key value in the dictionary of cid
        print(f"input eq:       {eq}")
        # break the eq
        lhs, rhs = eq_to_dict(eq)
        # print(lhs, rhs)
        # Loop over the data_bad_id items
        for val in data_bad_id:
            print(f"Replacing {val}")
            for j in range(len(cids_dict[val])):
                print(f"Replacing {val} with {cids_dict[val][j]}")
                eq_re = eq_re.replace(val, cids_dict[val][j])
                print(f"output eq_re:   {eq_re}")


if __name__ == "__main__":
    print("Program started", flush=True)
    cids_dict, smis_dict = fix_halogen_compounds()
    fix_halogen_reactions(cids_dict, smis_dict)
    print("Program finished", flush=True)
