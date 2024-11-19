import os

import pandas as pd
from rdkit import Chem as Chem
from rdkit import RDLogger

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

from .tools_mols import standardize_mol
from .tools_eq import eq_to_dict


def load_bad_entries(bad_file, target_str="molless"):
    with open(bad_file, 'r') as file:
        return [line.split(',')[0].strip() for line in file if target_str in line]


def get_reactions_with_substring(reactions_df, substring):
    return reactions_df[reactions_df['reaction'].str.contains(substring, case=False, na=False)]


def fix_halogen_compounds(c_id_bad_file="../data/C_IDs_bad.dat",
                          target_dir_C=r"../../data/kegg_data_C",
                          ):
    # Load the bad compound IDs
    data_bad_id = load_bad_entries(os.path.abspath(c_id_bad_file), target_str="X group")
    print(f"bad files ID: {data_bad_id}")

    # Prepare the halogen list to expand over
    hal_exp = ['F', 'Cl', 'Br', 'I']

    # Xe is not a halogen this ID is not a halogen and should be removed
    try:
        data_bad_id.remove("C13373")
    except ValueError:
        pass
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
        tmp_file = os.path.join(os.path.abspath(target_dir_C), data_bad_id[i], data_bad_id[i] + ".mol")
        with open(os.path.abspath(tmp_file), 'r') as f:
            file_data = f.read()
        # Initialize the list for the current data[i]
        cids_dict[data_bad_id[i]] = []
        smis_dict[data_bad_id[i]] = []
        # replace the X with the halogen
        for j, hal in enumerate(hal_exp):
            idx = num_range[i][j]
            mol = Chem.MolFromMolBlock(file_data.replace("X", hal))
            mol = standardize_mol(mol)
            smi = Chem.MolToSmiles(mol, allHsExplicit=True)
            # Determine the compound id from the ones already given
            cid = {
                "C00462": {"F": "C16487", "Cl": "C01327", "Br": "C13645", "I": "C05590"},
                "C01812": {"F": "C06108", "Cl": "C06755"}
            }.get(data_bad_id[i], {}).get(hal)

            # Generate the compound id if not found
            if cid is None:
                cid = f"C{int(99000 + idx):05d}"
            # print(f"Compound {data[i]} with halogen {hal}, idx {idx} -> {cid}")
            cids_dict[data_bad_id[i]].append(cid)
            smis_dict[data_bad_id[i]].append(smi)
    return cids_dict, smis_dict


def fix_halogen_reactions():
    c_id_file = "../data/kegg_data_C.csv.zip"
    r_id_file = "../data/kegg_data_R.csv.zip"
    r_id_file = "../data/atlas_data_kegg_R.csv.zip"

    cids_dict, smis_dict = fix_halogen_compounds()
    print(f"cids_dict {cids_dict}")
    print(f"smis_dict {smis_dict}")

    # Load the bad compound IDs
    # data_bad_id = load_bad_entries(os.path.abspath(c_id_bad_file), target_str="X group")

    # load the compounds data
    compounds = pd.read_csv(c_id_file, compression='zip')
    # # Get the compounds with the halogens
    # compounds = compounds[compounds['compound_id'].isin(data)]
    # print(compounds['compound_id'].values)
    target_dir = r"../../data/kegg_data_C"
    target_dir = os.path.abspath(target_dir)
    # Load the data
    reactions = pd.read_csv(r_id_file, compression='zip')
    re_set = set()
    for i in range(len(data)):
        # Get the equations
        equations = get_reactions_with_substring(reactions, data[i])
        re_set.update(equations['id'].values)
    re_set = list(re_set)
    print(f"Reactions with the halogens {re_set}")
    print(cids_dict)
    # loop over the reactions
    for i in range(len(re_set)):
        # Get the reaction
        reaction = reactions[reactions['id'] == re_set[i]]
        eq = reaction['reaction'].values
        # for each key value in the dictionary of cid
        print(eq)
        # break the eq
        lhs, rhs = eq_to_dict(eq[0])
        print(lhs, rhs)
        keys_list = list(lhs.keys())
        values_list = list(lhs.values())

        for i, key in enumerate(keys_list):
            if key in data:
                for j, val in enumerate(cids_dict[key]):
                    # # Get the index of the key
                    # idx = data.index(key)
                    # # Get the new value
                    new_value = cids_dict[key][j]
                    print(f"Replacing {key} with {new_value}")
                    # construct the new equation
                    lhs[new_value] = lhs.pop(key)
        #         print(f"Replacing {key} with {new_value}")
        #         # Update the lhs
        #         lhs[new_value] = lhs.pop(key)

        # # loop over the lhs
        # for key, value in lhs.items():
        #     for j, hal in enumerate(hal_exp):
        #         if key in data:
        #             # Get the index of the key
        #             idx = data.index(key)
        #             # Get the new value
        #             new_value = cids_dict[key][j]
        #             print(f"Replacing {key} with {new_value}")
        #             # Update the lhs
        #             lhs[new_value] = lhs.pop(key)
        #             # # Update the lhs
        #             # lhs[new_value] = lhs.pop(key)
        #     # if key in data:
        #     #     # Get the index of the key
        #     #     idx = data.index(key)
        #     #     # Get the new value
        #     #     new_value = cids_dict[key][hal_exp.index(value)]
        #     #     # Update the lhs
        #     #     lhs[new_value] = lhs.pop(key)


if __name__ == "__main__":
    print("Program started", flush=True)
    fix_halogen_compounds()
    fix_halogen_reactions()
    print("Program finished", flush=True)
