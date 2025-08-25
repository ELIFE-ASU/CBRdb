import os

import pandas as pd
from rdkit import Chem as Chem
from rdkit import RDLogger

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

from .tools_mols import standardize_mol, check_for_x_group, get_properties, check_for_r_group, replace_r_group
from .tools_files import make_custom_id, reaction_csv
from .tools_mp import tp_calc
from .tools_eq import convert_formula_to_dict, standardise_eq


def load_bad_entries(target_dir_c):
    """
    Loads and filters bad entries from a target directory.

    Parameters:
    target_dir_c (str): The path to the target directory containing the files.

    Returns:
    list: A list of filenames that are flagged as bad entries.
    """
    target_dir_c = os.path.abspath(target_dir_c)
    files = sorted([f for f in os.listdir(target_dir_c) if not f.startswith('.')])
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
    tuple: A tuple containing two dictionaries and a DataFrame:
        - cids_dict (dict): A dictionary mapping original compound IDs to new compound IDs.
        - smis_dict (dict): A dictionary mapping original compound IDs to their SMILES representations.
        - specific_halogens (DataFrame): A DataFrame summarizing the above halogen compounds.
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
    idx = 0
    smis_dict = {}
    cids_dict = {}
    for i in range(len(data_bad_id)):
        # Get the full path of the file
        tmp_file = os.path.join(target_dir_c, data_bad_id[i], data_bad_id[i] + ".mol")

        # Check for and standardize R groups. as in compound_super_safe_load
        flag_r = check_for_r_group(tmp_file)
        if flag_r:
            f_load_r = tmp_file.split(".")[0] + "_r.mol"
            replace_r_group(tmp_file, f_load_r)
            tmp_file = f_load_r

        # Open the file
        with open(tmp_file, 'r') as f:
            file_data = f.read()

        # Remove temporary file generated if R groups were found
        if flag_r:
            os.remove(tmp_file)
        
        # Initialize the list for the current data[i]
        cids_dict[data_bad_id[i]] = []
        smis_dict[data_bad_id[i]] = []
        # Replace the X with the halogen
        # First, use a sample halogen F, since these load correctly + no X+F combos exist
        file_data = file_data.replace("X", "F")
        # Create the molecule; use same read-in specs as compound_super_safe_load
        base_mol = Chem.MolFromMolBlock(file_data, sanitize=False, removeHs=False)
        # Standardize that molecule
        base_mol = standardize_mol(base_mol)
        # Iterate through halogens
        for j, hal in enumerate(hal_exp):
            mol = Chem.ReplaceSubstructs(base_mol,
                                          Chem.MolFromSmiles('F'),
                                          Chem.MolFromSmiles(hal),
                                          replaceAll=True)[0]
            # Standardize the molecule
            mol = standardize_mol(mol)
            smi = Chem.MolToSmiles(mol, allHsExplicit=True)

            # Ensure R group's presence is registered by Chem.MolFromSmiles in merger step
            if flag_r:
                smi = smi.replace('[H:0]', '*')

            # Determine the compound id from the ones already given
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
                idx += 1
            if f_print:
                print(f"Compound {data_bad_id[i]} with halogen {hal} -> {cid}: {smi}", flush=True)
            # Add the data to the dictionaries
            cids_dict[data_bad_id[i]].append(cid)
            smis_dict[data_bad_id[i]].append(smi)
    if f_print:
        print(f"cids_dict {cids_dict}", flush=True)
        print(f"Added {idx} new compounds", flush=True)
    specific_halogens = (pd.DataFrame({'compound_id': cids_dict, 'smiles': smis_dict})
                         .explode(['compound_id', 'smiles'])
                         .assign(elem=lambda x: x['smiles'].str.findall('F|Cl|Br|I').explode(),
                                 is_new=lambda x: x['compound_id'].str.startswith('C9'))
                         .reset_index(names='generic_id'))
    return cids_dict, smis_dict, specific_halogens


def merge_halogen_compounds(cids_dict,
                            smis_dict,
                            c_id_file="../data/kegg_data_C.csv",
                            int_file=None,
                            out_file=None):
    """
    Merges halogen compound data with existing compound data, removing duplicates and saving the result.

    Parameters:
    cids_dict (dict): A dictionary mapping original compound IDs to new compound IDs.
    smis_dict (dict): A dictionary mapping original compound IDs to their SMILES representations.
    c_id_file (str): The path to the file containing existing compound data. Default is '../data/kegg_data_C.csv'.
    int_file (str, optional): The path to the intermediate file where merged data will be saved. Default is None.
    out_file (str, optional): The path to the output file where final merged data will be saved. Default is the same as c_id_file.

    Returns:
    None
    """
    print("Merging halogen compounds", flush=True)
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

    mols = [Chem.MolFromSmiles(smi) for smi in smis_list]
    # Get the properties
    df_dict = {"compound_id": cids_list}
    # Add properties to the dictionary
    df_dict.update(get_properties(mols))

    # Create a dataframe
    df = pd.DataFrame(data=df_dict)
    if int_file is not None:
        df.to_csv(int_file, encoding='utf-8', index=False, float_format='%.3f')

    # Load the compounds data
    df_old = pd.read_csv(c_id_file)

    # Merge the dataframes
    df = pd.concat([df_old, df], ignore_index=True)
    # Drop the duplicates
    df = df.drop_duplicates(subset="compound_id")
    # Sort the dataframe by the compound ID
    df = df.sort_values(by="compound_id")
    # Save the dataframe
    df.to_csv(out_file, index=False, encoding='utf-8', float_format='%.3f')
    print("Done!", flush=True)
    return df


def merge_halogen_compounds_pd(C_main, specific_halogens, out_file="../data/kegg_data_C.csv", int_file=None):
    """
    Merges halogen compound data with existing compound data, saving the result.

    Parameters:
    C_main (pd.DataFrame): A pandas DataFrame containing existing compound data.
    specific_halogens (pd.DataFrame): A pandas DataFrame containing specific halogen compounds extrapolated from generic ones, from fix_halogen_compounds.
    out_file (str, optional): File name to which output csv should be written. Defaults to "../data/kegg_data_C.csv".
    int_file (str, optional): Path to the intermediate file where new halogen compounds will be saved. Default is None.

    Returns:
    pd.DataFrame: The updated DataFrame with new halogen compounds merged in.
    """
    # Filter out existing halogens
    new_halogens = specific_halogens.query('is_new').set_index('compound_id')
    # Calculate molecular properties
    property_names = ['smiles',
                      'smiles_capped',
                      'inchi_capped',
                      'formula',
                      'molecular_weight',
                      'n_heavy_atoms',
                      'n_chiral_centers']
    property_values = new_halogens['smiles'].map(Chem.MolFromSmiles).map(get_properties)
    halogen_properties = property_values.map(lambda x: dict(zip(property_names, x))).apply(pd.Series).reset_index()
    halogen_properties['smiles'] = new_halogens['smiles'].reset_index(drop=True)
    # Merge into compound database
    C_main = pd.concat([C_main, halogen_properties], ignore_index=True)
    # Write output file(s)
    if int_file is not None:
        halogen_properties.to_csv(int_file, index=False, encoding='utf-8', float_format='%.3f')
    C_main.sort_values(by='compound_id').to_csv(out_file, index=False, encoding='utf-8', float_format='%.3f')
    return C_main


def fix_halogen_reactions(cids_dict,
                          r_id_file="../data/atlas_data_R.csv",
                          int_file=None,
                          out_file=None,
                          f_print=True):
    """
    Fixes halogen reactions by replacing halogen placeholders in reaction equations and generating new reaction IDs.

    Parameters:
    cids_dict (dict): A dictionary mapping original compound IDs to new compound IDs.
    r_id_file (str): The path to the file containing existing reaction data. Default is '../data/atlas_data_R.csv'.
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
    data_bad_id = sorted(list(cids_dict.keys()))
    n_compounds = len(data_bad_id)
    if f_print:
        print(f"{n_compounds} bad compound files ID with halogens: {data_bad_id}", flush=True)

    # Load the reactions data
    df_reactions = pd.read_csv(r_id_file)
    if f_print:
        print(df_reactions.columns, flush=True)

    # Find the reactions with the halogens
    reactions_halogen_set = set()
    for i in range(len(data_bad_id)):
        # Get the equations
        equations = get_reactions_with_substring(df_reactions, data_bad_id[i])
        reactions_halogen_set.update(equations['id'].values)
    reactions_halogen_set = sorted(list(reactions_halogen_set))
    n_reactions = len(reactions_halogen_set)
    if f_print:
        print(f"Found {n_reactions} reactions with the halogens, IDs: {reactions_halogen_set}", flush=True)

    # Only keep the reactions with the halogens in the reactions_set
    df_halogens = df_reactions[df_reactions['id'].isin(reactions_halogen_set)]
    print(df_halogens, flush=True)

    df_halogens_exp = pd.DataFrame()

    idx = 0
    # loop over the reactions
    for i in range(n_reactions):
        if f_print:
            print(f"{i}: Reaction {reactions_halogen_set[i]}", flush=True)
        # Get the reaction
        reaction = df_halogens[df_halogens['id'] == reactions_halogen_set[i]]
        if f_print:
            print(reaction, flush=True)

        eq = reaction['reaction'].values[0]
        eq_re = eq
        if f_print:
            print(f"input eq:       {eq}", flush=True)
        # Loop over the data_bad_id items for each key value in the dictionary of cid
        for k in range(n_compounds):
            val = data_bad_id[k]
            # Check if the value is in the equation
            if val in eq:
                if f_print:
                    print(f"{k}: Replacing compound {val}", flush=True)
                for j in range(len(cids_dict[val])):
                    eq_re = eq.replace(val, cids_dict[val][j])
                    new_id = make_custom_id(idx_base + idx, prefix=prefix)
                    if f_print:
                        print(f"{j}: Replacing Halogen {val} with {cids_dict[val][j]}", flush=True)
                        print(f"output eq_re:   {eq_re}", flush=True)
                        print(f"New ID {new_id}", flush=True)
                    idx += 1
                    # Add the data to the dataframe
                    df_copy = reaction.copy()
                    df_copy['id'] = new_id
                    df_copy['reaction'] = eq_re
                    df_halogens_exp = df_halogens_exp._append(df_copy, ignore_index=True)
    print(f"Added {idx} new reactions", flush=True)

    if int_file is not None:
        print(f"Saving intermediate file {int_file}", flush=True)
        reaction_csv(df_halogens_exp, int_file)

    # Merge the dataframes
    df = pd.concat([df_reactions, df_halogens_exp], ignore_index=True)
    # Drop the duplicates
    df = df.drop_duplicates(subset="id")
    # Remove the old reactions
    df = df[~df['id'].isin(reactions_halogen_set)]
    # Sort the dataframe by the ID
    df = df.sort_values(by="id").reset_index(drop=True)
    # Save the dataframe
    reaction_csv(df, out_file)
    return df


def fix_halogen_reactions_without_existing_halogens(df_R, C_main, specific_halogens):
    """
    Fixes halogen reactions by replacing generic halogen placeholders with specific halogen compounds.

    Parameters:
    df_R (pd.DataFrame): A pandas DataFrame containing reaction data.
    C_main (pd.DataFrame): A pandas DataFrame containing main compound data.
    specific_halogens (pd.DataFrame): A DataFrame summarizing specific halogen compounds.

    Returns:
    pd.DataFrame: The updated DataFrame with fixed halogen reactions.
    """
    # make sure all reactions have same ID convention
    if df_R['id'].str.startswith("R").all():
        prefix, idx_base = 'R', 99000
    elif df_R['id'].str.startswith("A").all():
        prefix, idx_base = 'A', 990000
    else:
        raise ValueError("multiple ID formats found; standardize ID format before attempting to assign new IDs")

    halogen_elements = ['F', 'Cl', 'Br', 'I']
    halogen_cps = (C_main.set_index('compound_id')['formula']
                   .map(lambda x: ' '.join(set(convert_formula_to_dict(x).keys()).intersection(halogen_elements)))
                   .replace('', None)).dropna().str.split()

    # subset reactions with generic halogens
    halogen_reactions = df_R.loc[df_R['reaction'].str.contains("|".join(specific_halogens['generic_id']))].copy(
        deep=True)
    halogen_reactions['generics_used'] = halogen_reactions['reaction'].str.findall(
        "|".join(specific_halogens['generic_id']) + '|<=>').map(' '.join)
    halogen_reactions['generic_both_sides'] = halogen_reactions['generics_used'].str.split('<=>').map(len) == 2
    halogen_reactions['no_extant_found'] = halogen_reactions['reaction'].str.findall("|".join(halogen_cps.index)).map(
        len) == 0

    # alert user if existing halogen compounds are found - element(s) of these compounds should inform element(s) of specific compounds.
    if not halogen_reactions['no_extant_found'].all():
        print('existing halogen compounds found in: ',
              ' '.join(halogen_reactions[~halogen_reactions['no_extant_found']]['id']))

    # most of the reactions with generic halogens have only generics. enumerate these by element, keeping element consistent within a given reaction.
    halogen_grps = halogen_reactions.query('generic_both_sides & no_extant_found').assign(
        generics_used=lambda x: x.generics_used.str.findall(r'C\d{5}'))
    halogen_defs = specific_halogens.set_index(['generic_id', 'elem'])['compound_id']
    halogen_grps['specifics_used'] = halogen_grps['generics_used'].map(
        lambda x: [dict(zip(x, halogen_defs.loc[x, el])) for el in halogen_elements])
    halogen_grps = halogen_grps.explode('specifics_used').reset_index(drop=True)

    # replace generics with specifics.
    for i, row in halogen_grps.iterrows():
        for gen, spec in row['specifics_used'].items():
            halogen_grps.at[i, 'reaction'] = halogen_grps.at[i, 'reaction'].replace(gen, spec)

    # standardize reaction format and remove reactions already found elsewhere in the reaction database.
    halogen_grps['reaction'] = halogen_grps['reaction'].map(standardise_eq)
    halogen_grps = halogen_grps.query('~reaction.isin(@df_R.reaction)').reset_index(drop=True)

    # remove generic halogen reactions from the database
    df_R = df_R.query('~id.isin(@halogen_grps.id)').reset_index(drop=True)

    # assign new IDs to generic halogen reactions
    halogen_grps['id'] = prefix + (halogen_grps.index + idx_base).astype(str)

    # append these new reactions to the existing reaction database.
    df_R = pd.concat([df_R, halogen_grps[df_R.columns]], ignore_index=True).sort_values(by='id').reset_index(drop=True)

    return df_R
