import csv
import os

import pandas as pd
from rdkit import Chem as Chem
from rdkit import RDLogger

from .tools_mols import compound_super_safe_load, get_properties

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

from .tools_files import file_list_all, delete_files_substring
from .tools_eq import standardise_eq
from .tools_mp import tp_calc, mp_calc


def load_csv_to_dict(file_path):
    """
    Loads a CSV file and converts it to a dictionary.

    Parameters:
    file_path (str): The path to the CSV file.

    Returns:
    dict: A dictionary where the keys are the first column values and the values are the second column values.
    """
    result_dict = {}
    try:
        with open(file_path, mode='r', encoding='utf-8') as file:
            reader = csv.reader(file)
            for row in reader:
                if len(row) >= 2:  # Ensure there are at least two columns
                    key, value = row[0].strip(), row[1].strip()
                    result_dict[key] = value
    except FileNotFoundError:
        print(f"File not found: {file_path}, using an empty dictionary instead.", flush=True)

    return result_dict


def preprocess_kegg_c(target_dir, man_dict, outfile="kegg_data_C.csv.zip"):
    """
    Preprocesses KEGG compound data and saves it to a specified output file.

    Parameters:
    target_dir (str): The directory containing the KEGG compound data files.
    man_dict (dict): A dictionary of manual compound ID fixes.
    outfile (str, optional): The output file path for the preprocessed data. Defaults to "kegg_data_C.csv.zip".

    Returns:
    None
    """
    # Clean up the files
    delete_files_substring(target_dir, "_r")
    delete_files_substring(target_dir, "_p")
    # Get a list of all mol files in the directory
    files = [f for f in file_list_all(target_dir) if f.endswith('.mol')]

    print(f"Number of files: {len(files)}", flush=True)
    # Get the compound IDs
    arr_cid = [os.path.basename(file).split(".")[0] for file in files]
    # Load the mols in parallel
    mols = list(tp_calc(compound_super_safe_load, files))
    print(f"Number of mols: {len(mols)}", flush=True)
    # Get the ids of the failed mols
    ids_x = [cid for cid, mol in zip(arr_cid, mols) if mol is None]
    print(f"Number of X group files removed: {len(ids_x)}", flush=True)
    print(f"X group compounds {ids_x}", flush=True)
    # Remove the ids_x from the arr_cid
    arr_cid = [cid for cid in arr_cid if cid not in ids_x]
    # Remove the mols that are None
    mols = [mol for mol in mols if mol is not None]

    # Load the manual fixes
    mols += [Chem.MolFromSmiles(smiles) for cid, smiles in man_dict.items()]
    arr_cid += [cid for cid, smiles in man_dict.items()]

    # Get the properties
    properties = mp_calc(get_properties, mols)

    # Unpack the properties into the arrays
    arr_smiles, arr_formula, arr_mw, arr_n_heavy, arr_nc = zip(*properties)

    # Create a dataframe
    df = pd.DataFrame(data={
        "compound_id": arr_cid,
        "smiles": arr_smiles,
        "formula": arr_formula,
        "molecular_weight": arr_mw,
        "n_heavy_atoms": arr_n_heavy,
        "n_chiral_centers": arr_nc})
    # Sort the dataframe by the compound ID
    df = df.sort_values(by="compound_id").reset_index(drop=True).drop_duplicates().rename_axis(None, axis=1)
    # Save the dataframe
    df.to_csv(outfile, compression='zip', encoding='utf-8', index=False)
    print('Compound structural-info path: ' + outfile, flush=True)
    return df


def preprocess_kegg_c_metadata(target_dir='../../data/kegg_data_C_full',
                               outfile='data/kegg_data_C_metadata.csv.zip'):
    """
    Preprocesses KEGG compound metadata and saves it to a specified output file.

    Parameters:
    target_dir (str, optional): The directory containing the KEGG compound metadata files. Defaults to '../../data/kegg_data_C_full'.
    outfile (str, optional): The output file path for the preprocessed metadata. Defaults to 'data/kegg_data_C_metadata.csv.zip'.

    Returns:
    pd.DataFrame: A DataFrame containing the preprocessed compound metadata.
    """
    outfile = os.path.abspath(outfile)
    target_dir = os.path.abspath(target_dir)
    print('Importing compound metadata...', flush=True)

    paths = [os.path.join(root, file) for root, _, files in os.walk(target_dir) for file in files if
             file.endswith('.data')]

    df = pd.DataFrame({
        os.path.basename(path).split(".")[0]:
            pd.read_fwf(path, colspecs=[(0, 12), (12, -1)], header=None, names=['id', 'line'])
            .dropna(subset=['line']).ffill().set_index('id')['line'].str.strip().groupby(level=0).apply('|'.join)
        for path in paths
    }).drop('///', errors='ignore').T
    tar_list = ['name', 'remark', 'comment', 'sequence', 'type', 'brite']
    df = df.set_axis(df.columns.str.strip().str.lower(), axis=1).loc[:, tar_list].sort_index()
    df['glycan_ids'] = df['remark'].fillna('').str.extractall(r'(G\d{5})').groupby(level=0).agg(' '.join).replace('',
                                                                                                                  float(
                                                                                                                      'nan'))
    df.drop(columns='remark', inplace=True)
    df = df.reset_index().rename(columns={'index': 'compound_id'}).rename_axis(None, axis=1)
    df.to_csv(outfile, compression='zip', encoding='utf-8', index=False)
    print(f'Compound metadata path: {outfile}', flush=True)
    return df


def preprocess_kegg_r(target_dir, outfile, rm_gly=True):
    """
    Preprocesses KEGG reaction data and saves it to a specified output file.

    Parameters:
    target_dir (str): The directory containing the KEGG reaction data files.
    outfile (str): The output file path for the preprocessed data.
    rm_gly (bool, optional): Whether to remove reactions with glycan IDs. Defaults to True.

    Returns:
    None
    """
    print('Importing reaction data...', flush=True)
    # Get a list of all files in the directory
    paths = [m for n in [[f'{i}/{k}' for k in j] for i, _, j in list(os.walk(target_dir))[1:]] for m in n]

    # Import reaction data
    df = pd.DataFrame({os.path.basename(path).split(".")[0]:  # for each reaction ID
                           pd.read_fwf(path, colspecs=[(0, 12), (12, -1)], header=None,
                                       names=['id', 'line'])  # read file
                      .ffill(axis=0).set_index('id')  # indented lines relate to last-appearing header
                           ['line'].str.strip().groupby(level=0).apply('  ;  '.join)
                       # combine all lines for each header
                       for path in paths}).drop('///').T  # indexes are reaction IDs; cols are info types
    df = df.set_axis(df.columns.str.strip().str.lower(), axis=1).drop(  # remove columns not needed currently
        ['reference', 'authors', 'journal', 'title', 'brite', 'definition'], axis=1)
    if rm_gly:
        # Remove reactions with glycan IDs mixed in. "remark" column tags their equivalent reactions.
        df = df.loc[df['equation'].str.count(r"(\bG\d{5}\b)") == 0]
    for col in df.columns.difference(['comment']):
        df[col] = df[col].str.replace('  ;  ', ' ')  # ensure comment field structure is retained for parsing
    df['comment'] = df['comment'].str.replace('  ;  ', ';')
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
    df.rename(columns={'dblinks': 'rhea', 'entry': 'overall'}, inplace=True)
    df['overall'] = df['overall'].replace('', float('nan'))
    df = (df.loc[:, df.count().sort_values(ascending=False).index].drop(columns=['enzyme', 'equation'])
          .reset_index().rename({'index': 'id'}, axis=1).rename_axis(None, axis=1))
    df.to_csv(outfile, compression='zip', encoding='utf-8', index=False)
    print('Reaction info path: ' + outfile, flush=True)
    return df


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
        # gets compound metadata
        df_meta = preprocess_kegg_c_metadata(target_dir + '_full',
                                             outfile=out_file.replace('.csv.zip', '_metadata.csv.zip'))
        # Defines a list of bad CIDs to skip
        df_main = preprocess_kegg_c(target_dir, man_dict, outfile=out_file)
        print("C preprocessing done", flush=True)
        return df_meta, df_main
    elif target == "R":
        df = preprocess_kegg_r(target_dir, out_file)
        print("R preprocessing done", flush=True)
        return df
