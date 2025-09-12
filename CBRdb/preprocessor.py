import csv
import os

import pandas as pd
from rdkit import Chem as Chem
from rdkit import RDLogger

from .tools_mols import compound_super_safe_load, get_properties
from .merge_data_sets import id_indexed

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

from .tools_files import file_list_all, delete_files_substring, reaction_csv, compound_csv
from .tools_eq import standardise_eq
from .tools_mp import tp_calc


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


def preprocess_kegg_c(target_dir, man_dict):
    """
    Preprocesses KEGG compound data and saves it to a specified output file.

    Parameters:
    target_dir (str): The directory containing the KEGG compound data files.
    man_dict (dict): A dictionary of manual compound ID fixes.
    outfile (str, optional): The output file path for the preprocessed data. Defaults to "kegg_data_C.csv".

    Returns:
    None
    """
    # Clean up the files
    delete_files_substring(target_dir, "_r.mol")
    delete_files_substring(target_dir, "_p.mol")
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
    df_dict = {"compound_id": arr_cid}
    # Add properties to the dictionary
    df_dict.update(get_properties(mols))

    # Create a dataframe
    df = pd.DataFrame(data=df_dict)
    # Sort the dataframe by the compound ID
    df = df.sort_values(by="compound_id").reset_index(drop=True).drop_duplicates().rename_axis(None, axis=1)
    print('Finished importing compound structures.', flush=True)
    return df


def preprocess_kegg_c_metadata(target_dir='../../data/kegg_data_C_full',
                               valid_cids=None):
    """
    Preprocesses KEGG compound metadata.

    Parameters:
    target_dir (str, optional): The directory containing the KEGG compound metadata files. Defaults to '../../data/kegg_data_C_full'.
    valid_cids (iterable, optional): If provided, a list of compound IDs to keep. Defaults to None (keep all entries).

    Returns:
    pd.DataFrame: A DataFrame containing the preprocessed compound metadata.
    """
    target_dir = os.path.abspath(target_dir)

    # List metadata files in target_dir
    paths = pd.Series({file.split('.')[0]: os.path.join(root, file)
                       for root, _, files in os.walk(target_dir)
                       for file in files if file.endswith('.data')})
    # If list is empty, notify user of problem
    if len(paths) == 0:
        raise FileNotFoundError(f'No .data files in {target_dir} or any subdirectories.')
    # If user provided a list of compound IDs, import only those.
    if valid_cids is not None and hasattr(valid_cids, '__iter__'):
        if len(paths.index.intersection(valid_cids)) > 0:
            paths.drop(paths.index.difference(valid_cids), inplace=True, errors='ignore')

    print(f'Importing metadata for {len(paths.index)} compounds...', flush=True)

    # Import raw metadata
    df = pd.DataFrame({
        cid: pd.read_fwf(path, colspecs=[(0, 12), (12, -1)], header=None, names=['id', 'line'])
        .dropna(subset=['line']).ffill().set_index('id')['line'].str.strip().groupby(level=0).apply('~'.join)
        for cid, path in paths.items()}).drop('///', errors='ignore').T
    df = df.set_axis(axis=1, labels=df.columns.str.strip().str.lower()).sort_index()

    print(f'Formatting metadata...', flush=True)

    # Remove structural data since we get these elsewhere
    structure_cols = ['atom', 'bond', 'original', 'repeat', 'bracket']
    df.drop(columns=structure_cols + ['entry'], errors='ignore', inplace=True)

    # Process individual columns
    if 'brite' in df.columns:
        df['brite_full'] = df['brite'].copy(deep=True)
    if 'type' in df.columns:
        df['type'] = df['type'].str.lower()
    if 'exact_mass' in df.columns:
        df['exact_mass'] = df['mol_weight'].astype(float)
    if 'mol_weight' in df.columns:
        df['mol_weight'] = df['mol_weight'].astype(float)
    if 'name' in df.columns:
        df['nickname'] = df['name'].fillna('').map(lambda x: x.split(';~')[0])

    # Ease cross-referencing of databases by making database-specific fields
    for col in df.columns.intersection(['remark', 'dblinks']):
        icol = df[col].dropna().str.split('~').explode()
        col_df = pd.pivot(icol.str.split(': ', expand=True), columns=0, values=1)
        df = df.join(col_df, how='left')
        df.drop(columns=col, inplace=True)

    # Ease cross-referencing of KEGG GLYCAN and KEGG DRUG databases specifically
    if 'Same as' in df.columns:
        # prepare to extract glycan + drug IDs from field
        df['glycan'] = df['Same as'].copy(deep=True)
        df['drug'] = df['Same as'].copy(deep=True)
        # Future-proofing: check for presence of other ID types
        RetainSameAs = df['Same as'].str.split().explode().str.fullmatch(r'(G\d{5}|D\d{5})')
        RetainSameAs = RetainSameAs.eq(False).groupby(level=0).any()
        df['Same as'] = df['Same as'].mask(~RetainSameAs)

    # Standardize and compress representation of refs to KEGG databases
    string_of_ids = lambda l: ' '.join(sorted(list(set(l)))) if len(l) > 0 else float('nan')
    col_pat_mapping = {'module': r'(M\d{5})', 'glycan': r'(G\d{5})', 'drug': r'(D\d{5})',
                       'pathway': r'(map\d{5})', 'network': r'nt\d{5}(?:\(G\d{5}\))?',
                       'enzyme': r'(\d+\.\d+\.\d+\.\d+)', 'reaction': r'(R\d{5})', 'brite': r'(br\d{5})'}
    for col, pat in col_pat_mapping.items():
        if col in df.columns:
            df[col] = df[col].str.findall(pat).map(string_of_ids, na_action='ignore')

    print('Formatting metadata labels...', flush=True)
    kegg_cols = ['mol_weight', 'exact_mass', 'brite_full', 'sequence', 'type', 'formula', 'gene', 'organism']
    col_name_mapping = {'index': 'compound_id'} | {i: 'kegg_' + i for i in kegg_cols}
    col_name_mapping = col_name_mapping | {i: 'kegg_' + i for i in col_pat_mapping.keys()}
    df = df.sort_index().reset_index().rename(columns=col_name_mapping).rename_axis(None, axis=1)
    df.columns = df.columns.str.replace('-', '_').str.replace(' ', '_')
    df.dropna(axis=1, how='all', inplace=True)

    print('Finished importing compound metadata.', flush=True)
    print(f'Columns: {list(df.columns)}', flush=True)

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
                       for path in paths}).drop('///',
                                                errors='ignore').T  # indexes are reaction IDs; cols are info types
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
    ko_defs.sort_index().str.strip().to_csv(outfile.replace('.csv', '_kodefs.csv'), encoding='utf-8', header=None)
    del ko_defs

    # Extract reaction attributes and linkages
    df['reaction'] = df['equation'].str.replace('(side 1', '').str.replace('(side 2', '').apply(
        standardise_eq)  # standardize reaction formatting
    df['ec'] = df['enzyme'].fillna(' ').str.split().map(
        lambda x: ' '.join(sorted(list(x))))  # combine all ECs, including partials

    patterns = {'orthology': r"(\bK\d{5}\b)", 'pathway': r"(\brn\d{5}\b)", 'module': r"(\bM\d{5}\b)",
                'rclass': r"(\bRC\d{5}\b  \bC\d{5}_C\d{5})", 'dblinks': r"( \d{5})", 'entry': 'Overall'}
    [df.update(df[k].str.findall(v).map(lambda x: ' '.join(sorted(list(x))), na_action='ignore')) for k, v in
     patterns.items()]
    df.loc[:, 'rclass'] = df.loc[:, 'rclass'].str.replace('  ', '__')

    # Rename columns where appropriate
    df.rename(columns={'dblinks': 'rhea', 'entry': 'overall'}, inplace=True)
    df['overall'] = df['overall'].replace('', float('nan'))
    df = (df.loc[:, df.count().sort_values(ascending=False).index].drop(columns=['enzyme', 'equation'])
          .reset_index().rename({'index': 'id'}, axis=1).rename_axis(None, axis=1)).sort_values(by='id')
    reaction_csv(df, outfile)
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
    out_file = os.path.abspath(f"{out_file}_{target}.csv")
    cid_manual_file = os.path.abspath(cid_manual_file)

    if target == "C":
        # Defines a dictionary of manual fixes whose mol files are not found in target_dir
        man_dict = load_csv_to_dict(cid_manual_file)
        # converts compound mol files to smiles strings; defines a list of CIDs to skip
        df_main = preprocess_kegg_c(target_dir, man_dict)
        # gets compound metadata e.g. names + classifications
        df_meta = preprocess_kegg_c_metadata(target_dir + '_full',
                                             valid_cids=list(df_main['compound_id'].sort_values()))
        # merges the compound data
        df = df_main.merge(df_meta, on='compound_id', how='outer').sort_values(by='compound_id').reset_index(drop=True)
        # tag compounds whose structures are missing
        missing = log_missing_structures(df)
        # generates output file with compound data
        compound_csv(df, out_file)
        print("C preprocessing done. Compound info path:" + out_file, flush=True)
        return df
    elif target == "R":
        df = preprocess_kegg_r(target_dir, out_file)
        print("R preprocessing done", flush=True)
        return df


def log_missing_structures(df):
    """
    Logs compounds with metadata, but without SMILES strings or mol files. Uses metadata to infer priority compounds.

    Parameters:
    df (pd.DataFrame): A DataFrame containing compound data with 'compound_id' column.

    Returns:
    pd.DataFrame: A DataFrame of compounds that are missing but promising for follow-up.
    """
    # Get the set of existing compound IDs
    extant = set(df['compound_id'])

    # Get the list of compound IDs from the metadata files that are not in the existing set
    k = [i.split('/')[-2] for i in file_list_all('../../data/kegg_data_C_full') if
         i.endswith('.data') & (i.split('/')[-2] not in extant)]

    # Preprocess the KEGG compound metadata for compounds missing structures
    missing = preprocess_kegg_c_metadata(valid_cids=sorted(k))
    # Omit generic halogen compounds as we enumerate these elsewhere
    missing = missing.query('~ kegg_formula.str.contains("X", na=False)').reset_index(drop=True)
    # Tag compounds prioritized for manual addition
    missing['priority'] = missing.index.isin(filter_missing_structures(missing).index)
    # Save compounds to a file
    compound_csv(df_C=missing, file_address='../data/C_IDs_good.dat')

    return missing


def filter_missing_structures(df):
    """
    Filters compounds for manual addition based on specific criteria.

    Parameters:
    df (pd.DataFrame): A DataFrame containing compound data with 'compound_id' column.

    Returns:
    pd.DataFrame: A DataFrame of compounds prioritized for manual addition in C_IDs_manual.dat
    """
    # Define the query to identify compounds prioritized for manual addition
    query_dict = {'kegg_sequence': 'kegg_sequence.isna()',
                    'kegg_type': 'kegg_type.isna()',
                    'remark': '~remark.str.count(r"Same as: 	G").eq(0)',
                    'kegg_reaction': 'kegg_reaction.notna()',
                    'name': '~name.str.lower().str.contains("protein|globin|doxin|glycan|lase|peptide|rna|dna|steroid|lipid|lignin", na=False)',
                    'comment': '~comment.fillna("").str.lower().str.contains("peptide|protein|[KO:", na=False, regex=False)',
                    'kegg_formula': '~kegg_formula.str.contains("X", na=False)',
                    'kegg_brite': '~kegg_brite.str.contains("08009|08005", na=False)'}

    # Combine all eligible queries based on the fields present
    combined_query = ' & '.join([v for k, v in query_dict.items() if k in df.columns])

    # Apply the query to identify promising compounds
    priority_compounds = df.query(combined_query)

    return priority_compounds
