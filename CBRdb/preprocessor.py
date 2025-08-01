import csv
import os

import pandas as pd
from rdkit import Chem as Chem
from rdkit import RDLogger

from .tools_mols import compound_super_safe_load, get_properties

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

from .tools_files import file_list_all, delete_files_substring, reaction_csv
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
                               valid_cids=None,
                               keep_only=None):
    """
    Preprocesses KEGG compound metadata.

    Parameters:
    target_dir (str, optional): The directory containing the KEGG compound metadata files. Defaults to '../../data/kegg_data_C_full'.
    valid_cids (iterable, optional): If provided, a list of compound IDs to keep. Defaults to None (keep all entries).
    keep_only (list, optional): A list of metadata fields to keep. Defaults to None; if none, apply standard filtering.

    Returns:
    pd.DataFrame: A DataFrame containing the preprocessed compound metadata.
    """
    target_dir = os.path.abspath(target_dir)

    if valid_cids is not None:
        print(f'Importing metadata for {len(valid_cids)} compounds...', flush=True)
        paths = [os.path.join(root, file) for root, _, files in os.walk(target_dir) for file in files if
                 file.endswith('.data') and file.split('.')[0] in valid_cids]
    else:
        print(f'Importing compound metadata...', flush=True)
        paths = [os.path.join(root, file) for root, _, files in os.walk(target_dir) for file in files if
                 file.endswith('.data')]

    df = pd.DataFrame({
        os.path.basename(path).split(".")[0]:
            pd.read_fwf(path, colspecs=[(0, 12), (12, -1)], header=None, names=['id', 'line'])
            .dropna(subset=['line']).ffill().set_index('id')['line'].str.strip().groupby(level=0).apply('~'.join)
        for path in paths}).drop('///', errors='ignore').T
    df = df.set_axis(df.columns.str.strip().str.lower(), axis=1).sort_index()

    if keep_only is not None:
        if 'brite' in df.columns and 'brite' in keep_only:
            df['brite'] = df['brite'].str.lower().str.findall('protein|peptide|enzyme').fillna('').map(
                lambda x: ' '.join(sorted(list(set(x)))))
            if 'type' in df.columns:
                try:
                    df['type'] = (df['type'].fillna('').str.lower() + ' ' + df['brite']).str.strip().str.split().map(
                        lambda x: ' '.join(sorted(list(set(x)))))
                except:
                    df['type'] = df['type'].notna() or df['brite'].notna()
        if set(keep_only).issubset(set(df.columns)):
            return df.drop(columns=df.columns.difference(keep_only), errors='ignore')
        elif len(df.columns.intersection(keep_only)) > 0:
            print(f'Requested field not found: {set(keep_only).difference(set(df.columns))}', flush=True)
            print(f'Keeping other fields: {set(keep_only).intersection(set(df.columns))}', flush=True)
            return df.drop(columns=set(df.columns).difference(keep_only), errors='ignore')
        else:
            print(f'No requested fields found. Defaulting to keeping all metadata.', flush=True)
            return df

    if 'remark' in df.columns:
        df['glycan_ids'] = (df['remark'].fillna('').str.extractall(r'(G\d{5})')
                            .groupby(level=0).agg(' '.join).replace('', float('nan')))
        df['drug_ids'] = (df['remark'].fillna('').str.extractall(r'(D\d{5})')
                          .groupby(level=0).agg(' '.join).replace('', float('nan')))

    df = df.sort_index().reset_index().rename(
        columns={'index': 'compound_id', 'formula': 'kegg_formula', 'mol_weight': 'kegg_mol_weight'}).rename_axis(None,
                                                                                                                  axis=1)
    default_keep_cols = 'compound_id,comment,dblinks,kegg_mol_weight,kegg_formula,name,glycan_ids,drug_ids'.split(',')
    df.drop(columns=df.columns.difference(default_keep_cols), inplace=True, errors='ignore')
    df['kegg_mol_weight'] = df['kegg_mol_weight'].astype(float)
    df['nickname'] = df['name'].fillna('').map(lambda x: x.split(';~')[0])

    if valid_cids is not None:
        if hasattr(valid_cids, '__iter__') and len(set(df['compound_id'].values).intersection(valid_cids)) > 0:
            df = df.query('compound_id.isin(@valid_cids)')
        else:
            print(f'valid CIDs must be iterable and overlapping with IDs in metadata. Keeping all entries.', flush=True)
    print('Finished importing compound metadata.', flush=True)
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
        # generates output file with compound data
        df.to_csv(out_file, encoding='utf-8', index=False, float_format='%.3f')
        print("C preprocessing done. Compound info path:" + out_file, flush=True)
        return df
    elif target == "R":
        df = preprocess_kegg_r(target_dir, out_file)
        print("R preprocessing done", flush=True)
        return df


def log_compounds_for_followup(df):
    """
    Logs compounds without SMILES strings or mol files, whose metadata suggests we might want to seek structures elsewhere.

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

    # Define the list of metadata fields to keep
    keep_only = ['name', 'sequence', 'type', 'formula', 'remark', 'reaction', 'comment', 'brite', 'dblinks']

    # Define the query to identify promising compounds for manual addition
    compounds_manual_add_query = """sequence.isna() & type.isna() & ~remark.str.count(r"Same as: 	G").eq(0) & reaction.notna() \
                                & ~name.str.lower().str.contains("protein|globin|doxin|glycan|lase|peptide|rna|dna|steroid|lipid|lignin", na=False) \
                                & ~comment.fillna('').str.lower().str.contains("peptide|protein|[KO:", na=False, regex=False) \
                                & ~formula.str.contains("X", na=False) & ~brite.str.contains("rotein|nzyme|eptide", na=False)"""

    # Preprocess the KEGG compound metadata and apply the query to identify promising compounds
    missing_promising = preprocess_kegg_c_metadata(valid_cids=sorted(k), keep_only=keep_only).query(
        compounds_manual_add_query).rename_axis('compound_id').sort_index().reset_index()

    # Extract keywords from the 'name' field and drop rows with names ending in 'ase'
    kwds = missing_promising['name'].fillna('').str.strip('[|]|(|)').str.split().explode().dropna()
    missing_promising = missing_promising.drop(index=kwds[kwds.str.endswith('ase')].index, errors='ignore').drop(
        ['sequence', 'type', 'remark', 'brite'], axis=1)
    cols = 'compound_id,name,formula,reaction,comment,dblinks'.split(',')
    missing_promising = missing_promising.loc[:, cols]

    # Save the promising compounds to a file
    missing_promising.to_csv('../data/C_IDs_good.dat', encoding='utf-8', index=False)

    return missing_promising
