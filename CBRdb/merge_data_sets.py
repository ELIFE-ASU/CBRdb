import pandas as pd
import os

from .tools_files import reaction_csv, compound_csv, space_sep_str_cols_cps

out_fmt = {'encoding': 'utf-8', 'index': False}


def id_indexed(df: pd.DataFrame) -> pd.DataFrame:
    """ Convenience fxn: Sets 'id' or 'compound_id' column as a pd.DataFrame's index. Returns a copy (not modified inplace) """
    if df.index.name in ['id', 'compound_id']:
        return df
    elif 'id' in df.columns and df.index.name != 'id':
        return df.set_index('id')
    elif 'compound_id' in df.columns and df.index.name != 'compound_id':
        return df.set_index('compound_id')
    else:
        return df


def identify_duplicate_compounds(C_main):
    """
    Identifies duplicate compounds in a DataFrame and maps them to new unique IDs.

    Parameters:
    C_main (pd.DataFrame): The main DataFrame containing compound data with 'compound_id' and 'smiles' columns.

    Returns:
    pd.DataFrame: A DataFrame mapping old compound IDs to new unique compound IDs.
    """
    # ID the number to start counting back from when instantiating new IDs
    id_num = C_main.reset_index()['compound_id'].str.lstrip('C').astype(int)
    count_back_from = id_num.loc[id_num.diff().idxmax()] - 1
    # Flag duplicates among well-defined structures
    possible_dupes = C_main.dropna(subset='smiles').query(
        'smiles.duplicated(keep=False) & ~smiles.str.contains("*", regex=False)')
    # Group duplicate structures together
    C_dupemap = possible_dupes.groupby(by='smiles')['compound_id'].apply(sorted)
    # Assign a new ID for each dupe-group
    C_dupemap = C_dupemap.reset_index(drop=True).explode().rename('old_id')
    C_dupemap.index = ('C' + (count_back_from - C_dupemap.index).astype(str)).set_names('new_id')
    # Sort the duplicate compound groups (within and between)
    C_dupemap = C_dupemap.reset_index().sort_values(by=['new_id', 'old_id']).set_index('old_id')
    # Generate dupe-map output file for easy conversion from old IDs to new IDs
    C_dupemap.to_csv('../data/kegg_data_C_dupemap.csv', encoding='utf-8')
    return C_dupemap


def merge_duplicate_compounds(C_main: pd.DataFrame, C_dupemap: pd.DataFrame) -> pd.DataFrame:
    """
    Merges duplicate compounds in a compound DataFrame using the duplication-map DataFrame.
    Combines the lists of alternate IDs and relationships for all duplicated compounds.

    Parameters:
    C_main (pd.DataFrame): A DataFrame with compound data and compound_id column or index.name
    C_dupemap (pd.DataFrame): A DataFrame mapping duplicate compounds (old_id index) to their merged new_id.

    Returns:
    pd.DataFrame: a copy of the main compound DataFrame, but with the duplicate compounds merged.
    """
    # define functions for combining entries
    line_set_str = lambda x: '; '.join(sorted(set(x)))  # unique strings, sorted
    lines_set_str = lambda x: ';~'.join(dict.fromkeys(x.str.split(';~').sum()).keys())
    name_funcs = {'kegg_name': lines_set_str, 'nickname': line_set_str}
    tripleslash_join = lambda x: '///'.join(sorted(set(x.dropna()))) if not x.dropna().empty else None
    list_unique = lambda x: ' '.join(sorted(set(x.dropna().str.split(' ').explode()))) if not x.dropna().empty else None

    # standardize data format
    C_main_copy = id_indexed(C_main.copy(deep=True))
    # Assign new IDs to dupes
    duped_compounds = C_main_copy.loc[C_dupemap.index].assign(compound_id=C_dupemap).reset_index()
    # Flag presence/absence of repeated subunit indicator "n" in kegg_formula
    duped_compounds['n'] = duped_compounds['kegg_formula'].str.contains('n', na=False, regex=False)
    # Group by new IDs
    dupes = duped_compounds.groupby('compound_id')
    # Identify which fields do/don't have conflicting values
    n_vals = dupes.nunique()
    # Take the first non-NaN entry in fields/IDs without conflicting values
    group_attrs = dupes.first()[n_vals < 2].drop(columns=['old_id', 'n'])
    # For fields where values reflect space-separated lists...
    sscol = group_attrs.columns.intersection(set(space_sep_str_cols_cps.keys()))
    # Take the union of everything listed
    dupes_ss = dupes[sscol].agg(list_unique)[n_vals > 1]
    group_attrs.fillna(dupes_ss, inplace=True)

    # Procedurally combine nickname + name fields
    group_attrs.fillna(dupes.agg(name_funcs), inplace=True)

    # For the kegg_brite_full and comment columns...
    label_and_merge = ['comment', 'kegg_brite_full']
    for col in label_and_merge:
        # Store provenance to ease interpretation
        duped_compounds.update((duped_compounds['old_id'] + ': ' + duped_compounds[col]).rename(col))
        # Then merge those fields
        group_attrs.fillna(duped_compounds.groupby('compound_id')[col].apply(tripleslash_join), inplace=True)

    # The only columns left w/conflicts should be 'kegg_mol_weight', 'kegg_exact_mass', and 'kegg_formula'
    one_val = n_vals.columns[n_vals.lt(2).all()]
    cols_done = sscol.union(label_and_merge).union(one_val).union(name_funcs.keys())
    # Only one dupe group (C98917) has multiple values for mass/weight.
    # For self-consistency WRT kegg_formula (also multiple values), get its first entry whole-cloth.
    first_weight = id_indexed(dupes[['compound_id', 'kegg_mol_weight', 'kegg_exact_mass', 'kegg_formula']].head(1))[
        n_vals > 1].dropna(how='any')
    group_attrs.fillna(first_weight, inplace=True)

    # Prioritize retaining kegg_formula values with "n" in them
    sorted_formulas = duped_compounds.sort_values(by=['n', 'old_id'], ascending=[False, True]).groupby('compound_id')
    sorted_formulas = sorted_formulas[['kegg_formula']].first().drop(first_weight.index)[n_vals > 1].dropna()
    group_attrs.fillna(sorted_formulas, inplace=True)

    # Ensure no columns remain un-merged
    unmerged = C_main_copy.columns.difference(
        cols_done.union(['n', 'old_id', 'kegg_formula', 'kegg_mol_weight', 'kegg_exact_mass']))
    if len(unmerged) > 0:
        print("Warning: The following columns were not merged:", unmerged.tolist())
        print("Assuming first non-NaN entry is acceptable.")
        group_attrs.fillna(dupes[unmerged].first(), inplace=True)

    # Now, remove the duplicate entries from the main DataFrame and add the merged entries
    C_main_copy = pd.concat([C_main_copy.drop(C_dupemap.index), group_attrs])
    # sort by compound ID
    C_main_copy.sort_index(inplace=True)
    # return index to initial state
    C_main_copy.reset_index(inplace=True)

    return C_main_copy


def add_R_col_to_C_file(reaction_df : pd.DataFrame, C_file='../CBRdb_C_metadata.csv.zip'):
    """
    Adds a column to the compound metadata indicating which reactions each compound is involved in.
    Parameters:
    reaction_df (pd.DataFrame): Reaction DataFrame.
    C_file (str): File path for the compound metadata DataFrame.
    Returns:
    None: The function modifies the compound and reaction DataFrames in place and saves them to the specified files. 
    """
    print("Adding column in compound metadata to indicate associated reactions", flush=True)
    c2r = list_reactions_per_compound(reaction_df).map(' '.join)
    join_col_to_csv(col=c2r, file=C_file, col_name='CBRdb_R_ids')
    print("Finished appending reaction lists to compound metadata", flush=True)
    return None


def merge_hpc_calculations(final_output_Cs_fp: str|pd.DataFrame = '../CBRdb_C.csv.zip',
                           formation_energies_fp='../hpc/CBRdb_C_formation_energies.csv.gz',
                           mace_spectrum_fp='../hpc/CBRdb_C_mace_spectrum.csv.gz',
                           assembly_index_fp='../hpc/CBRdb_C_assembly_index.csv.zip',
                           c_dupemap_fp='../data/kegg_data_C_dupemap.csv'):
    """ 
    Merges HPC calculations into the main compounds data, overwriting those columns if present.

    """
    print("Merging HPC calculations into compound file or DataFrame", flush=True)

    # Specify file reading and writing parameters
    f_params_in = dict(index_col=0, low_memory=False, dtype=space_sep_str_cols_cps)
    f_params_out = dict(encoding='utf-8', index=True, compression='infer')

    # Import the datasets
    if not isinstance(final_output_Cs_fp, (str, pd.DataFrame)):
        raise ValueError("final_output_Cs_fp must be one of (str, pd.DataFrame)")
    elif isinstance(final_output_Cs_fp, str):
        data_c = pd.read_csv(os.path.abspath(final_output_Cs_fp), **f_params_in)
    else:
        not_id_indexed = final_output_Cs_fp.index.astype(str).str.isnumeric().all()
        data_c = id_indexed(final_output_Cs_fp)

    data_c_dupes = pd.read_csv(os.path.abspath(c_dupemap_fp), **f_params_in)
    data_c_formation = pd.read_csv(os.path.abspath(formation_energies_fp), **f_params_in)
    data_c_mace = pd.read_csv(os.path.abspath(mace_spectrum_fp), **f_params_in)
    data_c_aix = pd.read_csv(os.path.abspath(assembly_index_fp), **f_params_in)

    # For MACE and AI, prepare to merge on SMILES string.
    smi2cp = {v: k for k,v in data_c['smiles'].items()}
    data_c_mace.set_index('smiles', inplace=True)
    data_c_aix.set_index('smiles', inplace=True)
    mace_eligible = data_c_mace.filter(items=data_c['smiles'], axis=0).rename(smi2cp)
    aix_eligible = data_c_aix[['assembly_index']].filter(items=data_c['smiles'], axis=0).rename(smi2cp)
    data_c_to_add = mace_eligible.join(aix_eligible, how='outer').rename_axis('compound_id')

    # For formation energies, prepare to merge on ID by syncing IDs.
    data_c_formation.rename(index=data_c_dupes.iloc[:,0], inplace=True)

    # Remove previous compound calculations, if present.
    data_c.drop(columns=data_c_to_add.columns, errors='ignore', inplace=True)
    data_c.drop(columns=data_c_formation.columns, errors='ignore', inplace=True)

    # Merge the datasets
    data_c = data_c.join(data_c_to_add, how='left').join(data_c_formation, how='left')

    # Save the updated compounds file (note that floats are not rounded unlike in CBRdb.compound_csv)
    data_c.to_csv('../CBRdb_C.csv.zip', **f_params_out)

    print("HPC calculation merger complete", flush=True)

    if isinstance(final_output_Cs_fp, pd.DataFrame):
        if not_id_indexed:
            return data_c.reset_index()
        else:
            return data_c
    else:
        return None


def merge_hpc_thermo_params(final_output_Rs_fp='../CBRdb_R.csv.zip',
                             thermo_params_fp='../hpc/CBRdb_R_reaction_energies.csv.gz'):
    """ 
    Merges reaction thermodynamic parameters into the main reactions data file, overwriting it.

    """
    print("Merging HPC thermodynamic parameters into reaction file", flush=True)

    # Import the datasets
    reactions = pd.read_csv(final_output_Rs_fp, index_col=0, low_memory=False)
    thermo_params = pd.read_csv(thermo_params_fp, index_col=0, compression='gzip')

    # Merge the datasets
    reactions = reactions.join(thermo_params.add_prefix('thermo_'), how='left')

    # Save the updated reactions file (note that floats are not rounded unlike in CBRdb.reaction_csv)
    reactions.to_csv(final_output_Rs_fp)
    reactions.to_csv(final_output_Rs_fp + '.zip', compression='zip')

    print("HPC thermodynamic parameter merger complete", flush=True)

    return None


def separate_compound_metadata(final_output_Cs_fp="../CBRdb_C.csv",
                               Cs_metadata_fp="../CBRdb_C_metadata.csv",
                               union=False):
    """ Gets metadata from compound file, places it in a separate metadata file, removes from original. """
    params = {'encoding': 'utf-8', 'index': True, 'float_format': '%.3f'}
    C_main = pd.read_csv(final_output_Cs_fp, index_col=0, low_memory=False)
    general = ["comment", "CBRdb_R_ids"]
    metadata = (C_main.filter(like='xref_')
                .join(C_main.filter(like='kegg_'))
                .join(C_main.filter(items=general))
                .copy(deep=True))
    C_main.drop(columns=metadata.columns, inplace=True)

    if union:
        try:
            metadata_orig = pd.read_csv(Cs_metadata_fp, index_col=0, low_memory=False)
            metadata_orig.drop(columns=metadata.columns, inplace=True)
            metadata = metadata.join(metadata_orig, how='left')
        except Exception as e:
            pass

    metadata.to_csv(Cs_metadata_fp, **params)
    metadata.to_csv(Cs_metadata_fp+'.zip', **params)
    C_main.to_csv(final_output_Cs_fp, **params)
    C_main.to_csv(final_output_Cs_fp+'.zip', **params)

    return None




def list_reactions_per_compound(data_r : pd.DataFrame):
    cps = list_compounds_per_reaction(data_r)
    r2c = cps.explode().reset_index(name=0)
    c2r = r2c.groupby(0)[r2c.columns[0]].apply(set).map(sorted)
    return c2r


def list_compounds_per_reaction(data_r : pd.DataFrame):
    cps = id_indexed(data_r)['reaction'].str.findall(r'(C\d+)')
    r2c = cps.map(set).map(sorted)
    return r2c


def join_col_to_csv(col : pd.Series,
                    file : str,
                    col_name: str|None = None,
                    write_params=None):
    
    df = pd.read_csv(file, index_col=0, low_memory=False)

    icol = col.copy(deep=True)
    if isinstance(col_name, str):
        icol.rename(col_name, inplace=True)

    df.drop(columns=icol.name, inplace=True, errors='ignore')

    if df.index.intersection(icol.index).any():
        df = df.join(icol, how='left')
    elif len(df.index) == len(icol.index):
        df.loc[:,icol.name] = icol.tolist()
    else:
        raise ValueError('column and destination file must share indices or index length')
        return None

    if write_params is None:
        write_params = {'encoding': 'utf-8'}

    df.to_csv(file, index=True, **write_params)
    
    return None


