import pandas as pd
from .tools_files import reaction_csv, compound_csv

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


def merge_duplicate_reactions(df, r_dupemap):
    """
    Merges duplicate reactions in a DataFrame based on a user-provided duplicate map.

    Parameters:
    df (pd.DataFrame): The DataFrame containing reaction data with KEGG fields.
    r_dupemap (pd.Series): A Series mapping old reaction IDs to new unique reaction IDs.

    Returns:
    pd.DataFrame: A DataFrame with merged duplicate reactions.
    """
    # Replace reaction IDs with the duplicate map
    df['eqn_set'] = df['id'].replace(r_dupemap)

    # Define functions for aggregating data
    unique_str_sep = lambda x: ' '.join(sorted(list(set([i for i in ' '.join(x).split()]))))
    all_provided = lambda x: ' | '.join(sorted(list(x[x != ''].unique())))

    # Dictionary of aggregation functions for each column
    func_dict = {
        'reaction': lambda x: x.iloc[0],
        'id': unique_str_sep,
        'ec': unique_str_sep,
        'pathway': unique_str_sep,
        'orthology': unique_str_sep,
        'rhea': unique_str_sep,
        'module': unique_str_sep,
        'name': all_provided,
        'comment': all_provided,
        'rclass': lambda x: ' '.join(sorted(list(set(x.str.findall(r'(RC\d{5}__C\d{5}_C\d{5})').sum())))),
    }

    # Group by the new reaction IDs and aggregate the data
    deduped_df = (df.fillna('').groupby(by='eqn_set').aggregate(func_dict)
                  ).replace('', float('nan')).reset_index().rename(
        columns={'id': 'id_orig', 'eqn_set': 'id'})

    return deduped_df


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
    sum_entry_strs = lambda x: x.dropna().str.split(' ').sum()
    sort_union_join = lambda x: ' '.join(sorted(set(x)))
    line_set_str = lambda x: '; '.join(sorted(set(x))) # unique strings, sorted
    lines_set_str = lambda x: ';~'.join(dict.fromkeys(x.str.split(';~').sum()).keys())
    name_funcs = {'name': lines_set_str, 'nickname': line_set_str}

    # standardize data format
    C_main_copy = id_indexed(C_main.copy(deep=True))
    # combine names for each duplicate compound group
    combo_names = C_dupemap.join(C_main_copy).groupby('new_id').agg(name_funcs)
    # identify columns for which the value should reflect the union of values
    unify_col_options = ['kegg_reaction', 'kegg_enzyme', 'kegg_pathway', 'kegg_brite', 'kegg_module', 'kegg_glycan',
                         'PDB_CCD', 'ATC_code' , 'Drug_group', 'kegg_type', 'kegg_network',
                         'kegg_drug', 'PubChem', 'ChEBI', 'CAS', 'NIKKAJI', 'KNApSAcK', 'LIPIDMAPS']
    cols2unify = C_main_copy.columns.intersection(unify_col_options)
    # consider only those columns
    to_combine = C_dupemap.join(C_main_copy[cols2unify])
    # format values appropriately - where present, should be strings
    to_combine['PubChem'] = to_combine['PubChem'].map(lambda x: str(x).replace(',0', ''), na_action='ignore')
    C_main_copy['PubChem'] = C_main_copy['PubChem'].map(lambda x: str(x).replace(',0', ''), na_action='ignore')
    # combine each entry's (list of) values
    to_combine = to_combine.reset_index().groupby(by='new_id').agg(sum_entry_strs)
    # de-duplicate and sort each entry's (list of) values; cast as a string
    to_combine = to_combine.replace(0, pd.NA).map(sort_union_join, na_action='ignore')
    # rename the index
    to_combine.rename_axis('compound_id', inplace=True)
    # replace duplicate compound IDs with their new IDs
    C_main_copy.rename(index=C_dupemap['new_id'], inplace=True)
    # sort by compound ID
    C_main_copy.sort_index(inplace=True)
    # update with combined values
    C_main_copy.update(to_combine)
    # update with combined names
    C_main_copy.update(combo_names)
    # return index to initial state
    C_main_copy.reset_index(inplace=True)
    # de-duplicate compound entries
    C_main_copy.drop_duplicates(subset='compound_id', keep='first', inplace=True)

    return C_main_copy


def add_R_col_to_C_file(final_output_Cs_fp='../CBRdb_C.csv', final_output_Rs_fp='../CBRdb_R.csv'):
    """
    Adds a column to the compound DataFrame indicating which reactions each compound is involved in, and vice-versa.
    Parameters:
    final_output_Cs_fp (str): File path for the final output compound DataFrame.
    final_output_Rs_fp (str): File path for the final output reaction DataFrame.  
    Returns:
    None: The function modifies the compound and reaction DataFrames in place and saves them to the specified files. 
    """
    final_output_Rs = id_indexed(pd.read_csv(final_output_Rs_fp, index_col=0, dtype=object))
    final_output_Cs = id_indexed(pd.read_csv(final_output_Cs_fp, index_col=0, dtype=object))

    rid2cid = final_output_Rs['reaction'].str.findall(r'C\d{5}').map(set).map(sorted).rename('compound_id')
    cid2rid = rid2cid.explode().reset_index().groupby('compound_id')['id'].apply(sorted)

    final_output_Rs['CBRdb_C_ids'] = rid2cid.map(' '.join)
    final_output_Cs['CBRdb_R_ids'] = cid2rid.map(' '.join)

    compound_csv(df_C=final_output_Cs, file_address=final_output_Cs_fp)
    reaction_csv(df_R=final_output_Rs, file_address=final_output_Rs_fp)

    return None
