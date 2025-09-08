import pandas as pd

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

    sum_entry_strs = lambda x: x.dropna().str.split(' ').sum()
    sort_union_join = lambda x: ' '.join(sorted(set(x)))

    # standardize data format
    C_main_copy = id_indexed(C_main.copy(deep=True))
    combo_names = C_dupemap.join(C_main_copy[['name']]).groupby('new_id')['name'].agg(';~'.join)
    combo_names = combo_names.str.split(';~').map(lambda x: ';~'.join(dict.fromkeys(x).keys())).to_frame()
    # identify columns for which the value should reflect the union of values
    unify_col_options = ['kegg_reaction', 'kegg_enzyme', 'kegg_pathway', 'kegg_brite', 'kegg_module', 'kegg_glycan', 
                  'kegg_drug', 'PubChem', 'ChEBI', 'CAS', 'NIKKAJI', 'KNApSAcK', 'LIPIDMAPS']
    cols2unify = C_main_copy.columns.intersection(unify_col_options)
    # consider only those columns
    to_combine =  C_dupemap.join(C_main_copy[cols2unify])
    # format values appropriately - where present, should be strings
    to_combine.update(to_combine['PubChem'].dropna().astype(int).astype(str))
    # combine each entry's (list of) values
    to_combine = to_combine.reset_index().groupby(by='new_id').agg(sum_entry_strs).replace(0, pd.NA)
    # de-duplicate and sort each entry's (list of) values; cast as a string
    to_combine = to_combine.map(sort_union_join, na_action='ignore')
    # rename the index
    to_combine.rename_axis('compound_id', inplace=True)
    # replace duplicate compound IDs with their new IDs
    C_main_copy.rename(index=C_dupemap['new_id'], inplace=True)
    # sort by compound ID
    C_main_copy.sort_index(inplace=True)
    # update with combined values
    C_main_copy.update(to_combine)
    C_main_copy.update(combo_names)
    # return index to initial state
    C_main_copy.reset_index(inplace=True)
    # de-duplicate compound entries
    C_main_copy.drop_duplicates(subset='compound_id', keep='first', inplace=True)

    return C_main_copy


def add_R_col_to_C_file(final_output_Cs_fp = '../CBRdb_C.csv', final_output_Rs_fp = '../CBRdb_R.csv'):
    """
    Adds a column to the compound DataFrame indicating which reactions each compound is involved in.
    Parameters:
    final_output_Cs_fp (str): File path for the final output compound DataFrame.
    final_output_Rs_fp (str): File path for the final output reaction DataFrame.  
    Returns:
    None: The function modifies the compound DataFrame in place and saves it to the specified file  
    """
    final_output_Rs = pd.read_csv(final_output_Rs_fp, index_col=0, usecols=[0,1,2], dtype=object)
    final_output_Cs = pd.read_csv(final_output_Cs_fp, index_col=0)

    if 'id' in final_output_Rs.columns:
        final_output_Rs = final_output_Rs.set_index('id')
    if 'compound_id' in final_output_Cs.columns:
        final_output_Cs = final_output_Cs.set_index('compound_id')
    if 'CBRdb_R_ids' in final_output_Cs.columns:
        final_output_Cs.drop(columns=['CBRdb_R_ids'], inplace=True)    

    cid2rid = pd.Series(final_output_Rs.reaction.str.findall(r'C\d{5}').explode().to_frame().groupby('reaction').groups)
    cid2rid_printable = cid2rid.map(lambda x: ' '.join(sorted(list(set(x)), reverse=True)))

    final_output_Cs['CBRdb_R_ids'] = cid2rid_printable
    final_output_Cs.to_csv(final_output_Cs_fp, index=True)
    return None
