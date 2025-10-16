import pandas as pd
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


def merge_duplicate_reactions(df, r_dupemap):
    """
    Merges duplicate reactions in a DataFrame based on a user-provided duplicate map.

    Parameters:
    df (pd.DataFrame): The DataFrame containing reaction data with KEGG fields.
    r_dupemap (pd.Series): A Series mapping old reaction IDs to new unique reaction IDs.

    Returns:
    pd.DataFrame: A DataFrame with merged duplicate reactions.
    # TODO: add a check for any columns not yet handled here
    """
    all_rns = id_indexed(df.copy(deep=True))
    # Standardize CBRdb_C_ids column to space-separated strings
    if 'CBRdb_C_ids' in all_rns.columns:
        if not isinstance(all_rns['CBRdb_C_ids'].dropna().iloc[0], str):
            all_rns.update(all_rns['CBRdb_C_ids'].map(lambda x: ' '.join(sorted(set(x)))))
    # Focus on duplicate IDs
    df = all_rns.loc[r_dupemap.index].assign(eqn_set=r_dupemap)
    # Since bridgit scores are source-specific, tie to most_sim_kegg
    str_scores = (df['bridgit_score'].map(lambda x: f'({x:.3f}) ', na_action='ignore'))
    df['most_sim_kegg'] = (str_scores + df['most_sim_kegg']).dropna()
    # Group by duplicate reaction IDs
    dupes = df.groupby(by='eqn_set')
    # Label cases where more than one value exists
    to_combine = dupes.nunique().gt(1)
    # First, for each entry, capture the properties that lack multiple values.
    group_attrs = dupes.first(skipna=True).mask(to_combine)

    # In some columns, no combinations occur.
    single_val_cols = to_combine.columns[to_combine.any().eq(False)].to_list()
    # In some columns, the values are structured as space-separated lists.
    known_space_sep_cols = ['ec', 'orthology', 'pathway', 'rclass', 'rhea', 'kegg_id',
                            'msk_ecs', 'msk_metacyc', 'msk_mnxr', 'msk_rhea', 'msk_rns', 'CBRdb_C_ids']
    # For each dupe-group, concatenate its constituent lists.
    space_sep_entries = (df[known_space_sep_cols]
                         .apply(lambda x: x.str.split())
                         .groupby(df['eqn_set'])
                         .sum()[to_combine]
                         .replace(0, float('NaN')))
    # Now sort and dedupe those lists...
    space_sep_entries = space_sep_entries.map(lambda x: ' '.join(sorted(set(x))), na_action='ignore')
    # ... and update the group attributes with the non-nan values
    group_attrs.update(space_sep_entries)

    # Next, isolate the remaining attributes w/multiple values.
    cols = to_combine.columns.difference(single_val_cols + known_space_sep_cols)
    rows = to_combine.index[to_combine[cols].any(axis=1)]
    df_r = dupes.filter(lambda x: x.name in rows)[cols.union(['eqn_set'])]
    # For the comment column, store its provenance to ease interpretation
    df_r.update((df_r.index + ': ' + df_r['comment']).rename('comment'))
    # For the comment + name columns, combine the entries
    tilde_sep_cols = ['comment', 'name']
    for col in tilde_sep_cols:
        tilde_sep_entries = df_r.dropna(subset=col).groupby('eqn_set')[col].apply(
            lambda x: ' ~ '.join(sorted(set(x))))
        tilde_sep_entries = tilde_sep_entries[to_combine[col]]
        group_attrs.update(tilde_sep_entries)

    # For the most_sim_kegg column, if multiple values exist, label w/source and combine.
    new_most_sim_kegg = ((df_r.index + ': ' + df_r['most_sim_kegg'])
                         .rename(df_r['eqn_set']).dropna()
                         .groupby(level=0).apply(' ~ '.join))
    new_most_sim_kegg = new_most_sim_kegg.rename('most_sim_kegg')[to_combine['most_sim_kegg']]
    group_attrs.update(new_most_sim_kegg)
    # bridgit_score is NaN in group_attrs where multiple scores exist; keep it that way for now.
    group_attrs['id_orig'] = {i: ' '.join(sorted(j)) for i,j in dupes.groups.items()}
    # Add id_orig column for all entries.
    all_rns['id_orig'] = all_rns.index
    # Merge group_attrs into main reaction dataframe.
    all_rns.drop(r_dupemap.index, inplace=True)
    all_rns.drop(columns='eqn_set', inplace=True, errors='ignore')
    all_rns = pd.concat([all_rns, group_attrs]).reset_index(names='id')

    return all_rns


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
    cols2unify = C_main_copy.columns.intersection(space_sep_str_cols_cps)
    # consider only those columns
    to_combine = C_dupemap.join(C_main_copy[cols2unify])
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
    str_cols_dict = {k: str for k in space_sep_str_cols_cps}
    final_output_Rs = id_indexed(pd.read_csv(final_output_Rs_fp, index_col=0, dtype=str_cols_dict, low_memory=False))
    final_output_Cs = id_indexed(pd.read_csv(final_output_Cs_fp, index_col=0, dtype=str_cols_dict, low_memory=False))

    rid2cid = final_output_Rs['reaction'].str.findall(r'C\d{5}').map(set).map(sorted).rename('compound_id')
    cid2rid = rid2cid.explode().reset_index().groupby('compound_id')['id'].apply(sorted)

    final_output_Rs['CBRdb_C_ids'] = rid2cid.map(' '.join)
    final_output_Cs['CBRdb_R_ids'] = cid2rid.map(' '.join)

    compound_csv(df_C=final_output_Cs, file_address=final_output_Cs_fp)
    reaction_csv(df_R=final_output_Rs, file_address=final_output_Rs_fp)

    return None
