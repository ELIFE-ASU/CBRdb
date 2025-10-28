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
    group_attrs['id_orig'] = {i: ' '.join(sorted(j)) for i, j in dupes.groups.items()}
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
    line_set_str = lambda x: '; '.join(sorted(set(x)))  # unique strings, sorted
    lines_set_str = lambda x: ';~'.join(dict.fromkeys(x.str.split(';~').sum()).keys())
    name_funcs = {'name': lines_set_str, 'nickname': line_set_str}
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
    sscol = group_attrs.columns.intersection(space_sep_str_cols_cps)
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
