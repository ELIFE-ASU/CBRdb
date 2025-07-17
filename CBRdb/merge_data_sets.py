import pandas as pd

out_fmt = {'encoding': 'utf-8', 'index': False}


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
    id_num = C_main.reset_index()['compound_id'].str.lstrip('C').astype(int)
    count_back_from = id_num.loc[id_num.diff().idxmax()] - 1
    possible_dupes = (
        C_main.query(
            '~smiles.str.contains("*", regex=False) & smiles.duplicated(keep=False)').reset_index().sort_values(
            by=['smiles', 'compound_id'])
        .groupby('smiles')['compound_id'].apply(list).apply(sorted).reset_index(drop=True)
        .explode().reset_index(name='id_old').rename({'index': 'id_new'}, axis=1))
    possible_dupes['id_new'] = 'C' + (count_back_from - possible_dupes['id_new']).astype(str)
    possible_dupes = dict(zip(possible_dupes['id_old'], possible_dupes['id_new']))
    compound_mapping = pd.Series(possible_dupes).reset_index().groupby(by=0)['index'].apply(
        lambda x: ' '.join(sorted(list(x))))
    compound_mapping = (compound_mapping.str.split().explode().rename('old_id')
                        .reset_index().rename({0: 'new_id'}, axis=1).set_index('old_id'))
    compound_mapping.to_csv('../data/kegg_data_C_dupemap.csv', encoding='utf-8')
    return compound_mapping


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
