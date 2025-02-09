import os

import pandas as pd

out_fmt = {'encoding': 'utf-8', 'index': False}
from .tools_files import add_suffix_to_file


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
    unique_str_sep = lambda x: ' '.join(set([i for i in ' '.join(x).split()]))
    all_provided = lambda x: ' | '.join(x[x != ''].unique())

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
        'rclass': lambda x: ' | '.join(set(x.str.findall(r'(RC\d+  C\d+_C\d+)').sum())),
    }

    # Group by the new reaction IDs and aggregate the data
    deduped_df = (df.fillna('').groupby(by='eqn_set').aggregate(func_dict)
                  ).replace('', float('nan')).reset_index().rename(
        columns={'id': 'id_orig', 'eqn_set': 'id'})

    return deduped_df


def dedupe_compounds(data_folder='../data'):
    """
    Deduplicates compound files by replacing duplicate compound IDs with unique ones.

    Parameters:
    data_folder (str): The folder path where the compound, reaction, and dupe-map files are located. Defaults to '../data'.

    Returns:
    dict: dictionary of all reaction and compound datasets with de-duplicated compound IDs, plus the compound dupe-map itself
    """
    dupemap_file = f'{data_folder}/kegg_data_C_dupemap.csv'
    C_meta_file = f'{data_folder}/kegg_data_C_metadata.csv'
    C_main_file = f'{data_folder}/kegg_data_C.csv'
    atlas_data_R_file = f'{data_folder}/atlas_data_R.csv'
    kegg_data_R_file = f'{data_folder}/kegg_data_R.csv'

    # Read the duplicate map file
    dupemap = pd.read_csv(dupemap_file, header=0, index_col=0).iloc[:, 0]

    # Read and process the metadata file
    C_meta = (pd.read_csv(C_meta_file, header=0).assign(
        compound_id=lambda x: x['compound_id'].replace(dupemap))
              .drop_duplicates(subset='compound_id', keep='first')
              .sort_values(by='compound_id'))

    # Read and process the main compound file
    C_main = (pd.read_csv(C_main_file, header=0).assign(
        compound_id=lambda x: x['compound_id'].replace(dupemap))
              .drop_duplicates(subset='compound_id', keep='first')
              .sort_values(by='compound_id'))

    # Read and process the ATLAS reaction file
    atlas_data_R = pd.read_csv(atlas_data_R_file, header=0)
    atlas_data_R['reaction'] = (atlas_data_R['reaction'].str.split(expand=True).replace(dupemap)
                                .fillna('').apply(lambda x: ' '.join(x), axis=1).str.strip())

    # Read and process the KEGG reaction file
    kegg_data_R = pd.read_csv(kegg_data_R_file, header=0)
    kegg_data_R['reaction'] = (kegg_data_R['reaction'].str.split(expand=True).replace(dupemap)
                               .fillna('').apply(lambda x: ' '.join(x), axis=1).str.strip())

    # Save the deduped compound files and reaction files, and tag them as de-duped
    dd_suf = lambda x: add_suffix_to_file(x, '_deduped')
    C_meta.to_csv(dd_suf(C_meta_file), **out_fmt)
    C_main.to_csv(dd_suf(C_main_file), **out_fmt)

    dd_suf = lambda x: add_suffix_to_file(x, '_dedupedCs')
    atlas_data_R.to_csv(dd_suf(atlas_data_R_file), **out_fmt)
    kegg_data_R.to_csv(dd_suf(kegg_data_R_file), **out_fmt)
    
    CBRdb_C = C_main.merge(C_meta, on='compound_id', how='outer')
    CBRdb_C.to_csv(data_folder+'/CBRdb_C.csv', encoding='utf-8', index=False)

    datasets = dict(zip('CBRdb_C C_meta C_main atlas_data_R kegg_data_R dupemap'.split(),
                        [CBRdb_C, C_meta, C_main, atlas_data_R, kegg_data_R, dupemap]))
    return datasets

