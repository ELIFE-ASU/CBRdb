import os

import pandas as pd

out_fmt = {'compression': 'zip', 'encoding': 'utf-8', 'index': False}


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
    dupemap_file = f'{data_folder}/kegg_data_C_dupemap.csv.zip'
    C_meta_file = f'{data_folder}/kegg_data_C_metadata.csv.zip'
    C_main_file = f'{data_folder}/kegg_data_C.csv.zip'
    atlas_data_R_file = f'{data_folder}/atlas_data_R.csv.zip'
    kegg_data_R_file = f'{data_folder}/kegg_data_R.csv.zip'

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

    # Save the deduped compound files and reaction files
    C_meta.to_csv(C_meta_file, **out_fmt)
    C_main.to_csv(C_main_file, **out_fmt)
    atlas_data_R.to_csv(atlas_data_R_file, **out_fmt)
    kegg_data_R.to_csv(kegg_data_R_file, **out_fmt)

    datasets = dict(zip('C_meta C_main atlas_data_R kegg_data_R dupemap'.split(),
                        [C_meta, C_main, atlas_data_R, kegg_data_R, dupemap]))
    return datasets


def merge_and_create_unique_db(df1, df2, column_name, f_keep='first'):
    """
    Merges two DataFrames and drops duplicates based on a specified column.

    Parameters:
    df1 (pd.DataFrame): The first DataFrame to merge.
    df2 (pd.DataFrame): The second DataFrame to merge.
    column_name (str): The column name to check for duplicates.
    f_keep (str): Which duplicates (if any) to keep. 'first' keeps the first occurrence. Defaults to 'first'.

    Returns:
    pd.DataFrame: The merged DataFrame with duplicates removed.
    """
    # Merge the two DataFrames
    merged_df = pd.concat([df1, df2], ignore_index=True)
    # Drop duplicates, keeping the first occurrence (from df1)
    unique_df = merged_df.drop_duplicates(subset=[column_name], keep=f_keep)
    return unique_df


def get_nonunique_entries(df1, df2, column_name):
    """
    Identifies non-unique entries in the first DataFrame based on a specified column after merging with the second DataFrame.

    Parameters:
    df1 (pd.DataFrame): The first DataFrame to check for non-unique entries.
    df2 (pd.DataFrame): The second DataFrame to merge with the first DataFrame.
    column_name (str): The column name to check for duplicates.

    Returns:
    pd.DataFrame: A DataFrame containing the non-unique entries from the first DataFrame.
    """
    # Merge the two DataFrames
    merged_df = pd.concat([df1, df2], ignore_index=True)

    # Identify rows that are duplicated based on the specified column
    # non_unique_df = merged_df[merged_df.duplicated(subset=[column_name], keep=False)]
    non_unique_df = df1[merged_df.duplicated(subset=[column_name], keep=False)]
    return non_unique_df


def merge_data(merge_col='reaction',
               f_keep='first',
               kegg_file="../data/kegg_data_R_processed.csv.zip",
               atlas_file="../data/atlas_data_R_processed.csv.zip",
               out_file="../data/kegg_atlas_processed_merged.csv.zip"):
    """
    Merges KEGG and ATLAS data sets, drops duplicates based on a specified column, and saves the merged data.

    Parameters:
    merge_col (str): The column name to check for duplicates. Defaults to 'reaction'.
    f_keep (str): Which duplicates (if any) to keep. 'first' keeps the first occurrence. Defaults to 'first'.
    kegg_file (str): The file path to the KEGG data CSV file. Defaults to "../data/kegg_data_R_processed.csv.zip".
    atlas_file (str): The file path to the ATLAS data CSV file. Defaults to "../data/atlas_data_R_processed.csv.zip".
    out_file (str): The file path to save the merged data CSV file. Defaults to "../data/kegg_atlas_processed_merged.csv.zip".

    Returns:
    None
    """
    # Convert the file paths to absolute paths
    kegg_file = os.path.abspath(kegg_file)
    atlas_file = os.path.abspath(atlas_file)
    out_file = os.path.abspath(out_file)

    print("Merging the KEGG and ATLAS data sets", flush=True)
    print(f"KEGG file: {kegg_file}", flush=True)
    print(f"ATLAS file: {atlas_file}", flush=True)
    # Load the kegg reactions list
    kegg_data = pd.read_csv(kegg_file, index_col=0)
    n_kegg = kegg_data.shape[0]

    # Load the atlas reactions list
    atlas_data = pd.read_csv(atlas_file, index_col=0)
    # Drop the kegg id column
    if 'kegg_id' in atlas_data.columns:
        atlas_data = atlas_data.drop(columns=['kegg_id'])
    n_atlas = atlas_data.shape[0]
    # Merge the two databases
    merged_database = merge_and_create_unique_db(kegg_data, atlas_data, merge_col, f_keep=f_keep)
    n_merged = merged_database.shape[0]
    # print the difference in the number of reactions
    print(f"Number of KEGG reactions: {n_kegg}", flush=True)
    print(f"Number of ATLAS reactions: {n_atlas}", flush=True)
    print(f"Number of merged reactions: {n_merged}", flush=True)
    print(f"Number of non-unique reactions: {n_kegg + n_atlas - n_merged}", flush=True)
    # Save the merged database
    merged_database.to_csv(out_file, compression='zip', encoding='utf-8')
    print("Merged data saved! \n", flush=True)

    tmp = get_nonunique_entries(kegg_data, atlas_data, merge_col)
    if tmp.shape[0] > 0:
        print("Non-unique entries that were removed:", flush=True)
        print(tmp, flush=True)


def merge_data_retain_sources(kegg_file="../data/kegg_data_R_processed.csv.zip",
                              atlas_file="../data/atlas_data_R_processed.csv.zip",
                              out_file="../data/atlas_kegg_processed_merged_deduped.csv.zip",
                              ):
    """
    Merges KEGG and ATLAS data sets, standardizes compound order and EC numbers, and saves the merged data.

    Parameters:
    kegg_file (str): The file path to the KEGG data CSV file. Defaults to "../data/kegg_data_R_processed.csv.zip".
    atlas_file (str): The file path to the ATLAS data CSV file. Defaults to "../data/atlas_data_R_processed.csv.zip".
    out_file (str): The file path to save the merged data CSV file. Defaults to "../data/atlas_kegg_processed_merged_deduped.csv.zip".

    Returns:
    None
    """
    # Convert the file paths to absolute paths
    kegg_file = os.path.abspath(kegg_file)
    atlas_file = os.path.abspath(atlas_file)
    out_file = os.path.abspath(out_file)

    print("Merging the KEGG and ATLAS data sets", flush=True)
    print(f"KEGG file: {kegg_file}", flush=True)
    print(f"ATLAS file: {atlas_file}", flush=True)
    # Load the kegg reactions list
    kegg_data = pd.read_csv(kegg_file, index_col=0)
    # Load the atlas reactions list
    atlas_data = pd.read_csv(atlas_file, index_col=0)

    print('Adding coefficients where missing, then standardizing compound order', flush=True)
    add_ones_and_sort_compounds = lambda x: ' + '.join(
        sorted(['1 ' + i if i[0] == 'C' else i for i in x.split(' + ')]))
    atlas_data['reaction_sorted'] = atlas_data['reaction'].str.split(' <=> ').apply(
        lambda x: ' <=> '.join(sorted([add_ones_and_sort_compounds(i) for i in x])))
    kegg_data['reaction_sorted'] = kegg_data['reaction'].str.split(' <=> ').apply(
        lambda x: ' <=> '.join(sorted([add_ones_and_sort_compounds(i) for i in x])))

    print('Standardizing EC# rules between datasets', flush=True)
    single_sep = lambda x: x.strip().replace(' ', '|').replace(',', '|')  # standardizes separators
    extract_ec_serial_numbers = lambda x: '|'.join(sorted(list(set(
        [i for i in single_sep(x).replace('(rev)', '').split('|')
         if '-' not in i and i.count('.') == 3 and i.split('.')[-1].isnumeric()]))))  # some end in e.g. B##
    atlas_data['ec'] = atlas_data['ec'].apply(extract_ec_serial_numbers)

    print('Combining the datasets and preserving EC# sources', flush=True)
    concat_db = pd.concat([kegg_data, atlas_data], axis=0).reset_index()
    concat_db['ec'] = concat_db['ec'].str.split('|')
    concat_db = concat_db.drop(columns='reaction').set_index('reaction_sorted').explode('ec')
    concat_db['id2ec'] = (concat_db['index'] + ':' + concat_db['ec']) * (concat_db['ec'] != '')

    print('Merging the data', flush=True)
    merged_db = (concat_db.applymap(lambda x: [x]).groupby(level=0).sum()
                 .applymap(lambda x: '|'.join(sorted(list(set(x)))[::-1])).reset_index()
                 .rename(columns={'reaction_sorted': 'reaction'}).set_index('index'))
    merged_db.to_csv(out_file, compression='zip', encoding='utf-8')
    print("Merged data saved! \n", flush=True)
