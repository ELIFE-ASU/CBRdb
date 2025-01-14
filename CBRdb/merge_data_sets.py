import os
import time
from io import StringIO

import pandas as pd
import requests

from .tools_requests import prepare_session


def dedupe_compound_files(data_folder='../data'):
    """
    Deduplicates compound files by replacing duplicate compound IDs with unique ones.

    Parameters:
    data_folder (str): The folder path where the compound data files are located. Defaults to '../data'.

    Returns:
    pd.Series: A Series mapping duplicate compound IDs to unique ones.
    """
    dupemap_file = f'{data_folder}/kegg_data_C_dupemap.csv.zip'
    C_meta_file = f'{data_folder}/kegg_data_C_metadata.csv.zip'
    C_main_file = f'{data_folder}/kegg_data_C.csv.zip'

    # Read the duplicate map file
    dupemap = pd.read_csv(dupemap_file, header=0, index_col=0).iloc[:, 0]

    # Read and process the metadata file
    C_meta = (pd.read_csv(C_meta_file, header=0).assign(
        compound_id=lambda x: x['compound_id'].replace(dupemap))
              .drop_duplicates(subset='compound_id', keep='first')
              .sort_values(by='compound_id').reset_index(drop=True))

    # Read and process the main compound file
    C_main = (pd.read_csv(C_main_file, header=0).assign(
        compound_id=lambda x: x['compound_id'].replace(dupemap))
              .drop_duplicates(subset='compound_id', keep='first')
              .sort_values(by='compound_id').reset_index(drop=True))

    # Save the processed metadata and main compound files
    C_meta.to_csv(C_meta_file, index=False, compression='zip')
    C_main.to_csv(C_main_file, index=False, compression='zip')

    return dupemap


def get_ec_ids(session, kegg_website=r"https://rest.kegg.jp/link/enzyme/reaction", request_sleep=0.2):
    """
    Retrieves enzyme-reaction pairs from the KEGG website and returns them as a DataFrame.

    Parameters:
    session (requests.Session): The session object to use for making the request.
    kegg_website (str): The URL of the KEGG website to fetch the enzyme-reaction pairs from. Defaults to "https://rest.kegg.jp/link/enzyme/reaction".
    request_sleep (float): The time to sleep between requests to avoid overloading the server. Defaults to 0.2 seconds.

    Returns:
    pd.DataFrame: A DataFrame containing enzyme-reaction pairs with reactions as the index and enzymes as the values.
    """
    # Get the data
    try:
        # Get the response
        response = session.get(f"{kegg_website}", timeout=10.0)
        # Limit the number of requests
        time.sleep(request_sleep)
    except requests.exceptions.RequestException as e:
        # Some error in the connection
        print(f"Error connection exception {e}", flush=True)
        exit()
    # Check if the response is ok
    if response.ok:
        # Get all enzyme-reaction pairs from KEGG
        df = (pd.read_table(StringIO(response.text), header=None, index_col=None, names=['index', 'ec'])
              .applymap(lambda x: x.split(':')[1]))
        # Make enzyme list for each reaction
        df = df.groupby(by='index')['ec'].apply(lambda x: '|'.join(list(x))).to_frame()
        return df
    else:
        # Some error in the response
        print(f"Error response {response.status_code}")
        exit()


def replace_entries(df1, df2):
    """
    Replaces entries in df1 with corresponding entries from df2.

    Parameters:
    df1 (pd.DataFrame): The first DataFrame to be updated.
    df2 (pd.DataFrame): The second DataFrame with replacement values.

    Returns:
    pd.DataFrame: The updated DataFrame.
    """
    # Ensure the indices and columns match
    df1.update(df2)
    return df1


def fix_ec_ids(file_ec_ids="../data/ec_ids.csv.zip",
               input_file="../data/kegg_atlas_processed_merged.csv.zip",
               output_file=None):
    """
    Fixes EC IDs in the reaction data by replacing them with updated EC IDs from a specified file.

    Parameters:
    file_ec_ids (str): The file path to the EC IDs CSV file. Defaults to "../data/ec_ids.csv.zip".
    input_file (str): The file path to the input reaction data CSV file. Defaults to "../data/kegg_atlas_processed_merged.csv.zip".
    output_file (str, optional): The file path to save the updated reaction data CSV file. If None, the input file path is used. Defaults to None.

    Returns:
    None
    """
    # Check if the output file is specified
    if output_file is None:
        output_file = input_file
    # Convert the file paths to absolute paths
    file_ec_ids = os.path.abspath(file_ec_ids)
    input_file = os.path.abspath(input_file)
    output_file = os.path.abspath(output_file)

    # Grab the EC data from the KEGG website
    if not os.path.exists(file_ec_ids):
        print("Getting the EC data from the KEGG website", flush=True)
        session = prepare_session()
        data = get_ec_ids(session)
        # Write the data to a file
        data.to_csv(file_ec_ids, compression='zip', encoding='utf-8')
        print("EC data saved! \n", flush=True)
    # Read the data from the file
    data_ec = pd.read_csv(file_ec_ids, compression='zip')
    # Read the reaction data
    data_reactions = pd.read_csv(input_file, compression='zip')
    # Replace the entries in the reaction data
    print("Replacing the EC entries in the reaction data", flush=True)
    updated_df = replace_entries(data_reactions, data_ec)
    updated_df.to_csv(output_file, compression='zip', encoding='utf-8')
    print("EC data replaced! \n", flush=True)


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
    print(f"Number of KEGG reactions: {n_kegg}")
    print(f"Number of ATLAS reactions: {n_atlas}")
    print(f"Number of merged reactions: {n_merged}")
    print(f"Number of non-unique reactions: {n_kegg + n_atlas - n_merged}")
    # Save the merged database
    merged_database.to_csv(out_file, compression='zip', encoding='utf-8')
    print("Merged data saved! \n", flush=True)

    tmp = get_nonunique_entries(kegg_data, atlas_data, merge_col)
    if tmp.shape[0] > 0:
        print("Non-unique entries that were removed:")
        print(tmp)


def merge_data_retain_sources(kegg_file="../data/kegg_data_R_processed.csv.zip",
                              atlas_file="../data/atlas_data_R_processed.csv.zip",
                              out_file="../data/atlas_kegg_processed_merged_deduped.csv.zip",
                              ec_file="../data/ec_ids.csv.zip"):
    """
    Merges KEGG and ATLAS data sets, standardizes compound order and EC numbers, and saves the merged data.

    Parameters:
    kegg_file (str): The file path to the KEGG data CSV file. Defaults to "../data/kegg_data_R_processed.csv.zip".
    atlas_file (str): The file path to the ATLAS data CSV file. Defaults to "../data/atlas_data_R_processed.csv.zip".
    out_file (str): The file path to save the merged data CSV file. Defaults to "../data/atlas_kegg_processed_merged_deduped.csv.zip".
    ec_file (str): The file path to the EC IDs CSV file. Defaults to "../data/ec_ids.csv.zip".

    Returns:
    None
    """
    # Convert the file paths to absolute paths
    kegg_file = os.path.abspath(kegg_file)
    atlas_file = os.path.abspath(atlas_file)
    out_file = os.path.abspath(out_file)
    ec_file = os.path.abspath(ec_file)

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
    kegg_ecs = pd.read_csv(ec_file, compression='zip', index_col=0)['ec']
    kegg_data['ec'] = kegg_data.index.map(kegg_ecs).fillna('')  # until kegg_data_R_processed.csv.zip is fixed

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
