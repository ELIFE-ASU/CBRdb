import os
import time

import pandas as pd
import requests
from requests import Session
from requests.adapters import HTTPAdapter
from urllib3.util import Retry


def prepare_session():
    # Make the session
    s = Session()
    # Add retries
    retries = Retry(
        total=5,
        backoff_factor=0.1,
        status_forcelist=[502, 503, 504],
        allowed_methods={'POST'},
    )
    # Mount the session
    s.mount('https://', HTTPAdapter(max_retries=retries))
    return s


def get_ec_ids(session, kegg_website=r"https://rest.kegg.jp/link/enzyme/reaction", request_sleep=0.2):
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
        # Split the response into lines
        lines = response.text.replace("rn:", "").replace("ec:", "").split("\n")[:-1]
        # Split each string by the tab character and create a list of lists
        split_data = [item.split('\t') for item in lines]
        # Create a DataFrame from the list of lists
        df = pd.DataFrame(split_data, columns=['Reaction', 'Enzyme'])
        return df

    else:
        # Some error in the response
        print(f"Error response {response.status_code}")
        exit()


def fix_ec_ids(file_ec_ids="Data/ec_ids.csv.zip",
               input_file="Data/full_processed_merged.csv.zip",
               output_file="Data/full_processed_merged.csv.zip"):
    if not os.path.exists(file_ec_ids):
        session = prepare_session()
        data = get_ec_ids(session)
        # write the data to a file
        data.to_csv(file_ec_ids, compression='zip', encoding='utf-8')

    # Read the data from the file
    data = pd.read_csv(file_ec_ids, compression='zip')
    print(f"Data: {data}", flush=True)

    # read the merged data

    data_merged = pd.read_csv(input_file, compression='zip')

    # # Merge the two DataFrames
    # merged_df = pd.concat([data, data_ec_ids], ignore_index=True)
    # # Drop duplicates, keeping the first occurrence (from df1)
    # unique_df = merged_df.drop_duplicates(subset=['Reaction'], keep='first')
    # return unique_df


def merge_and_create_unique_db(df1, df2, column_name, f_keep='first'):
    # Merge the two DataFrames
    merged_df = pd.concat([df1, df2], ignore_index=True)
    # Drop duplicates, keeping the first occurrence (from df1)
    unique_df = merged_df.drop_duplicates(subset=[column_name], keep=f_keep)
    return unique_df


def get_nonunique_entries(df1, df2, column_name):
    # Merge the two DataFrames
    merged_df = pd.concat([df1, df2], ignore_index=True)

    # Identify rows that are duplicated based on the specified column
    # non_unique_df = merged_df[merged_df.duplicated(subset=[column_name], keep=False)]
    non_unique_df = df1[merged_df.duplicated(subset=[column_name], keep=False)]
    return non_unique_df


def merge_data(merge_col='id',
               f_keep='last',
               kegg_file="Data/kegg_data_R.csv.zip",
               atlas_file="Data/atlas_data_kegg_R.csv.zip",
               out_file="Data/kegg_data_R_merged.csv.zip"):
    print("Merging the KEGG and ATLAS data sets", flush=True)
    print(f"KEGG file: {kegg_file}", flush=True)
    print(f"ATLAS file: {atlas_file}", flush=True)
    # Load the kegg reactions list
    kegg_data = pd.read_csv(kegg_file)
    n_kegg = kegg_data.shape[0]

    # Load the atlas reactions list
    atlas_data = pd.read_csv(atlas_file)
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


if __name__ == "__main__":
    print("Program started", flush=True)
    # merge_data(merge_col='id',
    #            f_keep='last',
    #            kegg_file="Data/kegg_data_R.csv.zip",
    #            atlas_file="Data/atlas_data_kegg_R.csv.zip",
    #            out_file="Data/kegg_data_R_merged.csv.zip")

    merge_data(merge_col='reaction',
               f_keep='first',
               kegg_file="Data/atlas_data_kegg_R_processed.csv.zip",
               atlas_file="Data/atlas_data_R_processed.csv.zip",
               out_file="Data/atlas_kegg_processed_merged.csv.zip")

    merge_data(merge_col='reaction',
               f_keep='first',
               kegg_file="Data/kegg_data_R_processed.csv.zip",
               atlas_file="Data/atlas_data_R_processed.csv.zip",
               out_file="Data/kegg_atlas_processed_merged.csv.zip")
    # fix_ec_ids()
    print("Program ended", flush=True)
