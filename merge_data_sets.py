import pandas as pd


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
        print("Non-unique entries:")
        print(tmp)


if __name__ == "__main__":
    print("Program started", flush=True)
    merge_data(merge_col='id',
               f_keep='last',
               kegg_file="Data/kegg_data_R.csv.zip",
               atlas_file="Data/atlas_data_kegg_R.csv.zip",
               out_file="Data/kegg_data_R_merged.csv.zip")

    merge_data(merge_col='reaction',
               f_keep='first',
               kegg_file="Data/kegg_data_R_processed.csv.zip",
               atlas_file="Data/atlas_data_R_processed.csv.zip",
               out_file="Data/full_processed_merged.csv.zip")
    print("Program ended", flush=True)
