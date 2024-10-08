import pandas as pd


def merge_and_keep_duplicates_from_df1(df1, df2, column_name):
    df2_unique = df2[~df2[column_name].isin(df1[column_name])]
    return pd.concat([df1, df2_unique], ignore_index=True)


def main(kegg_file="Data/kegg_data_R_processed.csv.zip", atlas_file="Data/atlas_data_R_processed..csv.zip"):
    print("Merging the KEGG and ATLAS data sets", flush=True)
    print(f"KEGG file: {kegg_file}", flush=True)
    print(f"ATLAS file: {atlas_file}", flush=True)
    # Load the kegg reactions list
    kegg_data = pd.read_csv(kegg_file, index_col=0)
    n_kegg = kegg_data.shape[0]

    # Load the atlas reactions list
    atlas_data = pd.read_csv(atlas_file, index_col=0)
    # Drop the kegg id column
    atlas_data = atlas_data.drop(columns=['kegg_id'])
    n_atlas = atlas_data.shape[0]
    # Merge the two databases
    merged_database = merge_and_keep_duplicates_from_df1(kegg_data, atlas_data, 'reaction')
    n_merged = merged_database.shape[0]
    # print the difference in the number of reactions
    print(f"Number of KEGG reactions: {n_kegg}")
    print(f"Number of ATLAS reactions: {n_atlas}")
    print(f"Number of merged reactions: {n_merged}")
    print(f"Number of non-unique reactions: {n_kegg + n_atlas - n_merged}")
    # Save the merged database
    merged_database.to_csv("Data/full_processed_merged.csv.zip", compression='zip', encoding='utf-8')


if __name__ == "__main__":
    print("Program started", flush=True)
    main()
    print("Program ended", flush=True)
