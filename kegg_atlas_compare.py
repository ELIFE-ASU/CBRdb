import os
import time
import numpy as np
import re
import pandas as pd

def extract_all_ids_from_eq(eq):
    return [id for id in eq.split() if id.startswith("C")]

def extract_ids_from_eq(eq):
    # Split the equation
    lhs, rhs = eq.split("<=>")
    # Extract the ids from the lhs and rhs
    lhs_ids = [id for id in lhs.split() if id.startswith("C")]
    rhs_ids = [id for id in rhs.split() if id.startswith("C")]
    return lhs_ids, rhs_ids


if __name__ == "__main__":
    print("Program started", flush=True)
    lhs, rhs = extract_ids_from_eq("C00001 + C00002 <=> C00003 + C00004")
    print(lhs)
    print(rhs)

    # Load the kegg reactions list
    kegg_data = pd.read_csv("Data/kegg_data_R.csv.zip", index_col=0)
    # Load the atlas reactions list
    atlas_data = pd.read_csv("Data/atlas_reactions.csv.zip", index_col=0)
    print(kegg_data.head())
    kegg_reactions = kegg_data[kegg_data['Reaction'].notna()]
    atlas_reactions = atlas_data[atlas_data['Reaction'].notna()]

    # loop over the kegg reactions
    for idx, row in kegg_reactions.iterrows():
        print(f"Processing {idx}", flush=True)
        # get the ids from the equation
        kegg_ids = set(extract_all_ids_from_eq(row['Reaction']))
        # loop over the atlas reactions
        for idx2, row2 in atlas_reactions.iterrows():
            # get the ids from the equation
            atlas_ids = set(extract_all_ids_from_eq(row2['Reaction']))
            # check if the kegg ids are in the atlas ids
            if all([id in atlas_ids for id in kegg_ids]):
                print(f"Matched {idx} {idx2}", flush=True)
                break


    print("Program ended", flush=True)
