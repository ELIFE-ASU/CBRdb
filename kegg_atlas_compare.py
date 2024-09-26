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


def main():
    # Load the kegg reactions list
    kegg_data = pd.read_csv("Data/kegg_data_R.csv.zip", index_col=0)
    print(kegg_data.head())

    # Load the atlas reactions list
    atlas_data = pd.read_csv("Data/atlas_data_R.csv.zip", index_col=0)
    print(atlas_data.head())

    # Filter out the reactions that don't have a reaction
    kegg_reactions = kegg_data[kegg_data['reaction'].notna()]
    atlas_reactions = atlas_data[atlas_data['reaction'].notna()]

    matched_ids = []

    # Loop over the kegg reactions
    for idx, row in kegg_reactions.iterrows():
        print(f"Processing {idx}", flush=True)
        # Get the ids from the equation
        kegg_ids = set(extract_all_ids_from_eq(row['reaction']))
        # Loop over the atlas reactions
        for idx2, row2 in atlas_reactions.iterrows():
            # Get the ids from the equation
            atlas_ids = set(extract_all_ids_from_eq(row2['reaction']))
            # Check if the kegg ids are in the atlas ids
            if all([id in atlas_ids for id in kegg_ids]):
                print(f"Matched {idx} {idx2} {kegg_data["id"][idx]}", flush=True)
                # remove the matched ids from the atlas ids
                matched_ids.append(idx2)

    print("matched_ids", matched_ids, flush=True)


if __name__ == "__main__":
    print("Program started", flush=True)
    main()
    print("Program ended", flush=True)
