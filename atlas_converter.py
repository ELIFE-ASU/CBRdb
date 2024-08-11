import pandas as pd


def convert_atlas(in_file, out_file):
    # Load the data using pandas
    data = pd.read_csv(in_file)
    # Get the set of the ids
    atlas_ids = sorted(list(set(data["rn"])))
    N = len(atlas_ids)
    atlas_list = []
    eq_list = []
    # Loop over the atlas ids
    for i, atlas_id in enumerate(atlas_ids):
        # print after 100 iterations
        if i % 100 == 0:
            print(f"Processing {i}/{N} atlas IDs")
        all_reactions = data[data['rn'] == atlas_id]
        # Get the forward parts
        f_react = all_reactions[all_reactions['direction'] == "forward"]
        # Get the reverse parts
        r_react = all_reactions[all_reactions['direction'] == "reverse"]
        # Try and use the reactant bits
        if len(f_react) != 0:
            item = f_react.values.tolist()
        # If there are no forward parts, try the reverse parts
        elif len(r_react) != 0:
            item = r_react.values.tolist()
        # If there are no reverse parts, print an error and exit
        else:
            print(f"Problem {atlas_id} no forward or reverse IDs!")
            exit()
        # Get the reactants and products
        lhs = []
        rhs = []
        # Loop over the entries
        for j, entry in enumerate(item):
            # If the entry is negative, it is a reactant
            if entry[4] < 0:
                lhs.append(f"{int(abs(entry[4]))} {entry[3]}")
            # If the entry is positive, it is a product
            else:
                rhs.append(f"{int(abs(entry[4]))} {entry[3]}")
        lhs = " + ".join(lhs)
        rhs = " + ".join(rhs)
        eq_list.append(f"{lhs} <=> {rhs}")
        atlas_list.append(atlas_id)

    # Make the dataframe from the eq_list and the atlas_ids
    df = pd.DataFrame({'ID': atlas_list, 'Reaction': eq_list})
    # Write the data to a file
    df.to_csv(out_file, compression='zip', encoding='utf-8')


if __name__ == "__main__":
    print("Program start", flush=True)
    infile = r"C:\Users\louie\skunkworks\data\atlas_reactions.csv"
    outfile = r"C:\Users\louie\skunkworks\data\atlas_reactions.csv.zip"
    convert_atlas(infile, outfile)
    print("Program end", flush=True)
