import pandas as pd
import numpy as np


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
            print(f"Processing {i}/{N} atlas IDs", flush=True)
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
            print(f"Problem {atlas_id} no forward or reverse IDs!", flush=True)
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
    return None


def cleanup_eq_line(eq_line):
    eq_line = eq_line.replace(")", " ")
    eq_line = eq_line.replace("(", "")
    eq_line = eq_line.replace("<==>", " <==> ")
    eq_line = eq_line.replace("+", " + ")
    return eq_line


def clean_kegg_atlas(in_file, out_file):
    # open the file
    with open(in_file, "r") as f:
        # read the data
        data = f.read()
    # split the data by new lines
    data = data.split("\n")
    # init the lists
    re_id = []
    re_eq = []
    re_chem_names = []
    re_ec = []

    # loop over the data
    for i, line in enumerate(data):
        # split the line by the delimiter
        line = line.split(";")
        # Get the reaction id
        re_id.append(line[0])
        # Get the reaction equation
        re_eq.append(cleanup_eq_line(line[1]))
        # Get the reaction name
        re_chem_names.append(line[2].replace("|", ""))
        # Get the reaction EC
        re_ec.append(line[3].replace("|", ""))
    # Store the data in a dataframe
    df = pd.DataFrame({'id': re_id, 'reaction': re_eq, 'chemical_names': re_chem_names, 'ec': re_ec})
    # Write the data to a file
    df.to_csv(out_file, compression='zip', encoding='utf-8', index=False)
    return None


def clean_atlas(in_file, out_file):
    # Open the file
    with open(in_file, "r") as f:
        # Read the data
        data = f.read()
    # split the data by new lines
    data = data.split("\n")
    # init the lists
    re_id = []
    re_kegg_id = []
    re_eq = []
    re_ec = []

    # loop over the data
    for i, line in enumerate(data):
        # split the line by the delimiter
        line = line.split(";")
        # Get the reaction id
        re_id.append(line[0])
        # Get the KEGG reaction id
        re_kegg_id.append(line[1])
        # Get the reaction equation
        re_eq.append(cleanup_eq_line(line[3]))
        # Get the reaction EC
        re_ec.append(line[4].replace("|", ""))
    # Store the data in a dataframe
    df = pd.DataFrame({'id': re_id, 'kegg_id': re_kegg_id, 'reaction': re_eq, 'ec': re_ec})
    # Write the data to a file
    df.to_csv(out_file, compression='zip', encoding='utf-8', index=False)
    return None


if __name__ == "__main__":
    print("Program start", flush=True)
    # infile = r"C:\Users\louie\skunkworks\data\atlas_reactions.csv"
    # outfile = r"C:\Users\louie\skunkworks\data\atlas_reactions.csv.zip"
    # convert_atlas(infile, outfile)

    infile = r"C:\Users\louie\Downloads\atlas_kegg_reactions.dat"
    outfile = r"Data\atlas_kegg_reactions.csv.zip"
    clean_kegg_atlas(infile, outfile)

    infile = r"C:\Users\louie\Downloads\atlas_reactions.dat"
    outfile = r"Data\atlas_reactions.csv.zip"
    clean_atlas(infile, outfile)

    print("Program end", flush=True)
