import os

import pandas as pd


def cleanup_eq_line(eq_line):
    eq_line = eq_line.replace(")", " ")
    eq_line = eq_line.replace("(", "")
    eq_line = eq_line.replace("<=>", " <=> ")
    eq_line = eq_line.replace("<==>", " <=> ")
    eq_line = eq_line.replace("+", " + ")
    return eq_line


def cleanup_ec_line(ec_line):
    # Select the second to last element and split by the delimiter
    outline = ec_line[-2].split("/")[-1]
    # If the outline is empty, use the 4th element
    if outline == "":
        outline = ec_line[4]
    return outline


def clean_kegg_atlas(in_file, out_file):
    in_file = os.path.abspath(in_file)
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
        re_ec.append(line[3])
    # Store the data in a dataframe
    df = pd.DataFrame({'id': re_id, 'reaction': re_eq, 'chemical_names': re_chem_names, 'ec': re_ec})
    # Write the data to a file
    df.to_csv(out_file, compression='zip', encoding='utf-8', index=False)
    print("Data written to file", flush=True)
    return None


def clean_atlas(in_file, out_file):
    in_file = os.path.abspath(in_file)
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
        # re_ec.append(line[4].replace("|", ""))
        re_ec.append(cleanup_ec_line(line))
    # Store the data in a dataframe
    df = pd.DataFrame({'id': re_id, 'kegg_id': re_kegg_id, 'reaction': re_eq, 'ec': re_ec})
    # Write the data to a file
    df.to_csv(out_file, compression='zip', encoding='utf-8', index=False)
    print("Data written to file", flush=True)
    return None


def main():
    infile = r"..\data\atlas_kegg_reactions.dat"
    outfile = r"Data\atlas_data_kegg_R.csv.zip"
    clean_kegg_atlas(infile, outfile)

    infile = r"..\data\atlas_reactions.dat"
    outfile = r"Data\atlas_data_R.csv.zip"
    clean_atlas(infile, outfile)
    return None


if __name__ == "__main__":
    print("Program start", flush=True)
    main()
    print("Program end", flush=True)
