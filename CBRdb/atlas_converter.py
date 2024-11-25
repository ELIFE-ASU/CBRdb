import os

import pandas as pd

from .tools_eq import standardise_eq
from .tools_files import make_custom_id


def cleanup_eq_line(eq_line):
    """
    Cleans up the equation line by replacing specific characters and patterns.

    Parameters:
    eq_line (str): The equation line to be cleaned.

    Returns:
    str: The cleaned equation line.
    """
    replacements = {"-0": "1",
                    "(0)": "1",
                    ")": " ",
                    "(": "",
                    "<=>": " <=> ",
                    "<==>": " <=> ",
                    "+": " + "}
    for old, new in replacements.items():
        eq_line = eq_line.replace(old, new)
    return eq_line


def cleanup_ec_line(ec_line):
    """
    Cleans up the EC line by selecting the appropriate element based on specific conditions.

    Parameters:
    ec_line (list): A list of strings representing the EC line components.

    Returns:
    str: The cleaned EC line.
    """
    # Select the second to last element and split by the delimiter
    outline = ec_line[-2].split("/")[-1]
    # If the outline is empty, use the 4th element
    if outline == "":
        outline = ec_line[4]
    return outline


def clean_kegg_atlas(in_file="../../data/atlas_kegg_reactions.dat",
                     out_file="../data/atlas_data_kegg_R.csv.zip"):
    """
    Cleans and processes KEGG atlas reaction data.

    Parameters:
    in_file (str): The path to the input file containing KEGG atlas reactions.
    out_file (str): The path to the output file where cleaned data will be saved.

    Returns:
    None
    """
    # Get the absolute paths
    in_file = os.path.abspath(in_file)
    out_file = os.path.abspath(out_file)

    # Check if the file exists
    if not os.path.exists(in_file):
        print("File does not exist", flush=True)
        print("You need to put the atlas data here", flush=True)
        print(in_file, flush=True)
        raise FileNotFoundError

    # Open the file
    with open(in_file, "r") as f:
        # Read the data
        data = f.read()
    # Split the data by new lines
    data = data.split("\n")
    # Initialize the lists
    re_id = []
    re_eq = []
    re_chem_names = []
    re_ec = []

    # Loop over the data
    for i, line in enumerate(data):
        # Split the line by the delimiter
        line = line.split(";")
        id = line[0].strip()

        # Cleanup the equation line
        eq_line = cleanup_eq_line(line[1].strip())
        # If the <=> is not present, replace the last + with <=>
        if "<=>" not in eq_line:
            eq_line = eq_line.rsplit("+", 1)[0] + " <=> " + eq_line.rsplit("+", 1)[1]

        # Standardise the equation
        eq_line = standardise_eq(eq_line)

        # Clean up the chemical names
        chem_names = line[2].replace("|", "").strip()

        # Get EC
        ec = line[3].strip()

        # Get the reaction id
        re_id.append(id)

        # Get the reaction equation
        re_eq.append(eq_line)
        # Get the reaction name
        re_chem_names.append(chem_names)
        # Get the reaction EC
        re_ec.append(ec)
    # Store the data in a dataframe
    df = pd.DataFrame({'id': re_id, 'reaction': re_eq, 'chemical_names': re_chem_names, 'ec': re_ec})
    # Write the data to a file
    df.to_csv(out_file, compression='zip', encoding='utf-8', index=False)
    print("data written to file", flush=True)
    return None


def clean_atlas(in_file="../../data/atlas_reactions.dat",
                out_file="../data/atlas_data_R.csv.zip",
                f_exclude_kegg=True):
    """
    Cleans and processes atlas reaction data, optionally excluding KEGG reactions.

    Parameters:
    in_file (str): The path to the input file containing atlas reactions.
    out_file (str): The path to the output file where cleaned data will be saved.
    f_exclude_kegg (bool): Flag to exclude KEGG reactions from the output. Default is True.

    Returns:
    None
    """
    # Get the absolute paths
    in_file = os.path.abspath(in_file)
    out_file = os.path.abspath(out_file)

    # Check if the file exists
    if not os.path.exists(in_file):
        print("File does not exist", flush=True)
        print("You need to put the atlas data here", flush=True)
        print(in_file, flush=True)
        raise FileNotFoundError

    # Open the file
    with open(in_file, "r") as f:
        # Read the data
        data = f.read()
    # Split the data by new lines
    data = data.split("\n")
    # Initialize the lists
    re_id = []
    re_kegg_id = []
    re_eq = []
    re_ec = []

    # Loop over the data
    for i, line in enumerate(data):
        # Split the line by the delimiter
        line = line.split(";")
        # Get the reaction id
        id = make_custom_id(line[0], prefix="A", digits=6)
        # Get the reaction KEGG id
        kegg_id = line[1]
        # Cleanup the equation line
        eq_line = cleanup_eq_line(line[3])
        # Standardise the equation
        eq_line = standardise_eq(eq_line)
        # Get the reaction EC
        ec = cleanup_ec_line(line)
        # Append the data to the lists
        re_id.append(id)
        re_kegg_id.append(kegg_id)
        re_eq.append(eq_line)
        re_ec.append(ec)

    # Create a DataFrame from the lists
    df = pd.DataFrame({'id': re_id, 'kegg_id': re_kegg_id, 'reaction': re_eq, 'ec': re_ec})
    # Fill in the missing values with NaN
    df = df.replace("", float("NaN"))

    if f_exclude_kegg:
        # Select only the data which has a NaN kegg_id
        df = df[df["kegg_id"].isna()]
        # Remove the kegg_id column
        df = df.drop(columns=["kegg_id"])

    # Write the data to a file
    df.to_csv(out_file, compression='zip', encoding='utf-8', index=False)
    print("data written to file", flush=True)
    return None
