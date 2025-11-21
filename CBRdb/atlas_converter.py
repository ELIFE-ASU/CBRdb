import os

import pandas as pd

from .tools_eq import standardise_eq
from .tools_files import reaction_csv


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


def clean_atlas(in_file="../../data/atlas_reactions.dat",
                out_file="../data/atlas_data_R.csv",
                f_exclude_kegg=False,
                extract_msk=False):
    """
    Cleans and processes atlas reaction data, optionally excluding KEGG reactions.

    Parameters:
    in_file (str): The path to the input file containing atlas reactions.
    out_file (str): The path to the output file where cleaned data will be saved.
    f_exclude_kegg (bool): Flag to exclude KEGG reactions from the output. Default is False.
    extract_msk (bool): Flag to extract external IDs from the most_sim_kegg field. Default is False.

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

    # The ATLAS database headers are listed on p.6 of this guide:
    # https://lcsb-databases.epfl.ch/pathways/atlas/files/ATLAS_UserGuide.pdf
    in_file_cols = ['id', 'kegg_id', 'x', 'reaction', 'reaction_rule',
                    'dG_kJ/mol', 'dG_err', 'most_sim_kegg', 'bridgit_score']

    # Open the file.
    df = pd.read_table(in_file, header=None, sep=';', names=in_file_cols)

    # Remove known KEGG entries if requested
    if f_exclude_kegg:
        df = df.query('kegg_id.isna()')
    # kegg_id column is misleadingly incomplete WRT current KEGG data; also we merge dupes later
    df.drop('kegg_id', axis=1, inplace=True)
    # Standardize the ID column format
    df['id'] = 'A' + df['id'].astype(str).str.zfill(6)
    # ATLAS only annotates ECs to the 3rd digit
    ecs = df['reaction_rule'].str.findall(r'(\d\.\d+\.\d+\.[-])')
    df['ec'] = ecs.map(lambda x: ' '.join(sorted(set(x))).strip(), na_action='ignore')
    # Standardize format of reaction equation.
    df['reaction'] = df['reaction'].apply(cleanup_eq_line).apply(standardise_eq)
    # Some 'most_sim_kegg' entries are a single placeholder character, replace w/NaN
    df['most_sim_kegg'] = df['most_sim_kegg'].mask(lambda x: x.str.len().lt(2))
    # Remove unused columns.
    df.drop(columns=['x', 'reaction_rule', 'dG_kJ/mol', 'dG_err'], inplace=True)
    # Extract "most_similar_kegg" tags if requested
    if extract_msk:
        df = extract_msk_tags(df)
    # Write the data to a file
    reaction_csv(df, out_file)
    print("data written to file", flush=True)
    return df


def extract_msk_tags(df):
    """
    Extracts various identifiers from the 'most_sim_kegg' column of the DataFrame.
    Parameters:
    df (pd.DataFrame): The input DataFrame containing the 'most_sim_kegg' column.
    Returns:
    pd.DataFrame: The updated DataFrame with new columns for extracted identifiers.
    """
    # Split the 'most_sim_kegg' column into individual tags and explode into rows
    tags = df['most_sim_kegg'].str.split(r'[|]|[/]').explode()

    # Extract RHEA IDs
    rheas = tags[tags.str.contains('RHEA', na=False)].str.split('RHEA:').explode()
    rheas = rheas.replace('', pd.NA).dropna().str.strip(' ')
    df['msk_rhea'] = rheas.groupby(level=0).apply(list)

    # extract other relevant IDs from the 'most_sim_kegg' column (abbreviated 'msk')
    df['msk_ecs'] = df['most_sim_kegg'].str.findall(r'(\d\.\d+\.\d+\.\d+)')
    df['msk_rns'] = df['most_sim_kegg'].str.replace('BR', '').str.findall(r'(R\d{5})')
    df['msk_metacyc'] = df['most_sim_kegg'].str.findall(r'(RXN-\d+)')
    df['msk_mnxr'] = df['most_sim_kegg'].str.findall(r'(MNXR\d+)')

    for fxn in [set, sorted, ' '.join]:
        df.update(df.filter(like='msk').map(fxn, na_action='ignore'))
    df = df.replace('', float('nan'))
    # dG signs in the guide above are reversed in the data we were granted.
    # TODO: Calculate dG for a few rxns with extreme (+/-) reported dG. Diff signs, y/n?
    # df['dG_kJ/mol'] = -1 * df['dG_kJ/mol'].astype(float)
    # TODO: If we are to keep ATLAS's dG-related fields, should find out why some dG_err are negative.
    return df
