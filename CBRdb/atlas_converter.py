import os
import warnings

import pandas as pd

from .tools_eq import standardise_eq


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
    in_file = os.path.abspath(in_file).replace('_kegg', '')
    out_file = os.path.abspath(out_file)

    warnings.showwarning(f'\n\tclean_kegg_atlas() will be removed soon. \
                  \n\tFor now, writing CSV of KEGG reactions in ATLAS using clean_atlas. \
                  \n\tFor more reaction info, join with KEGG_R on kegg_id.', FutureWarning, filename='', lineno='')

    df = clean_atlas(in_file=in_file, out_file=out_file, f_exclude_kegg=False).dropna(subset="kegg_id")
    df.to_csv(out_file, compression='zip', encoding='utf-8', index=False)

    return df


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

    # Open the file. For header description, see p.6: https://lcsb-databases.epfl.ch/pathways/atlas/files/ATLAS_UserGuide.pdf
    df = (pd.read_table(in_file, header=None, sep=';', usecols=[0, 1, 3, 4, 5, 6],
                        converters={5: lambda x: -1.0 * pd.to_numeric(x, errors='coerce')},
                        # sign reversed from pdf above
                        names=['id', 'kegg_id', 'reaction', 'reaction_rule', 'dG_kJ/mol', 'dG_err'])
          .assign(id=lambda x: 'A' + x.id.astype(str).str.zfill(6)))  # format ATLAS reactions as: AXXXXXX

    # Extract each reaction's list of 3rd-level EC#s. Remove non-conforming EC#s (e.g. not just numbers and -).
    rr = (df['reaction_rule'].str.replace('-rev)', '-(rev)').str.split('|').explode().str.rstrip('(rev)').to_frame()
          .assign(format_ok=lambda x: x.reaction_rule.str.replace('.', '').str.replace('-', '').str.isnumeric())
          .query('format_ok').groupby(level=0)['reaction_rule'].apply(lambda x: ' '.join(set(x))))
    df = df.join(rr.rename('ec'), how='left').drop('reaction_rule', axis=1)

    # Standardize format of reaction equation.
    df['reaction'] = df['reaction'].apply(cleanup_eq_line).apply(standardise_eq)

    if f_exclude_kegg:
        df = df.query('kegg_id.isna()').drop('kegg_id', axis=1)

    # Write the data to a file
    df.to_csv(out_file, compression='zip', encoding='utf-8', index=False)
    print("data written to file", flush=True)
    return df
