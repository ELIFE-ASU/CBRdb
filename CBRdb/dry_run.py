import os
import pandas as pd
from .preprocessor import preprocess

def dry_run():
    # Load the

    # Preprocess the reactions data
    preprocess(target="R",
               target_dir=r"../../data/kegg_data",
               out_file=r"../data/kegg_data",
               cid_manual_file=r"../data/C_IDs_manual.dat")

