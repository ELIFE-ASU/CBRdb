import numpy as np
import pandas as pd

def get_formulas_from_ids(ids, file_path='Data/kegg_data_C.csv.zip'):
    data = pd.read_csv(file_path)
    return data.loc[data["compound_id"].isin(ids), "formula"].tolist()

if __name__ == "__main__":
    # ['Unnamed: 0', 'compound_id', 'smiles', 'formula', 'molecular_weight',
    #        'n_heavy_atoms', 'n_chiral_centers']
    print("Program started", flush=True)
    file_path = 'Data/kegg_data_C.csv.zip'
    data = pd.read_csv(file_path)
    print("Data loaded", flush=True)
    print("Data columns", data.columns, flush=True)
    print("Data shape", data.shape, flush=True)
    print("Data head", data.head(4).values, flush=True)