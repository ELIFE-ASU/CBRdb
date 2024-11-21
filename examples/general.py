import pandas as pd
import numpy as np
from rdkit import Chem as Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import rdMolDescriptors
import CBRdb

if __name__ == "__main__":
    print("Program started", flush=True)
    # # load the data
    # data = pd.read_csv("../data/atlas_kegg_processed_merged.csv.zip", compression='zip') # , index_col=0
    # print("data loaded", flush=True)
    # print("data shape", data.shape, flush=True)
    # # print the data columns
    # print("data columns", data.columns, flush=True)
    #
    # # get the entry with the id R04254
    # print(data[data["id"] == "A109140"].values, flush=True)
    #
    #
    # load to compound data
    data_c = pd.read_csv("../data/kegg_data_C.csv.zip")
    # print("data columns", data_c.columns, flush=True)
    # print(data_c[data_c["compound_id"] == "C00462"].values, flush=True)

    eq_bad_1 = "-0 C03561 <=> C06143 + 0 C00010"
    eq_bad_2 = "1 C00009 + 7 C00620 + 3 C00076 + -5 C02352 <=> 4 C08136 + 1 C01889"

    eq = "-0 C00002 + 1 C00001 <=> 4 C00009 + 3 C00008"
    eq = "-1 C03561 <=> C06143 + 1 C00010"
    eq_sd = CBRdb.standardise_eq(eq)
    print(eq_sd, flush=True)

    reactants, products, react_ele, prod_ele = CBRdb.get_elements_from_eq(eq_sd, data_c)
    print("Reactants: ", reactants, flush=True)
    print("Products:  ", products, flush=True)

    print(CBRdb.check_eq_unbalanced(react_ele, prod_ele), flush=True)


    print("Program end", flush=True)