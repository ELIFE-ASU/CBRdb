import pandas as pd

import CBRdb


def merge_atlas_kegg_r(atlas_file="../data/atlas_data_R.csv", kegg_file="../data/kegg_data_R.csv"):
    atlas_data = pd.read_csv(atlas_file)
    atlas_data = atlas_data.sort_values(by="kegg_id")

    print("id", atlas_data["id"].tolist(), flush=True)
    print("data loaded", flush=True)
    print("data shape", atlas_data.shape, flush=True)
    print("data columns", atlas_data.columns, flush=True)

    # Get the max index
    max_index = atlas_data["kegg_id"].tolist()
    print(f"max index: {max_index}", flush=True)

    # Load the kegg data
    kegg_data = pd.read_csv(kegg_file)
    kegg_data = kegg_data.sort_values(by="id")
    print("data loaded", flush=True)
    print("data shape", kegg_data.shape, flush=True)
    print("data columns", kegg_data.columns, flush=True)

    # get the indexes that are not in the atlas data set but are in the kegg data set
    missing_indexes = kegg_data.index.difference(atlas_data.index)
    print(f"missing indexes: {missing_indexes}", flush=True)

    return None


if __name__ == "__main__":
    print("Program started", flush=True)
    merge_atlas_kegg_r()
    exit()

    # r = Reaction.from_string("H2O -> H+ + OH-")
    # # Balance the reaction
    # r_balanced = balance_stoichiometry(r)
    # print(r_balanced, flush=True)

    # ferricyanide = Substance.from_formula('Fe(CN)6-3')
    # print(ferricyanide.composition, flush=True)
    #
    #
    # a = chemparse.parse_formula('Fe(CN)6')
    # b = chemparse.parse_formula('Fe(CN)6-3')
    # c = chemparse.parse_formula('-3')
    #
    # reac, prod = balance_stoichiometry({'NH4ClO4', 'Al'}, {'Al2O3', 'HCl', 'H2O', 'N2'})
    # print(reac, prod, flush=True)

    # load the data
    data = pd.read_csv("../data/kegg_data_R.csv")
    print("data loaded", flush=True)
    print("data shape", data.shape, flush=True)
    # print the data columns
    print("data columns", data.columns, flush=True)
    # Select the reaction column
    print(data["reaction"].tolist(), flush=True)
    # Make a list of the reactions which contain the string "n" or "m" in them
    reactions = data["reaction"].tolist()
    reactions = [r for r in reactions if "n" in r or "m" in r or "x" in r]
    print(reactions, flush=True)
    for r in reactions:
        print(r, flush=True)

    for r in reactions:
        # split the reaction into the reactants and products
        r = r.split(" <=> ")
        for i in r:
            print(f'{i} => {CBRdb.side_to_dict(i)}', flush=True)
    exit()

    # load to compound data
    data_c = pd.read_csv("../data/kegg_data_C.csv")
    print("data columns", data_c.columns, flush=True)
    # print(data_c[data_c["compound_id"] == "C18091"].values, flush=True)
    # print(data_c[data_c["compound_id"] == "C00126"].values, flush=True)
    # print(data_c[data_c["compound_id"] == "C00125"].values, flush=True)

    # Case for n, is to check if it is balanced

    eq = "C11113 + 1 C11131 <=> 1 C11114 + 1 C11181 + 1 C11198"
    eq = "n C00001 + m C00404 <=> n+1 C02174"
    eq = "2n+1 C00001 + 1 C00404 + n+10 C00002 <=> n+1 C02174"
    eq = CBRdb.standardise_eq(eq)

    # # Convert the Eq into the dicts
    # reactants, products = CBRdb.eq_to_dict(eq)
    # print(f"reactants: {reactants}", flush=True)
    # print(f"products: {products}", flush=True)
    #
    # converted_reactants, converted_products = CBRdb.get_formulas_from_eq(eq, data_c)
    # print(converted_reactants, flush=True)
    # print(converted_products, flush=True)
    #
    # # convert all the dict values to strings
    # converted_reactants = {k: str(v) for k, v in converted_reactants.items()}
    # converted_products = {k: str(v) for k, v in converted_products.items()}
    #
    # print(CBRdb.contains_var_list(converted_reactants, converted_products), flush=True)
    # # in the reactants get all the values in the dict
    # reactants_values = list(converted_reactants.values())
    # print(reactants_values, flush=True)
    # print(CBRdb.solve_for(reactants_values), flush=True)

    # Convert the Eq into the formula dicts
    reactants, products = CBRdb.get_formulas_from_eq(eq, data_c)

    print(CBRdb.check_vars_eq_balanced(reactants, products), flush=True)

    print("Program end", flush=True)
