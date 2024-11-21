import pandas as pd

import CBRdb

if __name__ == "__main__":
    print("Program started", flush=True)
    # load the data
    data = pd.read_csv("../data/kegg_data_R.csv.zip", compression='zip')  # , index_col=0
    print("data loaded", flush=True)
    print("data shape", data.shape, flush=True)
    # print the data columns
    print("data columns", data.columns, flush=True)
    # Select the reaction column
    print(data["reaction"].tolist(), flush=True)

    # load to compound data
    data_c = pd.read_csv("../data/kegg_data_C.csv.zip")
    # print("data columns", data_c.columns, flush=True)
    # print(data_c[data_c["compound_id"] == "C00462"].values, flush=True)

    eq_bad_1 = "-0 C03561 <=> C06143 + 0 C00010"
    eq_bad_2 = "1 C00009 + 7 C00620 + 3 C00076 + -5 C02352 <=> 4 C08136 + 1 C01889"

    eq = "-0 C00002 + 1 C00001 <=> 4 C00009 + 3 C00008"
    eq = "-1 C03561 <=> C06143 + 1 C00010"
    # breaks!
    eq = "n-1 C03561 <=> n C06143 + 1 C00010 + C00011"

    # side_to_dict
    lhs = eq.split('<=>')[0]
    lhs = CBRdb.side_to_dict(lhs)
    rhs = eq.split('<=>')[1]
    rhs = CBRdb.side_to_dict(rhs)
    print(f'{lhs} <=> {rhs}', flush=True)

    eq = "1 C00027 + 2 C00126 + C00282 <=> 2 C00001 + 2 C00125"
    reactants, products, react_ele, prod_ele = CBRdb.get_elements_from_eq(eq, data_c)

    eq_sd = CBRdb.standardise_eq(eq)
    print(eq_sd, flush=True)

    reactants, products, react_ele, prod_ele = CBRdb.get_elements_from_eq(eq_sd, data_c)
    print("Reactants: ", reactants, flush=True)
    print("Products:  ", products, flush=True)

    print(CBRdb.check_eq_unbalanced(react_ele, prod_ele), flush=True)

    print("Program end", flush=True)
