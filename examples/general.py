import pandas as pd
import chemparse
import CBRdb
import re

import re
import chemparse

def strip_ionic_states(formula):
    # Regex to match trailing ionic states (+, -, +2, -4, etc.)
    return re.sub(r'[+-]\d*$', '', formula)

def convert_formula_to_dict(formula):
    # Count the occurrences of '*' followed by an optional number
    star_count = sum(int(num) if num else 1 for num in re.findall(r'\*(\d*)', formula))

    # Replace '*' followed by an optional number with an empty string
    formula_new = re.sub(r'\*\d*', '', formula)

    # Remove the ionic states
    formula_new = strip_ionic_states(formula_new)

    # Parse the formula into a dictionary
    formula_dict = {k: int(v) for k, v in chemparse.parse_formula(formula_new).items()}

    # Add the '*' count to the dictionary if there were any '*'
    if star_count > 0:
        formula_dict['*'] = star_count

    return formula_dict

if __name__ == "__main__":
    print("Program started", flush=True)
    # load the data
    data = pd.read_csv("../data/kegg_data_R.csv.zip", compression='zip')  # , index_col=0
    print("data loaded", flush=True)
    print("data shape", data.shape, flush=True)
    # print the data columns
    print("data columns", data.columns, flush=True)
    # # Select the reaction column
    # print(data["reaction"].tolist(), flush=True)

    # load to compound data
    data_c = pd.read_csv("../data/kegg_data_C.csv.zip")
    print("data columns", data_c.columns, flush=True)
    print(data_c[data_c["compound_id"] == "C18091"].values, flush=True)
    print(data_c[data_c["compound_id"] == "C00126"].values, flush=True)
    print(data_c[data_c["compound_id"] == "C00125"].values, flush=True)

    exit()

    #print(data_c["formula"].tolist(), flush=True)

    tmp = convert_formula_to_dict("C2H2*BrO2")
    print(tmp, flush=True)
    tmp = convert_formula_to_dict("Te+1")
    print(tmp, flush=True)
    tmp = CBRdb.convert_formula_to_dict("C2H2*BrO2")
    print(tmp, flush=True)
    tmp = CBRdb.convert_formula_to_dict("Te-")
    print(tmp, flush=True)
    tmp = CBRdb.convert_formula_to_dict("C2H4NO2-")
    print(tmp, flush=True)



    exit()

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
