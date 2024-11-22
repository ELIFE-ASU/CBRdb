import pandas as pd
import chemparse
import CBRdb
import re

import re
import chemparse
from chempy import Substance
from chempy import balance_stoichiometry
from sympy import symbols
from chempy import Reaction

import re
from collections import defaultdict

import re


def side_to_dict(s):
    pattern = r'([+-]?[\w*/+-]*)\s*(C\d{5})(?:\(([^)]*)\))?'
    result = {}
    for match in re.finditer(pattern, s):
        coeff_before, cnum, coeff_inside = match.groups()
        coeff_before = coeff_before.strip()
        if coeff_inside:
            coeff = coeff_inside.strip()
        elif coeff_before:
            coeff = coeff_before
        else:
            coeff = 1
        if cnum in result:
            result[cnum] += '+' + coeff
        else:
            result[cnum] = coeff
    # If the coefficient is just a + sign then replace it with 1
    for cnum, coeff in result.items():
        if coeff == '+':
            result[cnum] = 1
    # Convert the coefficients to integers
    for cnum, coeff in result.items():
        try:
            result[cnum] = int(coeff)
        except ValueError:
            pass
    return result


def clean_up_eq(eq):
    # Find 'n+1 C02616' and replace it with '(n+1) C02616'
    eq = re.sub(r'([A-Z]\d{5})\(([^)]+)\)', r'(\2) \1', eq)
    return eq


def side_to_dict(s):
    s = clean_up_eq(s)
    coeff = re.split(r"[A-Z]\d{5}", s)
    # strip the white spaces
    coeff = [c.strip() for c in coeff]

    # if the coeff is equal to "+" then replace it with 1
    coeff = [c if c != "+" else "1" for c in coeff]

    # if the coeff is equal to "+ " then replace it with ""
    coeff = [c.replace("+ ", "").strip() for c in coeff]

    # replace the empty strings with 1
    coeff = [1 if c == '' else c for c in coeff]

    # try to convert the strings to integers
    coeff_out = []
    for c in coeff:
        try:
            coeff_out.append(int(c))
        except ValueError:
            coeff_out.append(c)

    # Get the matches
    matches = re.findall(r"[A-Z]\d{5}", s)
    # Create the dictionary
    return dict(zip(matches, coeff_out))


if __name__ == "__main__":
    print("Program started", flush=True)
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
    data = pd.read_csv("../data/kegg_data_R.csv.zip", compression='zip')  # , index_col=0
    # print("data loaded", flush=True)
    # print("data shape", data.shape, flush=True)
    # # print the data columns
    # print("data columns", data.columns, flush=True)
    # # Select the reaction column
    # print(data["reaction"].tolist(), flush=True)
    # Make a list of the reactions which contain the string "n" or "m" in them
    reactions = data["reaction"].tolist()
    reactions = [r for r in reactions if "n" in r or "m" in r or "x" in r]
    # split the reactions into the reactants and products

    for r in reactions:
        # split the reaction into the reactants and products
        r = r.split(" <=> ")
        for i in r:
            print(f'{i} => {CBRdb.side_to_dict(i)}', flush=True)
    exit()

    # load to compound data
    data_c = pd.read_csv("../data/kegg_data_C.csv.zip")
    # print("data columns", data_c.columns, flush=True)
    # print(data_c[data_c["compound_id"] == "C18091"].values, flush=True)
    # print(data_c[data_c["compound_id"] == "C00126"].values, flush=True)
    # print(data_c[data_c["compound_id"] == "C00125"].values, flush=True)

    eq = "1 C00007 + 2 C00339 + C01438 + n C00001 + n + 1 C00002"
    print(eq, flush=True)
    print(side_to_dict(eq))
    exit()

    eq = CBRdb.standardise_eq(eq)
    print(eq, flush=True)
    # Convert the Eq into the dicts
    reactants, products = CBRdb.eq_to_dict(eq)
    print(reactants, flush=True)
    print(products, flush=True)

    # Get the conversion of the ids to formulas
    react_id_form_key = CBRdb.get_ids_to_formulas(reactants, data_c)
    prod_id_form_key = CBRdb.get_ids_to_formulas(products, data_c)

    print(react_id_form_key, flush=True)
    print(prod_id_form_key, flush=True)

    # Convert the reactants into formulas
    converted_reactants = CBRdb.convert_ids_to_formulas(reactants, react_id_form_key)
    converted_products = CBRdb.convert_ids_to_formulas(products, prod_id_form_key)

    print(converted_reactants, flush=True)
    print(converted_products, flush=True)

    # reactants, products, react_ele, prod_ele = CBRdb.get_elements_from_eq(eq, data_c)
    #
    # eq_sd = CBRdb.standardise_eq(eq)
    # print(eq_sd, flush=True)
    #
    # reactants, products, react_ele, prod_ele = CBRdb.get_elements_from_eq(eq_sd, data_c)
    # print("Reactants: ", reactants, flush=True)
    # print("Products:  ", products, flush=True)
    #
    # print(CBRdb.check_eq_unbalanced(react_ele, prod_ele), flush=True)

    print("Program end", flush=True)
