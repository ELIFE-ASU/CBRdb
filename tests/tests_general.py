import os

import numpy as np
import pandas as pd
from ase.io import read
from ase.visualize import view
from rdkit import Chem as Chem

import CBRdb


def assert_dicts_equal(d1, d2):
    """
    Asserts that two dictionaries are equal by comparing their keys and values.

    Parameters:
    d1 (dict): The first dictionary to compare.
    d2 (dict): The second dictionary to compare.

    Raises:
    AssertionError: If the dictionaries do not have the same keys or if any corresponding values differ.
    """
    # Check that both have the same set of keys
    assert d1.keys() == d2.keys(), f"Key sets differ: {d1.keys()} != {d2.keys()}"

    # Check each corresponding value
    for key in d1.keys():
        assert d1[key] == d2[key], f"Value mismatch for key '{key}': {d1[key]} != {d2[key]}"


def get_unique_elements():
    data_c = pd.read_csv(os.path.abspath("../data/kegg_data_C.csv"))
    smiles_list = data_c['smiles'].tolist()
    # Get the set of unique element symbols from the data
    element_symbols = set()
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            for atom in mol.GetAtoms():
                element_symbols.add(atom.GetSymbol())
    return element_symbols


def test_side_to_dict():
    print(flush=True)
    tmp = CBRdb.side_to_dict('C00001 + 1 C00002')
    assert_dicts_equal(tmp, {'C00001': 1, 'C00002': 1})

    tmp = CBRdb.side_to_dict('C00001 + C00002')
    assert_dicts_equal(tmp, {'C00001': 1, 'C00002': 1})

    tmp = CBRdb.side_to_dict('1 C00001 + 1 C00002')
    assert_dicts_equal(tmp, {'C00001': 1, 'C00002': 1})

    tmp = CBRdb.side_to_dict('+1 C00001 + 1 C00002')
    assert_dicts_equal(tmp, {'C00001': 1, 'C00002': 1})

    tmp = CBRdb.side_to_dict('-1 C00001 + 1 C00002')
    assert_dicts_equal(tmp, {'C00001': -1, 'C00002': 1})

    tmp = CBRdb.side_to_dict('n C00001 + 1 C00002')
    assert_dicts_equal(tmp, {'C00001': 'n', 'C00002': 1})

    tmp = CBRdb.side_to_dict('n+1 C00001 + 1 C00002')
    assert_dicts_equal(tmp, {'C00001': 'n+1', 'C00002': 1})

    tmp = CBRdb.side_to_dict('n-1 C00001 + 1 C00002')
    assert_dicts_equal(tmp, {'C00001': 'n-1', 'C00002': 1})

    tmp = CBRdb.side_to_dict('0 C00001 + 1 C00002')
    assert_dicts_equal(tmp, {'C00002': 1})

    tmp = CBRdb.side_to_dict("1 C00007 + 2 C00339 + C01438")
    assert_dicts_equal(tmp, {'C00007': 1, 'C00339': 2, 'C01438': 1})

    tmp = CBRdb.side_to_dict("C00024 + n C00083 + n C00005 + n C00004 + 2n C00080")
    assert_dicts_equal(tmp, {'C00024': 1, 'C00083': 'n', 'C00005': 'n', 'C00004': 'n', 'C00080': '2n'})

    tmp = CBRdb.side_to_dict("C00003 + C00039(n) + C02128(m)")
    assert_dicts_equal(tmp, {'C00003': 1, 'C00039': 'n', 'C02128': 'm'})

    tmp = CBRdb.side_to_dict("C00020 + C00455 + C00039(n+m)")
    assert_dicts_equal(tmp, {'C00020': 1, 'C00455': 1, 'C00039': 'n+m'})

    tmp = CBRdb.side_to_dict("C03323(m) + C03323(n)")
    assert_dicts_equal(tmp, {'C03323': 'm+n'})

    tmp = CBRdb.side_to_dict("C03323(m-1) + C03323(n+1)")
    assert_dicts_equal(tmp, {'C03323': 'm-1+n+1'})

    tmp = CBRdb.side_to_dict("n C00001 + 1 C00404")
    assert_dicts_equal(tmp, {'C00001': 'n', 'C00404': 1})


def test_convert_formula_to_dict():
    print(flush=True)
    tmp = CBRdb.convert_formula_to_dict("C2H2*BrO2")
    assert_dicts_equal(tmp, {'C': 2, 'H': 2, 'Br': 1, 'O': 2, '*': 1})

    tmp = CBRdb.convert_formula_to_dict("C2H2*32BrO2")
    assert_dicts_equal(tmp, {'C': 2, 'H': 2, 'Br': 1, 'O': 2, '*': 32})

    tmp = CBRdb.convert_formula_to_dict("Te+", strip_ionic=False)
    assert_dicts_equal(tmp, {'Te': 1, '+': 1})

    tmp = CBRdb.convert_formula_to_dict("Te-", strip_ionic=False)
    assert_dicts_equal(tmp, {'Te': 1, '-': 1})

    tmp = CBRdb.convert_formula_to_dict("Te+1", strip_ionic=False)
    assert_dicts_equal(tmp, {'Te': 1, '+': 1})

    tmp = CBRdb.convert_formula_to_dict("Te-1", strip_ionic=False)
    assert_dicts_equal(tmp, {'Te': 1, '-': 1})

    tmp = CBRdb.convert_formula_to_dict("C2H4*NO2-", strip_ionic=False)
    assert_dicts_equal(tmp, {'C': 2, 'H': 4, 'N': 1, 'O': 2, '-': 1, '*': 1})


def test_get_formulas_from_eq():
    print(flush=True)
    data_c = pd.read_csv(os.path.abspath("../data/kegg_data_C.csv"))
    eq = "1 C00001 <=> 1 C00007"
    reactants, products = CBRdb.get_formulas_from_eq(eq, data_c)
    print(reactants, flush=True)
    print(products, flush=True)
    assert_dicts_equal(reactants, {'H2O': 1})
    assert_dicts_equal(products, {'O2': 1})
    print(flush=True)

    # possible case of false positive
    eq = "2 C19610 + C00027 + 2 C00080 <=> 2 C19611 + 2 C00001"  # R00011
    # eq = CBRdb.standardise_eq(eq)
    print(eq, flush=True)

    reactants, products = CBRdb.get_formulas_from_eq(eq, data_c)
    print(reactants, flush=True)
    print(products, flush=True)
    print(flush=True)
    assert_dicts_equal(reactants, {'Mn': 2, 'H2O2': 1, 'H': 2})
    assert_dicts_equal(products, {'Mn': 2, 'H2O': 2})


def test_eq_n_solver():
    print(flush=True)
    expr = "n-1"
    result = CBRdb.find_min_integers(expr)
    print(result)
    assert_dicts_equal(result, {'n': 2})

    expr = "m-1+n+1"
    result = CBRdb.find_min_integers(expr)
    print(result)
    assert_dicts_equal(result, {'n': 1, 'm': 1})

    expr = "2*n"
    result = CBRdb.find_min_integers(expr)
    print(result)
    assert_dicts_equal(result, {'n': 1})

    expr = "2*n + 1"
    result = CBRdb.find_min_integers(expr)
    print(result)
    assert_dicts_equal(result, {'n': 1})

    # This is messed up but works
    expr = "2*(n - 1) + 2*(m - 1)+ (x-10)"
    result = CBRdb.find_min_integers(expr)
    print(result)
    assert_dicts_equal(result, {'x': 11, 'n': 6, 'm': 6})


def test_eq_to_dict():
    print(flush=True)
    eq = "2 C00027 <=> 2 C00001 + 1 C00007"
    reactants, products = CBRdb.eq_to_dict(eq)
    print(reactants, flush=True)
    print(products, flush=True)
    # convert back to string
    eq_out = CBRdb.dicts_to_eq(reactants, products)
    print(eq_out, flush=True)
    assert eq_out == eq


def test_eq_to_symbols():
    print(flush=True)
    data_c = pd.read_csv(os.path.abspath("../data/kegg_data_C.csv"))
    eq = "2 C00027 <=> 2 C00001 + 1 C00007"

    reactants, products, react_ele, prod_ele = CBRdb.get_elements_from_eq(eq, data_c)

    assert_dicts_equal(reactants, {'H2O2': 2})
    assert_dicts_equal(products, {'H2O': 2, 'O2': 1})

    assert_dicts_equal(react_ele, {'O': 4, 'H': 4})
    assert_dicts_equal(prod_ele, {'O': 4, 'H': 4})

    # Convert the dict back into eq form
    eq_out = CBRdb.get_eq(eq, reactants, products, data_c)
    assert eq_out == eq


def test_eq_standard():
    print(flush=True)
    # Check it works when the equation is already standard
    eq = "2 C00027 <=> 2 C00001 + 1 C00007"
    eq_out = CBRdb.standardise_eq(eq)
    assert eq_out == eq
    # Check it works when the equation is not standard, and they are out of order
    eq_wrong = "2 C00027 <=> 1 C00007 + 2 C00001"
    eq_out = CBRdb.standardise_eq(eq_wrong)
    assert eq_out == eq

    # This is a more complex example R00011
    eq = "2 C19610 + C00027 + 2 C00080 <=> 2 C19611 + 2 C00001"
    eq_out = CBRdb.standardise_eq(eq)
    print(eq_out, flush=True)
    assert eq_out == "1 C00027 + 2 C00080 + 2 C19610 <=> 2 C00001 + 2 C19611"


def test_eq_balanced():
    print(flush=True)
    data_c = pd.read_csv(os.path.abspath("../data/kegg_data_C.csv"))

    # # The equation is balanced
    # print("The equation is balanced", flush=True)
    # eq = "2 C00027 <=> 2 C00001 + 1 C00007"
    # reactants, products, react_ele, prod_ele = CBRdb.get_elements_from_eq(eq, data_c)
    # result1 = CBRdb.check_eq_unbalanced(react_ele, prod_ele)
    # result2 = CBRdb.full_check_eq_unbalanced(eq, data_c)
    # print(result1, flush=True)
    # print(result2, flush=True)
    # # assert result1 == False
    # # assert result2 == False
    # print(flush=True)
    #
    # # The equation is not balanced
    # print("The equation is not balanced", flush=True)
    # eq = "2 C00027 <=> 2 C00001 + 1 C00007 + 1 C00008"
    # reactants, products, react_ele, prod_ele = CBRdb.get_elements_from_eq(eq, data_c, strip_ionic=False)
    # result1 = CBRdb.check_eq_unbalanced(react_ele, prod_ele)
    # result2 = CBRdb.full_check_eq_unbalanced(eq, data_c)
    # print(result1, flush=True)
    # print(result2, flush=True)
    # # assert result1 == True
    # # assert result2 == True
    # print(flush=True)

    # possible case of false positive
    eq = "2 C19610 + C00027 + 2 C00080 <=> 2 C19611 + 2 C00001"  # R00011
    # eq = CBRdb.standardise_eq(eq)
    print(eq, flush=True)
    # converted_reactants, converted_products = get_formulas_from_eq(eq, c_data, strip_ionic=strip_ionic)
    #
    reactants, products, react_ele, prod_ele = CBRdb.get_elements_from_eq(eq, data_c)

    print(reactants, flush=True)
    print(products, flush=True)
    print(react_ele, flush=True)
    print(prod_ele, flush=True)

    # check if the equation is balanced
    res = CBRdb.full_check_eq_unbalanced(eq, data_c)
    print(res, flush=True)


def test_eq_difference():
    pass


def test_contains_var_list_check():
    print(flush=True)
    # Test which sees if there is variable in the equation
    data_c = pd.read_csv(os.path.abspath("../data/kegg_data_C.csv"))
    eq = "m C00404 + n C00001 <=> (n+1) C02174 + x C00001"
    reactants, products = CBRdb.get_formulas_from_eq(eq, data_c)
    vars = CBRdb.contains_var_list(reactants, products)
    assert vars == ['m', 'n', 'x']
    assert CBRdb.check_contains_var_list(eq, data_c) == True

    # Test which sees if there is no variable in the equation
    eq = "C00404 + C00001 <=> C02174"
    reactants, products = CBRdb.get_formulas_from_eq(eq, data_c)
    vars = CBRdb.contains_var_list(reactants, products)
    assert vars == []
    assert CBRdb.check_contains_var_list(eq, data_c) == False

    eq = "C00404 + n C00001 <=> (n+1) C02174"
    reactants, products = CBRdb.get_formulas_from_eq(eq, data_c)
    vars = CBRdb.contains_var_list(reactants, products)
    assert vars == ['n']
    assert CBRdb.check_contains_var_list(eq, data_c) == True


def test_vars_eq_balanced():
    ############################### DOUBLE CHECK THIS FUNCTION ########################################
    # print(flush=True)
    # data_c = pd.read_csv(os.path.abspath("../data/kegg_data_C.csv"))
    # eq = "2 C00027 <=> 2 C00001 + 1 C00007"
    # reactants, products, react_ele, prod_ele = CBRdb.get_elements_from_eq(eq, data_c)
    # result = CBRdb.check_eq_unbalanced(react_ele, prod_ele)
    # assert result == False  # The equation is balanced
    #
    # eq = "2 C00027 <=> 2 C00001 + 1 C00007 + 1 C00008"
    # reactants, products, react_ele, prod_ele = CBRdb.get_elements_from_eq(eq, data_c)
    # result = CBRdb.check_eq_unbalanced(react_ele, prod_ele)
    # assert result == True  # The equation is not balanced
    pass


def test_missing_formulas():
    # Check if there is a missing compound ID in the equation
    data_c = pd.read_csv(os.path.abspath("../data/kegg_data_C.csv"))
    eq = "2 C00027 <=> 2 C00001 + 1 C00007"
    assert CBRdb.check_missing_formulas(eq, data_c) == False
    eq = "2 C00027 <=> 2 C00001 + 1 C00007 + 1 C99998"
    assert CBRdb.check_missing_formulas(eq, data_c) == True


def test_strip_ionic():
    print(flush=True)
    data_c = pd.read_csv(os.path.abspath("../data/kegg_data_C.csv"))
    eq = "C05359 <=> C99999"
    reactants, products, react_ele, prod_ele = CBRdb.get_elements_from_eq(eq, data_c, strip_ionic=False)
    print(reactants, flush=True)
    print(products, flush=True)
    print(react_ele, flush=True)
    print(prod_ele, flush=True)
    # there should be a star in the reactants
    assert "*" in react_ele
    # there should be a star in the products
    assert '-' in prod_ele
    print(flush=True)

    reactants, products, react_ele, prod_ele = CBRdb.get_elements_from_eq(eq, data_c, strip_ionic=True)
    print(reactants, flush=True)
    print(products, flush=True)
    print(react_ele, flush=True)
    print(prod_ele, flush=True)
    # there should be no star in the reactants
    assert "*" not in react_ele
    # there should be no star in the products
    assert '-2' not in prod_ele
    print(flush=True)

    eq = "C18091 + C00007 + C01847 <=> C00084 + C00088 + C00061 + C00001"  # R00025
    reactants, products, react_ele, prod_ele = CBRdb.get_elements_from_eq(eq, data_c, strip_ionic=True)
    print(reactants, flush=True)
    print(products, flush=True)
    print(react_ele, flush=True)
    print(prod_ele, flush=True)
    assert "-" not in react_ele


def test_missing_elements():
    print(flush=True)
    # Check if there is a missing element in the equation
    data_c = pd.read_csv(os.path.abspath("../data/kegg_data_C.csv"))
    eq = "2 C00027 <=> 2 C00001 + 1 C00007"
    eq = CBRdb.standardise_eq(eq)
    reactants, products, react_ele, prod_ele = CBRdb.get_elements_from_eq(eq, data_c)
    res = CBRdb.check_full_missing_elements(eq, data_c)
    assert CBRdb.check_missing_elements(react_ele, prod_ele) == False  # The equation is balanced
    assert res == False

    eq = "2 C00027 <=> 2 C00001 + 1 C00007 + 1 C99999 + 1 C05359"
    eq = CBRdb.standardise_eq(eq)
    reactants, products, react_ele, prod_ele = CBRdb.get_elements_from_eq(eq, data_c)
    res = CBRdb.check_full_missing_elements(eq, data_c)
    assert CBRdb.check_missing_elements(react_ele, prod_ele) == True  # C99999 is missing
    assert res == True
    missing_in_react, missing_in_prod = CBRdb.get_missing_elements(react_ele, prod_ele)
    print("Missing in reactants:        ", missing_in_react, flush=True)
    print("Missing in products:         ", missing_in_prod, flush=True)

    eq = "2 C19610 + C00027 + 2 C00080 <=> 2 C19611 + 2 C00001"  # R00011
    eq = CBRdb.standardise_eq(eq)

    # check if the equation is balanced
    res = CBRdb.full_check_eq_unbalanced(eq, data_c)
    print(res, flush=True)

    reactants, products, react_ele, prod_ele = CBRdb.get_elements_from_eq(eq, data_c)
    res = CBRdb.check_full_missing_elements(eq, data_c)
    # assert CBRdb.check_missing_elements(react_ele, prod_ele) == True  # C99999 is missing
    # assert res == True
    missing_in_react, missing_in_prod = CBRdb.get_missing_elements(react_ele, prod_ele)
    print("Missing in reactants:        ", missing_in_react, flush=True)
    print("Missing in products:         ", missing_in_prod, flush=True)


def test_get_small_compounds():
    print(flush=True)
    # Function that loads the small compounds
    small_1 = CBRdb.get_small_compounds(n=1)
    assert len(small_1) == 64  # 1: 64  mostly metals
    small_2 = CBRdb.get_small_compounds(n=2)
    assert len(small_2) == 46  # 2: 46  small molecules
    small_3 = CBRdb.get_small_compounds(n=3)
    assert len(small_3) == 60  # 3: 60  medium molecules

    for item in small_1.values:
        print(item, flush=True)

    small_3_all = CBRdb.get_small_compounds_all(n=3)
    assert len(small_3_all) == 60 + 46 + 64


def test_rebalance_eq():
    print(flush=True)
    data_c = pd.read_csv(os.path.abspath("../data/kegg_data_C.csv"))

    # Attempt to rebalance the equation when there is no need to rebalance
    eq = CBRdb.standardise_eq("2 C00089 <=> C00031 + C03661")
    # Rebalance the equation
    eq_out = CBRdb.rebalance_eq(eq, data_c)
    assert eq_out == eq

    # Attempt to rebalance the equation when there is a need to rebalance
    eq = CBRdb.standardise_eq("4 C00089 <=> C00031 + C03661")
    # Rebalance the equation
    eq_out = CBRdb.rebalance_eq(eq, data_c)
    assert eq_out == "2 C00089 <=> 1 C00031 + 1 C03661"

    # Attempt to rebalance the equation when there is repeated compounds
    eq = CBRdb.standardise_eq("2 C00089 + 2 C00126 <=> C00031 + C03661 + 2 C00125")
    # Rebalance the equation
    eq_out = CBRdb.rebalance_eq(eq, data_c)
    assert eq_out == eq

    # Attempt to rebalance the equation when it is impossible to balance
    eq = CBRdb.standardise_eq("2 C00089 + 2 C00126 <=> C00031 + C03661 + 2 C00125 + 1 C99999")
    # Rebalance the equation
    eq_out = CBRdb.rebalance_eq(eq, data_c)
    assert eq_out == False


def test_get_compounds_with_elements():
    print(flush=True)
    data_c = pd.read_csv(os.path.abspath("../data/kegg_data_C.csv"))
    # Rebalancer would fail on this equation

    eq = CBRdb.standardise_eq("1 C00027 + 2 C00126 <=> 2 C00001 + 2 C00125")

    reactants, products, react_ele, prod_ele = CBRdb.get_elements_from_eq(eq, data_c)
    diff_ele_react, diff_ele_prod = CBRdb.compare_dict_values(react_ele, prod_ele)
    # Get the set of keys in react_ele and prod_ele
    element_symbols = list(set(diff_ele_react.keys()).union(set(diff_ele_prod.keys())))
    assert element_symbols == ['H']


def test_inject_compounds():
    print(flush=True)
    data_c_1 = CBRdb.get_small_compounds(n=1)

    data_c = pd.read_csv(os.path.abspath("../data/kegg_data_C.csv"))
    # Rebalancer would fail on this equation
    eq = CBRdb.standardise_eq("1 C00027 + 2 C00126 <=> 2 C00001 + 2 C00125")

    reactants, products, react_ele, prod_ele = CBRdb.get_elements_from_eq(eq, data_c)
    # Get the difference in the elements population
    diff_ele_react, diff_ele_prod = CBRdb.compare_dict_values(react_ele, prod_ele)
    print("Differences in reactants:    ", diff_ele_react, flush=True)
    print("Differences in products:     ", diff_ele_prod, flush=True)
    assert diff_ele_react == {'H': 90}
    assert diff_ele_prod == {'H': 92}
    # Get the difference in the elements population
    missing_in_react, missing_in_prod = CBRdb.get_missing_elements(react_ele, prod_ele)
    print("Missing in reactants:        ", missing_in_react, flush=True)
    print("Missing in products:         ", missing_in_prod, flush=True)

    # Get the eq information
    reactants, products, react_ele, prod_ele = CBRdb.get_elements_from_eq(eq, data_c)
    # Get the difference in the elements
    diff_ele_react, diff_ele_prod = CBRdb.compare_dict_values(react_ele, prod_ele)

    # Get the compounds that might match
    compounds = CBRdb.get_compounds_with_matching_elements(data_c_1, diff_ele_react, diff_ele_prod)

    print("Compounds that might match:  ", compounds, flush=True)

    # Inject the compounds into the equation
    eq_out = CBRdb.fix_imbalance_core(eq, diff_ele_react, diff_ele_prod, compounds[0])
    # Rebalance the equation
    eq_out = CBRdb.rebalance_eq(eq_out, data_c)

    assert eq_out == "1 C00027 + 2 C00080 + 2 C00126 <=> 2 C00001 + 2 C00125"


def test_kitchen_sink():
    print(flush=True)
    data_c_1 = CBRdb.get_small_compounds(n=1)
    data_c = pd.read_csv(os.path.abspath("../data/kegg_data_C.csv"))
    # Rebalancer would fail on this equation
    eq = CBRdb.standardise_eq("1 C00027 + 2 C00126 <=> 2 C00001 + 2 C00125")
    eq_out = CBRdb.kitchen_sink(eq, data_c, data_c_1, None)
    print(eq_out, flush=True)
    assert eq_out == "1 C00027 + 2 C00080 + 2 C00126 <=> 2 C00001 + 2 C00125"


def test_star_problem():
    print(flush=True)
    data_c = pd.read_csv(os.path.abspath("../data/kegg_data_C.csv"))
    # Rebalancer would fail on this equation
    eq = CBRdb.standardise_eq("1 C00001 + 1 C00454 <=> 1 C00009 + 1 C00215")
    print(eq, flush=True)
    reactants, products, react_ele, prod_ele = CBRdb.get_elements_from_eq(eq, data_c)
    print(reactants, flush=True)
    print(products, flush=True)
    print(react_ele, flush=True)
    print(prod_ele, flush=True)

    assert CBRdb.dict_ele_contains_star(react_ele, prod_ele)


def test_dict_values_problem():
    print(flush=True)
    data_c = pd.read_csv(os.path.abspath("../data/kegg_data_C.csv"))
    # Rebalancer would fail on this equation
    eq = CBRdb.standardise_eq("C00011 <=> 1 C07728")
    print(eq, flush=True)
    reactants, products, react_ele, prod_ele = CBRdb.get_elements_from_eq(eq, data_c)
    diff_ele_react, diff_ele_prod = CBRdb.compare_dict_values(react_ele, prod_ele)

    print(reactants, flush=True)
    print(products, flush=True)
    print(react_ele, flush=True)
    print(prod_ele, flush=True)
    print(diff_ele_react, flush=True)
    print(diff_ele_prod, flush=True)


def test_rebalancer_fail():
    print(flush=True)
    data_c_1 = CBRdb.get_small_compounds(n=1)
    data_c = pd.read_csv(os.path.abspath("../data/kegg_data_C.csv"))
    # Rebalancer would fail on this equation
    eq = CBRdb.standardise_eq("1 C00001 + 1 C05195 <=> 1 C19684")
    eq_out = CBRdb.kitchen_sink(eq, data_c, data_c_1, None)
    print(eq_out, flush=True)


def test_mol_replacer_smi():
    print(flush=True)

    smi = "*CC(=O)c1ccc(C(=O)O)cc1"
    smi_out = CBRdb.mol_replacer_smi(smi, target="[H]")
    print(smi_out, flush=True)
    assert smi_out == "[H]CC(=O)c1ccc(C(=O)O)cc1"


def test_mol_replacer():
    print(flush=True)
    smi = "*CC(=O)c1ccc(C(=O)O)cc1"
    mol = Chem.MolFromSmiles(smi)
    mol_out = CBRdb.mol_replacer(mol, target="[H]")
    smi_out = Chem.MolToSmiles(mol_out)
    inchi_out = Chem.MolToInchi(mol_out)
    print(smi_out, flush=True)
    print(inchi_out, flush=True)
    assert smi_out == "[H]CC(=O)c1ccc(C(=O)O)cc1"


def test_plot_eq_line():
    print(flush=True)
    eq = "1 C00001 + 1 C00454 <=> 1 C00009 + 1 C00215"
    data_c = pd.read_csv(os.path.abspath("../data/kegg_data_C.csv"))
    CBRdb.plot_eq_line(eq, data_c)
    pass


def test_plot_reaction_id():
    print(flush=True)
    data_c = pd.read_csv(os.path.abspath("../data/kegg_data_C.csv"))
    data_r = pd.read_csv(os.path.abspath("../data/kegg_data_R.csv"))
    CBRdb.plot_reaction_id('R00001', data_r, data_c)
    pass


def test_calculate_ccsd_energy():
    print(flush=True)
    smi = "O"  # SMILES representation of the molecule (water in this case)
    atoms, charge, multiplicity = CBRdb.smi_to_atoms(smi)
    energy = CBRdb.calculate_ccsd_energy(atoms,
                                         charge=charge,
                                         multiplicity=multiplicity)
    print(energy, flush=True)
    ref_energy = -2077.148240270791  # Reference Gibbs free energy value
    assert np.allclose(energy, ref_energy,
                       atol=1e-1), f"Calculated energy {energy} does not match reference {ref_energy}"


def test_calculate_free_energy():
    print(flush=True)
    smi = "O"  # SMILES representation of the molecule (water in this case)
    atoms, charge, multiplicity = CBRdb.smi_to_atoms(smi)
    energy, enthalpy, entropy = CBRdb.calculate_free_energy(atoms, charge=charge, multiplicity=multiplicity)
    print(energy, enthalpy, entropy, flush=True)
    ref_energy = -2079.599879755067  # Reference Gibbs free energy value
    assert np.allclose(energy, ref_energy,
                       atol=1e-3), f"Calculated energy {energy} does not match reference {ref_energy}"
    ref_enthalpy = -2079.017167243439
    assert np.allclose(enthalpy, ref_enthalpy,
                       atol=1e-3), f"Calculated enthalpy {enthalpy} does not match reference {ref_enthalpy}"
    ref_entropy = -0.5827125116277472
    assert np.allclose(entropy, ref_entropy,
                       atol=1e-3), f"Calculated entropy {entropy} does not match reference {ref_entropy}"

    energy, enthalpy, entropy = CBRdb.calculate_free_energy(atoms,
                                                            charge=charge,
                                                            multiplicity=multiplicity,
                                                            use_ccsd=True)
    print(energy, enthalpy, entropy, flush=True)
    ref_energy = -2077.1273362752813  # Reference Gibbs free energy value
    assert np.allclose(energy, ref_energy,
                       atol=1e-3), f"Calculated energy {energy} does not match reference {ref_energy}"
    ref_enthalpy = -2076.5446237636534
    assert np.allclose(enthalpy, ref_enthalpy,
                       atol=1e-3), f"Calculated enthalpy {enthalpy} does not match reference {ref_enthalpy}"
    ref_entropy = -0.5827125116277472
    assert np.allclose(entropy, ref_entropy,
                       atol=1e-3), f"Calculated entropy {entropy} does not match reference {ref_entropy}"

    energy, enthalpy, entropy = CBRdb.calculate_free_energy(atoms,
                                                            charge=charge,
                                                            multiplicity=multiplicity,
                                                            use_ccsd=True,
                                                            f_solv=True)
    print(energy, enthalpy, entropy, flush=True)
    ref_energy = -2077.0719458992658  # Reference Gibbs free energy value
    assert np.allclose(energy, ref_energy,
                       atol=1e-3), f"Calculated energy {energy} does not match reference {ref_energy}"
    ref_enthalpy = -2076.4889517497927
    assert np.allclose(enthalpy, ref_enthalpy,
                       atol=1e-3), f"Calculated enthalpy {enthalpy} does not match reference {ref_enthalpy}"
    ref_entropy = -0.5829941494730994
    assert np.allclose(entropy, ref_entropy,
                       atol=1e-3), f"Calculated entropy {entropy} does not match reference {ref_entropy}"


def test_optimise_atoms():
    print(flush=True)
    smi = "O"  # SMILES representation of the molecule (water in this case)
    atoms, charge, multiplicity = CBRdb.smi_to_atoms(smi)
    atoms = CBRdb.optimise_atoms(atoms, charge=charge, multiplicity=multiplicity)
    pos = atoms.get_positions()
    print("Optimized positions:", pos, flush=True)
    ref_pos = [[0.02628414, -2.70228834, 4.02632285],
               [0.83770207, -3.11418435, 3.72496985],
               [-0.50718621, -3.43022732, 4.3494073]]
    assert np.allclose(atoms.get_positions(), ref_pos)


def test_calculate_vib_spectrum():
    print(flush=True)
    smi = "O"  # SMILES representation of the molecule (water in this case)
    atoms, charge, multiplicity = CBRdb.smi_to_atoms(smi)

    data_ir, data_raman, data_vib = CBRdb.calculate_vib_spectrum(atoms,
                                                                 charge=charge,
                                                                 multiplicity=multiplicity)
    print(data_ir, flush=True)
    print(data_raman, flush=True)
    print(data_vib, flush=True)


def test_calculate_goat():
    print(flush=True)
    smi = "OCCCCC"  # SMILES representation of the molecule (water in this case)
    atoms, charge, multiplicity = CBRdb.smi_to_atoms(smi)

    atoms, data_goat = CBRdb.calculate_goat(atoms,
                                            charge=charge,
                                            multiplicity=multiplicity)

    print(data_goat, flush=True)

    view(atoms)


def test_calculate_free_energy_batch():
    print(flush=True)
    # t_list = [300, 400]  # List of temperatures in Kelvin
    # p_list = [1.0, 2.0]  # List of pressures in atm
    # smi = "OCCCCC"  # SMILES representation of the molecule (water in this case)
    # atoms, charge, multiplicity = CBRdb.smi_to_atoms(smi)
    # CBRdb.calculate_free_energy_batch(atoms, t_list, p_list, charge, multiplicity)
    atoms = read('data/orca.xyz')  # Load atoms from an XYZ file
    hessian = os.path.join(os.getcwd(), 'data/orca.hess')  # Load hessian from a numpy file
    energy, enthalpy, entropy = CBRdb.calculate_free_energy_batch(atoms, hessian, 400, 1.0)
    print(energy, enthalpy, entropy, flush=True)
    # 3.2295803570495636 4.874475172217753 -1.64489481516819
    # Values messed up
    ref_energy = -2077.0719458992658  # Reference Gibbs free energy value
    assert np.allclose(energy, ref_energy,
                       atol=1e-3), f"Calculated energy {energy} does not match reference {ref_energy}"
    ref_enthalpy = -2076.4889517497927
    assert np.allclose(enthalpy, ref_enthalpy,
                       atol=1e-3), f"Calculated enthalpy {enthalpy} does not match reference {ref_enthalpy}"
    ref_entropy = -0.5829941494730994
    assert np.allclose(entropy, ref_entropy,
                       atol=1e-3), f"Calculated entropy {entropy} does not match reference {ref_entropy}"


def test_to_smarts_rxn_line():
    """
    Test the `to_smarts_rxn_line` function from the `CBRdb` module.

    This test verifies that the function correctly converts a chemical equation
    into a SMARTS reaction line.

    Steps:
    1. Load compound data from a CSV file.
    2. Define a chemical equation with reactants and products.
    3. Convert the equation into a SMARTS reaction line using `to_smarts_rxn_line`.
    4. Assert that the generated SMARTS reaction line matches the expected value.

    Raises:
    -------
    AssertionError:
        If the generated SMARTS reaction line does not match the expected value.
    """
    print(flush=True)
    data_c = pd.read_csv(os.path.abspath("../data/kegg_data_C.csv"))  # Load compound data
    eq = "2 C19610 + C00027 + 2 C00080 <=> 2 C19611 + 2 C00001"  # Define chemical equation

    r_smarts = CBRdb.to_smarts_rxn_line(eq, data_c, add_stoich=True)  # Convert equation to SMARTS reaction line
    print(r_smarts, flush=True)
    assert r_smarts == '2[Mn+2].1OO.2[H+]>>2[Mn+3].2[H]O[H]'  # Verify the result

    r_smarts = CBRdb.to_smarts_rxn_line(eq, data_c, add_stoich=False)  # Convert equation to SMARTS reaction line
    print(r_smarts, flush=True)
    assert r_smarts == '[Mn+2].OO.[H+]>>[Mn+3].[H]O[H]'  # Verify the result


def test_get_properties():
    smis_list = [
        "O",  # Water
        "CC(=O)O",  # Acetic acid
        "C1=CC=CC=C1",  # Benzene
        "C1CCCCC1",  # Cyclohexane
        "C1=CC=C(C=C1)C(=O)O"  # Benzoic acid
    ]
    mols = [(Chem.MolFromSmiles(smi)) for smi in smis_list]

    # Get the properties
    properties = CBRdb.get_properties(mols)
    print(properties)


def get_formation_references(mol):
    from collections import defaultdict

    # Get elemental composition
    atom_counts = defaultdict(int)
    symbols = set()
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        atom_counts[symbol] += 1.0
        symbols.add(symbol)

    supported_set = {'H', 'Ba', 'Cu', 'S', 'K', 'Cl', 'Rb', 'B', 'N', 'Se', 'Te', 'O', 'Fe', 'Co', 'Mg', 'Ge', 'I',
                     'Tl', 'Pt', 'Xe', 'Zn', 'Gd', 'Cd', 'C', 'Al', 'F', 'Li', 'Ca', 'Ni', 'Th', 'Sr', 'Sn', 'Au', 'Ag',
                     'V', 'Pu', 'Sb', 'Mn', 'Cr', 'Mo', 'Ra', 'Pb', 'Bi', 'As', 'Hg', 'Si', 'Br', 'Rn', 'P', 'W', 'Be',
                     'Na'}
    if not symbols.issubset(supported_set):
        raise ValueError(f"Unsupported elements in the molecule: {symbols - supported_set}")

    # Standard reference molecules
    references = []
    # loop over the atom_counts dictionary
    for symbol, count in atom_counts.items():
        if count > 0:
            if symbol == 'H':
                references.append(("[H][H]", atom_counts['H'] / 2.0))
            elif symbol == 'Cl':
                references.append(("Cl-Cl", atom_counts['Cl'] / 2.0))
            elif symbol == 'N':
                references.append(("N#N", atom_counts['N'] / 2.0))
            elif symbol == 'O':
                references.append(("O=O", atom_counts['O'] / 2.0))
            elif symbol == 'I':
                references.append(("I-I", atom_counts['I'] / 2.0))
            elif symbol == 'F':
                references.append(("F-F", atom_counts['F'] / 2.0))
            elif symbol == 'Br':
                references.append(("Br-Br", atom_counts['Br'] / 2.0))
            else:
                references.append((f"{symbol}", atom_counts[f"{symbol}"]))

    return references


def test_free_energy_formation():
    print(flush=True)
    smi = "OO"
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)

    atoms, charge, multiplicity = CBRdb.smi_to_atoms(smi)
    energy, enthalpy, entropy = CBRdb.calculate_free_energy(atoms,
                                                            charge=charge,
                                                            multiplicity=multiplicity,
                                                            xc='pbe')
    print(energy, flush=True)

    energy_atoms = 0.0
    references = get_formation_references(mol)
    for ref_smi, ref_count in references:
        ref_atoms, ref_charge, ref_multiplicity = CBRdb.smi_to_atoms(ref_smi)
        ref_energy, _, _ = CBRdb.calculate_free_energy(ref_atoms,
                                                       charge=ref_charge,
                                                       multiplicity=ref_multiplicity,
                                                       xc='pbe')
        print(f"Reference: {ref_smi}, Count: {ref_count}, Energy: {ref_energy}", flush=True)
        energy_atoms += ref_energy * ref_count
    print(energy_atoms, flush=True)
