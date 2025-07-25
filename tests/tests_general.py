import os

import pandas as pd
from rdkit import Chem as Chem
import numpy as np
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


def test_calculate_free_energy():
    """
    Test the `calculate_free_energy` function from the `CBRdb` module.

    This test verifies that the calculated Gibbs free energy for a given molecule
    matches the expected reference value within a specified tolerance.

    Steps:
    1. Create a molecule object from its SMILES representation.
    2. Calculate the Gibbs free energy using the `calculate_free_energy` function.
    3. Compare the calculated energy to the reference value using `np.allclose`.

    Raises:
    -------
    AssertionError:
        If the calculated energy does not match the reference energy within the specified tolerance.
    """
    print(flush=True)
    smi = "O"  # SMILES representation of the molecule (water in this case)
    mol = Chem.MolFromSmiles(smi)  # Convert SMILES to an RDKit molecule object
    energy = CBRdb.calculate_free_energy(mol)  # Calculate the Gibbs free energy
    print(energy, flush=True)
    ref_energy = -2079.975658927851  # Reference Gibbs free energy value
    assert np.allclose(energy, ref_energy,
                       atol=1e-1), f"Calculated energy {energy} does not match reference {ref_energy}"


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

    r_smarts = CBRdb.to_smarts_rxn_line(eq, data_c)  # Convert equation to SMARTS reaction line
    print(r_smarts, flush=True)
    assert r_smarts == '2[Mn+2].1OO.2[H+]>>2[Mn+3].2[H]O[H]'  # Verify the result
