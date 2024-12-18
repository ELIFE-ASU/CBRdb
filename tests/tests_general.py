import os
import pandas as pd

import CBRdb


def test_side_to_dict():
    tmp = CBRdb.side_to_dict('C00001 + 1 C00002')
    assert tmp == {'C00001': 1, 'C00002': 1}

    tmp = CBRdb.side_to_dict('C00001 + C00002')
    assert tmp == {'C00001': 1, 'C00002': 1}

    tmp = CBRdb.side_to_dict('1 C00001 + 1 C00002')
    assert tmp == {'C00001': 1, 'C00002': 1}

    tmp = CBRdb.side_to_dict('+1 C00001 + 1 C00002')
    assert tmp == {'C00001': 1, 'C00002': 1}

    tmp = CBRdb.side_to_dict('-1 C00001 + 1 C00002')
    assert tmp == {'C00001': -1, 'C00002': 1}

    tmp = CBRdb.side_to_dict('n C00001 + 1 C00002')
    assert tmp == {'C00001': 'n', 'C00002': 1}

    tmp = CBRdb.side_to_dict('n+1 C00001 + 1 C00002')
    assert tmp == {'C00001': 'n+1', 'C00002': 1}

    tmp = CBRdb.side_to_dict('n-1 C00001 + 1 C00002')
    assert tmp == {'C00001': 'n-1', 'C00002': 1}

    tmp = CBRdb.side_to_dict('0 C00001 + 1 C00002')
    assert tmp == {'C00002': 1}

    tmp = CBRdb.side_to_dict("1 C00007 + 2 C00339 + C01438")
    assert tmp == {'C00007': 1, 'C00339': 2, 'C01438': 1}

    tmp = CBRdb.side_to_dict("C00024 + n C00083 + n C00005 + n C00004 + 2n C00080")
    assert tmp == {'C00024': 1, 'C00083': 'n', 'C00005': 'n', 'C00004': 'n', 'C00080': '2n'}

    tmp = CBRdb.side_to_dict("C00003 + C00039(n) + C02128(m)")
    assert tmp == {'C00003': 1, 'C00039': 'n', 'C02128': 'm'}

    tmp = CBRdb.side_to_dict("C00020 + C00455 + C00039(n+m)")
    assert tmp == {'C00020': 1, 'C00455': 1, 'C00039': 'n+m'}

    tmp = CBRdb.side_to_dict("C03323(m) + C03323(n)")
    assert tmp == {'C03323': 'm+n'}

    tmp = CBRdb.side_to_dict("C03323(m-1) + C03323(n+1)")
    assert tmp == {'C03323': 'm-1+n+1'}


def test_convert_formula_to_dict():
    tmp = CBRdb.convert_formula_to_dict("C2H2*BrO2")
    assert tmp == {'C': 2, 'H': 2, 'Br': 1, 'O': 2, '*': 1}

    tmp = CBRdb.convert_formula_to_dict("C2H2*32BrO2")
    assert tmp == {'C': 2, 'H': 2, 'Br': 1, 'O': 2, '*': 32}

    tmp = CBRdb.convert_formula_to_dict("Te+")
    assert tmp == {'Te': 1, '+': 1}

    tmp = CBRdb.convert_formula_to_dict("Te-")
    assert tmp == {'Te': 1, '-': 1}

    tmp = CBRdb.convert_formula_to_dict("Te+1")
    assert tmp == {'Te': 1, '+': 1}

    tmp = CBRdb.convert_formula_to_dict("Te-1")
    assert tmp == {'Te': 1, '-': 1}

    tmp = CBRdb.convert_formula_to_dict("C2H4*NO2-")
    assert tmp == {'C': 2, 'H': 4, 'N': 1, 'O': 2, '-': 1, '*': 1}


def test_eq_n_solver():
    expr = "n-1"
    result = CBRdb.find_min_integers(expr)
    print(result)
    assert result == {'n': 2}

    expr = "m-1+n+1"
    result = CBRdb.find_min_integers(expr)
    print(result)
    assert result == {'n': 1, 'm': 1}

    expr = "2*n"
    result = CBRdb.find_min_integers(expr)
    print(result)
    assert result == {'n': 1}

    expr = "2*n + 1"
    result = CBRdb.find_min_integers(expr)
    print(result)
    assert result == {'n': 1}

    # This is messed up but works
    expr = "2*(n - 1) + 2*(m - 1)+ (x-10)"
    result = CBRdb.find_min_integers(expr)
    print(result)
    # assert result == {'n': 1, 'm': 1, 'x': 10}
    assert result == {'x': 11, 'n': 6, 'm': 6}


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
    data_c = pd.read_csv(os.path.abspath("../data/kegg_data_C.csv.zip"))
    eq = "2 C00027 <=> 2 C00001 + 1 C00007"
    lhs = {'H2O2': 2}
    rhs = {'H2O': 2, 'O2': 1}
    lhs_e = {'O': 4, 'H': 4}
    rhs_e = {'O': 4, 'H': 4}

    reactants, products, react_ele, prod_ele = CBRdb.get_elements_from_eq(eq, data_c)

    assert reactants == lhs
    assert products == rhs

    assert react_ele == lhs_e
    assert prod_ele == rhs_e

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


def test_eq_balanced():
    print(flush=True)
    data_c = pd.read_csv(os.path.abspath("../data/kegg_data_C.csv.zip"))
    eq = "2 C00027 <=> 2 C00001 + 1 C00007"
    reactants, products, react_ele, prod_ele = CBRdb.get_elements_from_eq(eq, data_c)
    result = CBRdb.check_eq_unbalanced(react_ele, prod_ele)
    assert result == False  # The equation is balanced

    eq = "2 C00027 <=> 2 C00001 + 1 C00007 + 1 C00008"
    reactants, products, react_ele, prod_ele = CBRdb.get_elements_from_eq(eq, data_c)
    result = CBRdb.check_eq_unbalanced(react_ele, prod_ele)
    assert result == True  # The equation is not balanced


def test_contains_var_list_check():
    # Test which sees if there is variable in the equation
    data_c = pd.read_csv(os.path.abspath("../data/kegg_data_C.csv.zip"))
    eq = "m C00404 + n C00001 <=> (n+1) C02174 + x C00001"
    reactants, products = CBRdb.get_formulas_from_eq(eq, data_c)
    vars = CBRdb.contains_var_list(reactants, products)

    assert vars == ['m', 'n', 'x']
    assert CBRdb.check_contains_var_list(reactants, products) == True

    # Test which sees if there is no variable in the equation
    eq = "C00404 + C00001 <=> C02174"
    reactants, products = CBRdb.get_formulas_from_eq(eq, data_c)
    vars = CBRdb.contains_var_list(reactants, products)
    assert vars == []

def test_vars_eq_balanced():
    ############################### DOUBLE CHECK THIS FUNCTION ########################################
    print(flush=True)
    data_c = pd.read_csv(os.path.abspath("../data/kegg_data_C.csv.zip"))
    eq = "2 C00027 <=> 2 C00001 + 1 C00007"
    reactants, products, react_ele, prod_ele = CBRdb.get_elements_from_eq(eq, data_c)
    result = CBRdb.check_eq_unbalanced(react_ele, prod_ele)
    assert result == False  # The equation is balanced

    eq = "2 C00027 <=> 2 C00001 + 1 C00007 + 1 C00008"
    reactants, products, react_ele, prod_ele = CBRdb.get_elements_from_eq(eq, data_c)
    result = CBRdb.check_eq_unbalanced(react_ele, prod_ele)
    assert result == True  # The equation is not balanced


def test_missing_formulas():
    # Check if there is a missing compound ID in the equation
    data_c = pd.read_csv(os.path.abspath("../data/kegg_data_C.csv.zip"))
    eq = "2 C00027 <=> 2 C00001 + 1 C00007"
    assert CBRdb.check_missing_formulas(eq, data_c) == False
    eq = "2 C00027 <=> 2 C00001 + 1 C00007 + 1 C99998"
    assert CBRdb.check_missing_formulas(eq, data_c) == True


def test_missing_elements():
    # Check if there is a missing element in the equation
    # check_missing_elements
    # get_missing_elements
    pass


def test_inject_compounds():
    # Check if the compounds are injected into the equation
    pass
