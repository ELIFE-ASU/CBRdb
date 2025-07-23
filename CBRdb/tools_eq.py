import copy
import os
import re

import chemparse
import sympy as sp
from chempy import balance_stoichiometry
from rdkit import Chem
from rdkit.Chem import Draw


def strip_ionic_states(formula):
    """
    Removes trailing ionic states from a chemical formula.

    Parameters:
    formula (str): The chemical formula from which to remove ionic states.

    Returns:
    str: The chemical formula without trailing ionic states.
    """
    # Regex to match trailing ionic states (+, -, +2, -4, etc.)
    return re.sub(r'[+-]\d*$', '', formula)


def convert_formula_to_dict(formula, strip_ionic=True):
    """
    Converts a chemical formula into a dictionary with element symbols as keys and their counts as values.

    Parameters:
    formula (str): The chemical formula to be converted.
    strip_ionic (bool, optional): Flag to indicate whether to remove ionic states from the formula. Default is False.

    Returns:
    dict: A dictionary where keys are element symbols and values are their counts in the formula.
    """

    # Count the occurrences of '*' followed by an optional number
    star_count = sum(int(num) if num else 1 for num in re.findall(r'\*(\d*)', formula))

    # Replace '*' followed by an optional number with an empty string
    formula_new = re.sub(r'\*\d*', '', formula)

    # Remove the ionic states
    if strip_ionic:
        formula_new = strip_ionic_states(formula_new)

    # Parse the formula into a dictionary
    formula_dict = {k: int(v) for k, v in chemparse.parse_formula(formula_new).items()}

    # Add the '*' count to the dictionary if there were any '*'
    if star_count > 0:
        formula_dict['*'] = star_count

    return formula_dict


def multiply_dict(input_dict, multiplier):
    """
    Multiplies each value in the input dictionary by the given multiplier.

    Parameters:
    input_dict (dict): The dictionary with values to be multiplied.
    multiplier (int or float): The number by which to multiply each value in the dictionary.

    Returns:
    dict: A new dictionary with the same keys as input_dict and values multiplied by the multiplier.
    """
    return {key: value * multiplier for key, value in input_dict.items()}


def add_dicts(dict1, dict2):
    """
    Adds the values of two dictionaries for matching keys.

    Parameters:
    dict1 (dict): The first dictionary.
    dict2 (dict): The second dictionary.

    Returns:
    dict: A new dictionary with keys from both input dictionaries.
          The values are the sum of values from dict1 and dict2 for matching keys.
    """
    return {key: dict1.get(key, 0) + dict2.get(key, 0) for key in set(dict1) | set(dict2)}


def subtract_dicts(dict1, dict2):
    """
    Subtracts the values of the second dictionary from the first dictionary for matching keys.

    Parameters:
    dict1 (dict): The first dictionary.
    dict2 (dict): The second dictionary.

    Returns:
    dict: A new dictionary with keys from both input dictionaries.
          The values are the result of subtracting values in dict2 from dict1 for matching keys.
    """
    return {key: dict1.get(key, 0) - dict2.get(key, 0) for key in set(dict1) | set(dict2)}


def compare_dict_keys(dict1, dict2):
    """
    Compares the keys of two dictionaries and identifies keys that are missing in each dictionary.

    Parameters:
    dict1 (dict): The first dictionary.
    dict2 (dict): The second dictionary.

    Returns:
    tuple: A tuple containing two sets:
           - The first set contains keys that are in dict1 but not in dict2.
           - The second set contains keys that are in dict2 but not in dict1.
    """
    missing_in_dict2 = set(dict1) - set(dict2)
    missing_in_dict1 = set(dict2) - set(dict1)
    return missing_in_dict2, missing_in_dict1


def compare_dicts(dict1, dict2):
    """
    Compares two dictionaries to check if they have the same keys and corresponding values.

    Parameters:
    dict1 (dict): The first dictionary to compare.
    dict2 (dict): The second dictionary to compare.

    Returns:
    bool: True if both dictionaries have the same keys and corresponding values, False otherwise.
    """
    return dict1.keys() == dict2.keys() and all(dict1[key] == dict2[key] for key in dict1)


def compare_dict_values(dict1, dict2):
    """
    Compares the values of two dictionaries and identifies differences.
    If any values are None, they are treated as 0.

    Parameters:
    dict1 (dict): The first dictionary.
    dict2 (dict): The second dictionary.

    Returns:
    tuple: A tuple containing two dictionaries:
           - The first dictionary contains key-value pairs from dict1 where the values differ from dict2.
           - The second dictionary contains key-value pairs from dict2 where the values differ from dict1.
    """
    diff_in_dict1 = sort_dict_by_keys({key: (dict1.get(key) or 0) for key in set(dict1) | set(dict2) if
                                       (dict1.get(key) or 0) != (dict2.get(key) or 0)})
    diff_in_dict2 = sort_dict_by_keys({key: (dict2.get(key) or 0) for key in set(dict1) | set(dict2) if
                                       (dict1.get(key) or 0) != (dict2.get(key) or 0)})
    return diff_in_dict1, diff_in_dict2


def get_formulas_from_ids(ids, c_data):
    """
    Retrieves the chemical formulas corresponding to a list of compound IDs from a DataFrame.

    Parameters:
    ids (list): A list of compound IDs.
    c_data (DataFrame): A pandas DataFrame containing compound data with 'compound_id' and 'formula' columns.

    Returns:
    list: A list of chemical formulas corresponding to the given compound IDs.
    """
    return c_data.loc[c_data["compound_id"].isin(ids), "formula"].tolist()


def generate_compound_dict(c_data):
    """
    Generates a dictionary from the compound_id, with compound_id as the key and formula as the value.

    Parameters:
    c_data (DataFrame): A pandas DataFrame containing compound data with 'compound_id' and 'formula' columns.

    Returns:
    dict: A dictionary with compound_id as keys and formulas as values.
    """
    return dict(zip(c_data['compound_id'], c_data['formula']))


def clean_up_eq(eq):
    """
    Cleans up the chemical equation string by reformatting compound identifiers.

    Parameters:
    eq (str): The chemical equation string to be cleaned.

    Returns:
    str: The cleaned chemical equation string with compound identifiers reformatted.
    """
    # Replace 'n+1 C02616' with '(n+1) C02616'
    return re.sub(r'([A-Z]\d{5})\(([^)]+)\)', r'(\2) \1', eq)


def merge_duplicates(matches, coeff_out):
    """
    Merges duplicate entries in the matches list by summing their corresponding coefficients.

    Parameters:
    matches (list): A list of match identifiers.
    coeff_out (list): A list of coefficients corresponding to the matches.

    Returns:
    tuple: A tuple containing two lists:
           - The first list contains unique match identifiers.
           - The second list contains the summed coefficients for each unique match.
    """
    merged_coeff = {}
    for match, coeff in zip(matches, coeff_out):
        if match in merged_coeff:
            merged_coeff[match] = (
                merged_coeff[match] + coeff
                if isinstance(merged_coeff[match], int) and isinstance(coeff, int)
                else f"{merged_coeff[match]}+{coeff}"
            )
        else:
            merged_coeff[match] = coeff
    return list(merged_coeff.keys()), list(merged_coeff.values())


def side_to_dict(s):
    """
    Converts a chemical equation side into a dictionary with compounds as keys and their coefficients as values.

    Parameters:
    s (str): A string representing one side of a chemical equation.

    Returns:
    dict: A dictionary where keys are compound identifiers and values are their coefficients.
    """
    # Clean up the equation string
    s = clean_up_eq(s)

    # Split the string by compound identifiers and strip whitespace
    coeff = [c.strip() for c in re.split(r"[A-Z]\d{5}", s)]

    # Replace "+" with "1" and strip remaining whitespace
    coeff = [c if c != "+" else "1" for c in coeff]
    coeff = [c.replace("+ ", "").strip() for c in coeff]

    # Replace empty strings with 1
    coeff = [1 if c == '' else c for c in coeff]

    coeff_out = []
    for c in coeff:
        try:
            # Convert coefficients to integers if possible
            coeff_out.append(int(c))
        except ValueError:
            # Strip parentheses if conversion fails
            coeff_out.append(c.strip('(').strip(')'))
    # try:
    #     # Convert coefficients to integers if possible
    #     coeff_out = [int(c) for c in coeff]
    # except ValueError:
    #     coeff_out = [str(c).strip('(').strip(')') for c in coeff]

    # Find all compound identifiers in the string
    matches = re.findall(r"[A-Z]\d{5}", s)

    # Merge duplicates in matches and their corresponding coefficients
    matches, coeff_out = merge_duplicates(matches, coeff_out)

    # Create a dictionary with compound identifiers as keys and coefficients as values
    out = {k: v for k, v in zip(matches, coeff_out) if v != 0}
    return out


def eq_to_dict(eq):
    """
    Converts a chemical equation string into a tuple of dictionaries representing the reactants and products.

    Parameters:
    eq (str): A string representing a chemical equation, with reactants and products separated by '<=>'.

    Returns:
    tuple: A tuple containing two dictionaries:
           - The first dictionary represents the reactants.
           - The second dictionary represents the products.
    """
    return tuple(map(side_to_dict, eq.split('<=>')))


def dict_to_side(d):
    """
    Converts a dictionary of molecules and their counts into a string representation for one side of a chemical equation.

    Parameters:
    d (dict): A dictionary where keys are molecule identifiers and values are their counts.

    Returns:
    str: A string representation of the dictionary in the format 'count molecule + count molecule + ...'.
    """
    return ' + '.join(f"{v} {k}" for k, v in d.items())


def dicts_to_eq(reactants, products):
    """
    Converts dictionaries of reactants and products into a string representation of a chemical equation.

    Parameters:
    reactants (dict): A dictionary where keys are molecule identifiers and values are their counts for reactants.
    products (dict): A dictionary where keys are molecule identifiers and values are their counts for products.

    Returns:
    str: A string representation of the chemical equation in the format 'reactants <=> products'.
    """
    return f"{dict_to_side(reactants)} <=> {dict_to_side(products)}"


def strip_plus_x(input_string):
    """
    Removes any occurrences of '+<number>' from the input string and replaces remaining '+' characters with an empty string.

    Parameters:
    input_string (str): The string to be processed.

    Returns:
    str: The processed string with '+<number>' and '+' characters removed.
    """
    return re.sub(r'\+\d+', '', input_string).replace('+', '') if '+' in input_string else input_string


def get_ids_to_formulas(compound_dict, c_data):
    """
    Retrieves the chemical formulas corresponding to the compound IDs in the given dictionary.

    Parameters:
    compound_dict (dict): A dictionary where keys are compound IDs.
    c_data (DataFrame): A pandas DataFrame containing compound data with 'compound_id' and 'formula' columns.

    Returns:
    dict: A dictionary where keys are compound IDs and values are the corresponding chemical formulas with '+<number>' and '+' characters removed.
    """
    comp_dict = generate_compound_dict(c_data)
    return {id: strip_plus_x(comp_dict[id]) for id in compound_dict.keys()}


def convert_ids_to_formulas(in_dict, react_id_form):
    """
    Converts a dictionary of compound IDs to their corresponding chemical formulas.

    Parameters:
    in_dict (dict): A dictionary where keys are compound IDs and values are their counts.
    react_id_form (dict): A dictionary mapping compound IDs to their chemical formulas.

    Returns:
    dict: A dictionary where keys are chemical formulas and values are their counts.
    """
    return {react_id_form[id]: count for id, count in in_dict.items()}


def convert_formulas_to_ids(formulas_dict, react_id_form):
    """
    Converts a dictionary of chemical formulas to their corresponding compound IDs.

    Parameters:
    formulas_dict (dict): A dictionary where keys are chemical formulas and values are their counts.
    react_id_form (dict): A dictionary mapping compound IDs to their chemical formulas.

    Returns:
    dict: A dictionary where keys are compound IDs and values are their counts.
    """
    inverse_react_id_form = {v: k for k, v in react_id_form.items()}
    return {inverse_react_id_form[formula]: count for formula, count in formulas_dict.items()}


def remove_electron_cid(compound_dict):
    """
    Checks if the compound ID 'C05359' is in the dictionary and removes it if it is.

    Parameters:
    compound_dict (dict): The dictionary to check and modify.

    Returns:
    dict: The modified dictionary with 'C05359' removed if it was present.
    """
    if 'C05359' in compound_dict:
        del compound_dict['C05359']
    return compound_dict


def convert_form_dict_to_elements(form_dict, strip_ionic=True):
    """
    Converts a dictionary of chemical formulas and their counts into a dictionary of elements and their counts.

    Parameters:
    form_dict (dict): A dictionary where keys are chemical formulas and values are their counts.
    strip_ionic (bool, optional): Flag to indicate whether to remove ionic states from the formulas. Default is True.

    Returns:
    dict: A dictionary where keys are element symbols and values are their counts in the formulas.
    """
    # Create a deep copy of the input dictionary to avoid modifying the original
    form_dict = copy.deepcopy(form_dict)

    # Initialize an empty dictionary to store the element counts
    elements = {}

    # Iterate over each formula and its count in the input dictionary
    for formula, count in form_dict.items():
        # Convert the formula to a dictionary of elements and their counts, multiply by the count, and add to the elements dictionary
        elements = add_dicts(elements, multiply_dict(convert_formula_to_dict(formula, strip_ionic=strip_ionic), count))

    # Return the dictionary of elements and their counts
    return elements


def get_eq(old_eq, reactants, products, c_data):
    """
    Generates a new chemical equation string by converting reactants and products from formulas to IDs.

    Parameters:
    old_eq (str): The original chemical equation string.
    reactants (dict): A dictionary where keys are chemical formulas and values are their counts for reactants.
    products (dict): A dictionary where keys are chemical formulas and values are their counts for products.
    c_data (DataFrame): A pandas DataFrame containing compound data with 'compound_id' and 'formula' columns.

    Returns:
    str: A string representation of the new chemical equation in the format 'reactants <=> products'.
    """
    lhs, rhs = eq_to_dict(old_eq)
    l_key = get_ids_to_formulas(lhs, c_data)
    r_key = get_ids_to_formulas(rhs, c_data)

    reactants = convert_formulas_to_ids(reactants, l_key)
    products = convert_formulas_to_ids(products, r_key)

    eq_left = [f"{v} {k}" for k, v in reactants.items()]
    eq_right = [f"{v} {k}" for k, v in products.items()]
    return " + ".join(eq_left) + " <=> " + " + ".join(eq_right)


def get_formulas_from_eq(eq, c_data, strip_ionic=True):
    """
    Converts a chemical equation string into dictionaries of reactants and products with their chemical formulas.

    Parameters:
    eq (str): A string representing a chemical equation, with reactants and products separated by '<=>'.
    c_data (DataFrame): A pandas DataFrame containing compound data with 'compound_id' and 'formula' columns.
    strip_ionic (bool, optional): Flag to indicate whether to remove ionic states from the formulas. Default is True.

    Returns:
    tuple: A tuple containing two dictionaries:
           - The first dictionary contains the converted reactants with chemical formulas as keys and their counts as values.
           - The second dictionary contains the converted products with chemical formulas as keys and their counts as values.
    """
    # Convert the Eq into the dicts
    reactants, products = eq_to_dict(eq)

    if strip_ionic:
        reactants = remove_electron_cid(reactants)
        products = remove_electron_cid(products)

    # Get the conversion of the ids to formulas
    react_id_form_key = get_ids_to_formulas(reactants, c_data)
    prod_id_form_key = get_ids_to_formulas(products, c_data)

    # Convert the reactants into formulas
    converted_reactants = convert_ids_to_formulas(reactants, react_id_form_key)
    converted_products = convert_ids_to_formulas(products, prod_id_form_key)

    return converted_reactants, converted_products


def get_elements_from_eq(eq, c_data, strip_ionic=True):
    """
    Converts a chemical equation string into dictionaries of reactants and products,
    and then converts these dictionaries into dictionaries of elements and their counts.

    Parameters:
    eq (str): A string representing a chemical equation, with reactants and products separated by '<=>'.
    c_data (DataFrame): A pandas DataFrame containing compound data with 'compound_id' and 'formula' columns.
    strip_ionic (bool, optional): Flag to indicate whether to remove ionic states from the formulas. Default is True.

    Returns:
    tuple: A tuple containing four dictionaries:
           - The first dictionary contains the converted reactants with chemical formulas as keys and their counts as values.
           - The second dictionary contains the converted products with chemical formulas as keys and their counts as values.
           - The third dictionary contains the elements and their counts in the reactants.
           - The fourth dictionary contains the elements and their counts in the products.
    """
    # Convert the Eq into the formula dicts
    converted_reactants, converted_products = get_formulas_from_eq(eq, c_data, strip_ionic=strip_ionic)

    # Convert the formulas into reactants
    react_ele = convert_form_dict_to_elements(converted_reactants, strip_ionic=strip_ionic)
    prod_ele = convert_form_dict_to_elements(converted_products, strip_ionic=strip_ionic)

    return converted_reactants, converted_products, react_ele, prod_ele


def get_missing_elements(react_ele, prod_ele):
    """
    Identifies the elements that are missing in the reactants and products.

    Parameters:
    react_ele (dict): A dictionary where keys are element symbols and values are their counts in the reactants.
    prod_ele (dict): A dictionary where keys are element symbols and values are their counts in the products.

    Returns:
    tuple: A tuple containing two lists:
           - The first list contains elements that are missing in the reactants.
           - The second list contains elements that are missing in the products.
    """
    missing_in_prod, missing_in_react = compare_dict_keys(react_ele, prod_ele)
    return list(missing_in_react), list(missing_in_prod)


def check_missing_elements(react_ele, prod_ele):
    """
    Checks if there are any missing elements in the reactants or products.

    Parameters:
    react_ele (dict): A dictionary where keys are element symbols and values are their counts in the reactants.
    prod_ele (dict): A dictionary where keys are element symbols and values are their counts in the products.

    Returns:
    bool: True if there are missing elements in either reactants or products, False otherwise.
    """
    missing_in_react, missing_in_prod = get_missing_elements(react_ele, prod_ele)
    if len(missing_in_react + missing_in_prod) > 0:
        return True
    else:
        return False


def check_full_missing_elements(eq, c_data):
    """
    Checks if there are any missing elements in the reactants or products of a chemical equation.

    Parameters:
    eq (str): A string representing a chemical equation, with reactants and products separated by '<=>'.
    c_data (DataFrame): A pandas DataFrame containing compound data with 'compound_id' and 'formula' columns.

    Returns:
    bool: True if there are missing elements in either reactants or products, False otherwise.
    """
    _, _, react_ele, prod_ele = get_elements_from_eq(eq, c_data)
    return check_missing_elements(react_ele, prod_ele)


def check_missing_formulas(eq, c_data):
    """
    Checks if there are any missing chemical formulas for the compound IDs in the given chemical equation.

    Parameters:
    eq (str): A string representing a chemical equation, with reactants and products separated by '<=>'.
    c_data (DataFrame): A pandas DataFrame containing compound data with 'compound_id' and 'formula' columns.

    Returns:
    bool: True if there are missing formulas for any compound IDs, False otherwise.
    """
    reactants, products = eq_to_dict(eq)
    ids = list(reactants.keys()) + list(products.keys())
    formulas = get_formulas_from_ids(ids, c_data)
    return len(formulas) != len(set(ids))


def full_check_eq_unbalanced(eq, c_data):
    """
    Checks if a chemical equation is unbalanced by verifying if all element counts are positive and non-zero,
    and if the reactants and products have matching element counts.

    Parameters:
    eq (str): A string representing a chemical equation, with reactants and products separated by '<=>'.
    c_data (DataFrame): A pandas DataFrame containing compound data with 'compound_id' and 'formula' columns.

    Returns:
    bool: True if the equation is unbalanced, False otherwise.
    """
    _, _, react_ele, prod_ele = get_elements_from_eq(eq, c_data)
    return check_eq_unbalanced(react_ele, prod_ele)


def check_eq_unbalanced(react_ele, prod_ele):
    """
    Checks if a chemical equation is unbalanced by verifying if all element counts are positive and non-zero,
    and if the reactants and products have matching element counts.

    Parameters:
    react_ele (dict): A dictionary where keys are element symbols and values are their counts in the reactants.
    prod_ele (dict): A dictionary where keys are element symbols and values are their counts in the products.

    Returns:
    bool: True if the equation is unbalanced, False otherwise.
    """
    # Check if all values in the dictionaries are positive and non-zero
    all_positive_react = all(value > 0 for value in react_ele.values())
    all_positive_prod = all(value > 0 for value in prod_ele.values())

    # Check if the dictionaries are unbalanced
    diff_ele_react, diff_ele_prod = compare_dict_values(react_ele, prod_ele)
    unbalanced = len(diff_ele_react) + len(diff_ele_prod) > 0

    return not all_positive_react or not all_positive_prod or unbalanced


def sort_dict_by_keys(input_dict):
    """
    Returns a dictionary sorted by its keys.

    Parameters:
    input_dict (dict): The dictionary to be sorted.

    Returns:
    dict: A new dictionary sorted by its keys.
    """
    input_dict = copy.deepcopy(input_dict)
    return dict(sorted(input_dict.items()))


def ordered_reaction_series(reaction_series):
    """
    Sorts a pandas Series of reaction equations such that reaction directionality is not retained.

    Parameters:
    reaction_series (pd.Series): A Series containing reaction equations.

    Returns:
    pd.Series: A Series with sorted reaction equations.
    """
    return (reaction_series.str.split(' <=> ', expand=True)
            .map(lambda x: [' + '.join(sorted(x.split(' + ')))])
            .sum(axis=1).apply(lambda x: ' <=> '.join(sorted(x))))


def generate_reaction_dupemap(df, prefix="T"):
    """
    Transforms the 'reaction' column of a DataFrame to identify duplicate reactions and returns a map of oldID to newID with a user-specified prefix.

    Parameters:
    df (pd.DataFrame): The DataFrame containing reaction data with an 'id' column and a 'reaction' column.
    prefix (str): The prefix to use for the new IDs. Defaults to "T".

    Returns:
    pd.Series: A Series mapping old IDs to new IDs with the specified prefix.
    """
    # Group the reactions by sorting and standardizing them
    equation_groups = ordered_reaction_series(df.set_index('id')['reaction']).to_frame(name='eq_grp')

    # Identify duplicated groups and assign a group number
    duped_groups = equation_groups.query('eq_grp.duplicated(keep=False)').assign(
        grp_num=lambda x: x.eq_grp.factorize()[0])

    # Create a map of old IDs to new IDs with the specified prefix
    dupemap = prefix + duped_groups['grp_num'].astype(str).str.zfill(5)

    # Save the map to a CSV file
    dupemap.sort_index().to_csv(f'../data/{prefix}_reaction_dupemap.csv', encoding='utf-8')

    return dupemap


def standardise_eq(eq):
    """
    Standardise the equation by sorting the reactants and products.

    Parameters:
    eq (str): The equation to be standardised.

    Returns:
    str: The standardised equation.
    """
    reactants, products = eq_to_dict(eq)
    return dicts_to_eq(sort_dict_by_keys(reactants), sort_dict_by_keys(products))


def contains_var_list(reactants, products, var_list=None):
    """
    Checks if any of the values in the reactants or products dictionaries contain any of the variables in var_list.
    Returns a list of the variables that are present.

    Parameters:
    reactants (dict): A dictionary where keys are compound identifiers and values are their counts for reactants.
    products (dict): A dictionary where keys are compound identifiers and values are their counts for products.
    var_list (list, optional): A list of variable names to check for in the values. Default is ['n', 'm', 'x'].

    Returns:
    list: A list of variables that are present in any value.
    """
    if var_list is None:
        var_list = ['n', 'm', 'x']
    # Convert all the dict values to strings
    reactants = {k: str(v) for k, v in reactants.items()}
    products = {k: str(v) for k, v in products.items()}
    present_vars = []
    for value in reactants.values():
        for var in var_list:
            if var in value and var not in present_vars:
                present_vars.append(var)
    for value in products.values():
        for var in var_list:
            if var in value and var not in present_vars:
                present_vars.append(var)
    return present_vars


def check_contains_var_list(eq, data_c, var_list=None):
    """
    Checks if any of the values in the reactants or products dictionaries contain any of the variables in var_list.

    Parameters:
    eq (str): A string representing a chemical equation, with reactants and products separated by '<=>'.
    data_c (DataFrame): A pandas DataFrame containing compound data with 'compound_id' and 'formula' columns.
    var_list (list, optional): A list of variable names to check for in the values. Default is ['n', 'm', 'x'].

    Returns:
    bool: True if any of the values contain a variable from var_list, False otherwise.
    """
    reactants, products = get_formulas_from_eq(eq, data_c)
    return len(contains_var_list(reactants, products, var_list)) > 0


def solve_for(elements, var='n'):
    """
    Solves for the smallest value of the variable that makes each element expression greater than 0.

    Parameters:
    elements (list): A list of strings representing element expressions.
    var (str): The variable to solve for. Default is 'n'.

    Returns:
    int: The smallest value of the variable that makes each element expression greater than 0.
    """

    min_var = 1
    for element in elements:
        # Replace var with a symbolic variable
        expr = element  # .replace(var, var)
        # Solve for the smallest var that makes the expression greater than 0
        var_value = eval(expr.replace(var, '0'))
        if var_value <= 0:
            min_var = max(min_var, -var_value + 1)
    return min_var


def find_min_integers(expr_str):
    """
    Finds the smallest integer values of the variables that make the expression positive.

    Note this finds a valid solution, but it may not be the smallest possible value.

    Parameters:
    expr_str (str): The input expression as a string.

    Returns:
    dict: A dictionary with variable names as keys and their smallest positive integer values as values.
    """
    # Parse the expression string
    expr = sp.sympify(expr_str)

    # Extract the variables from the expression
    variables = expr.free_symbols

    # Initialize a dictionary to store the results
    result = {}

    # Solve iteratively for each variable
    for var in variables:
        # Set all other variables to 1 for simplicity
        assumptions = {v: 1 for v in variables if v != var}
        # Solve for the current variable
        sol = sp.solve(expr.subs(assumptions) - 1, var)
        # Find the smallest integer value that satisfies the condition
        min_val = max(1, sp.ceiling(sol[0])) if sol else 1
        result[var] = min_val

    # Return the results as a dictionary where the keys are strings
    return {str(k): v for k, v in result.items()}


def delete_pm_keys(dictionary):
    """
    Deletes any keys in the dictionary that are equal to '+' or '-'.

    Parameters:
    dictionary (dict): The dictionary from which to delete the keys.

    Returns:
    dict: The dictionary with the specified keys removed.
    """
    keys_to_remove = ['+', '-']
    for key in keys_to_remove:
        if key in dictionary:
            del dictionary[key]
    return dictionary


def fix_multiply_tar(expression, target_letter):
    """
    Replaces occurrences of a number followed by a target letter with the number followed by '*' and the target letter.

    Parameters:
    expression (str): The input string containing the expression.
    target_letter (str): The target letter to match after the number.

    Returns:
    str: The modified expression with replacements.
    """
    pattern = rf'(\d+){target_letter}'
    replacement = rf'\1*{target_letter}'
    return re.sub(pattern, replacement, expression)


def fix_multiply_tar_all(expression, target_letters=None):
    """
    Replaces occurrences of a number followed by any target letter in the target_letters list with
    the number followed by '*' and the target letter.

    Parameters:
    expression (str): The input string containing the expression.
    target_letters (list, optional): A list of target letters to match after the number. Defaults to ["n", "m", "x"].

    Returns:
    str: The modified expression with replacements.
    """
    if target_letters is None:
        target_letters = ["n", "m", "x"]
    for target_letter in target_letters:
        expression = fix_multiply_tar(expression, target_letter)
    return expression


def check_vars_eq_balanced(reactants, products):
    """
    Checks if a chemical equation with variables is balanced by solving for the variable values and evaluating the equation.

    Parameters:
    reactants (dict): A dictionary where keys are compound identifiers and values are their counts for reactants.
    products (dict): A dictionary where keys are compound identifiers and values are their counts for products.

    Returns:
    bool: True if the equation is balanced, False otherwise.
    """
    # Convert all the dict values to strings
    reactants = {k: str(v) for k, v in reactants.items()}
    products = {k: str(v) for k, v in products.items()}

    # Get the values in the reactants and products
    reactants_values = list(reactants.values())
    products_values = list(products.values())
    full_list = reactants_values + products_values

    # Solve for n
    n_val = solve_for(full_list)
    print(f"n = {n_val}", flush=True)

    # Substitute the n value into the reactants and products
    reactants = {k: fix_multiply_tar(v, 'n').replace('n', str(n_val)) for k, v in reactants.items()}
    products = {k: fix_multiply_tar(v, 'n').replace('n', str(n_val)) for k, v in products.items()}
    print(reactants, flush=True)
    print(products, flush=True)

    # eval the values in the reactants and products
    reactants = {k: eval(v) for k, v in reactants.items()}
    products = {k: eval(v) for k, v in products.items()}
    print(reactants, flush=True)
    print(products, flush=True)

    return check_eq_unbalanced(reactants, products)


def replace_substrings(s: str, replacements: dict):
    """
    Within a string, replaces all instances of each key in a dictionary with its corresponding value.
    Keys and values must all be strings. Values may be empty strings.

    Parameters:
    s (str): The string

    replacements (dict): A dictionary with format {"string_to_be_replaced" : "replacement_string"}

    Returns:
    s_out: The string with all replacements made.
    """
    s_out = s
    for k, v in replacements.items():
        s_out = s_out.replace(k, v)
    return s_out


def inject_compounds(eq_line, missing_r, missing_p, missing_dict):
    """
    Injects missing compounds into a reaction equation.

    Parameters:
    eq_line (str): The original reaction equation line.
    missing_r (list): A list of missing reactant compound IDs.
    missing_p (list): A list of missing product compound IDs.
    missing_dict (dict): A dictionary mapping missing compound names to their IDs.

    Returns:
    str: The updated reaction equation with missing compounds injected.
    """
    eq_left, eq_right = map(str.strip, eq_line.split("<=>"))
    eq_left += ''.join(f" + {missing_dict[item]}" for item in missing_r)
    eq_right += ''.join(f" + {missing_dict[item]}" for item in missing_p)
    return f"{eq_left} <=> {eq_right}"


def compare_and_delete_keys(dict1, dict2):
    """
    Compares two dictionaries, deletes keys that are the same in both dictionaries,
    and returns the two new dictionaries along with the deleted items.

    Parameters:
    dict1 (dict): The first dictionary.
    dict2 (dict): The second dictionary.

    Returns:
    tuple: A tuple containing four dictionaries:
           - The first dictionary with the common keys removed from dict1.
           - The second dictionary with the common keys removed from dict2.
           - A dictionary of the deleted items from dict1.
           - A dictionary of the deleted items from dict2.
    """
    dict1 = copy.deepcopy(dict1)
    dict2 = copy.deepcopy(dict2)
    common_keys = set(dict1.keys()) & set(dict2.keys())
    dict1_del = {key: dict1[key] for key in dict1 if key in common_keys}
    dict2_del = {key: dict2[key] for key in dict2 if key in common_keys}
    for key in common_keys:
        del dict1[key]
        del dict2[key]
    return dict1, dict2, dict1_del, dict2_del


def rebalance_eq_core(eq, data_c):
    """
    Rebalances a chemical equation by converting reactants and products from formulas to IDs,
    checking for duplicate keys, and balancing the stoichiometry.

    Parameters:
    eq (str): The original chemical equation string.
    data_c (DataFrame): A pandas DataFrame containing compound data with 'compound_id' and 'formula' columns.

    Returns:
    str: The rebalanced and standardized chemical equation.
    """
    # Get the reactants and products from the eq
    reactants, products = get_formulas_from_eq(eq, data_c)

    f_repeated = False
    reactants_del = {}
    products_del = {}
    # Check if any of the keys are duplicate in the reactants and products
    if len(set(reactants.keys()) & set(products.keys())) > 0:
        f_repeated = True
        # Get the unique reactants and products
        reactants, products, reactants_del, products_del = compare_and_delete_keys(reactants, products)

    # Balance the equation
    reactants, products = balance_stoichiometry(set(reactants.keys()),
                                                set(products.keys()),
                                                underdetermined=None)
    # Convert the formulas back to eq form
    reactants = dict(reactants)
    products = dict(products)

    # Check if the reactants and products were the same
    if f_repeated:
        # In case of repeated compounds, add the deleted compounds back
        reactants = add_dicts(reactants, reactants_del)
        products = add_dicts(products, products_del)

    # Convert the dict back into eq form and standardise
    return standardise_eq(get_eq(eq, dict(reactants), dict(products), data_c))


def rebalance_eq(eq, data_c):
    """
    Attempts to rebalance a chemical equation by calling the rebalance_eq_core function.
    If rebalancing fails due to a ValueError, it prints an error message and returns False.

    Parameters:
    eq (str): The original chemical equation string.
    data_c (DataFrame): A pandas DataFrame containing compound data with 'compound_id' and 'formula' columns.

    Returns:
    str or bool: The rebalanced chemical equation string if successful, otherwise False.
    """
    try:
        eq = rebalance_eq_core(eq, data_c)
        return eq
    except ValueError as e:
        print(f"Could not find stoichiometry on first attempt: {e}", flush=True)
        return False


def fix_imbalance_core(eq_line, diff_ele_react, diff_ele_prod, inject):
    """
    Fixes the imbalance in a reaction equation by injecting a specified compound.

    Parameters:
    eq_line (str): The original reaction equation line.
    diff_ele_react (dict): A dictionary of element differences on the reactant side.
    diff_ele_prod (dict): A dictionary of element differences on the product side.
    inject (str): The compound ID to inject to balance the equation.

    Returns:
    str: The updated reaction equation with the injected compound.
    """
    # Get the reactant and product sides
    eq_left, eq_right = map(str.strip, eq_line.split("<=>"))
    # Find which side has the lowest
    if sum(diff_ele_react.values()) < sum(diff_ele_prod.values()):
        eq_left += f" + {inject}"
    else:
        eq_right += f" + {inject}"
    # Update eq_line with the new equation
    return standardise_eq(f"{eq_left} <=> {eq_right}")


def plot_eq_line(eq, data_c, render_dir='render', size=(600, 600)):
    """
    Plots the reactants and products of a chemical equation.

    Parameters:
    eq (str): The chemical equation to plot.
    data_c (pd.DataFrame): A pandas DataFrame containing compound data with 'compound_id' and 'smiles' columns.
    render_dir (str, optional): The directory where the reaction images will be saved. Default is 'render'.
    size (tuple, optional): The size of the images to be rendered. Default is (600, 600).

    Returns:
    tuple: A tuple containing two lists:
           - The first list contains the keys of the reactants.
           - The second list contains the keys of the products.
    """
    print(f'Input eq: {eq}', flush=True)
    eq = standardise_eq(eq)

    # Get the reactants and products from the eq
    reactants, products = eq_to_dict(eq)

    # Get the keys of reactants as a list
    reactants_keys = list(reactants.keys())
    reactant_values = list(reactants.values())
    products_keys = list(products.keys())
    products_values = list(products.values())

    # Get the smiles for the reactants and products
    reactants_smi = data_c.loc[data_c['compound_id'].isin(reactants_keys), 'smiles'].tolist()
    products_smi = data_c.loc[data_c['compound_id'].isin(products_keys), 'smiles'].tolist()

    # Make a directory
    os.makedirs(render_dir, exist_ok=True)

    # Plot the reactants and products
    print('Reactants:', flush=True)
    for i, smi in enumerate(reactants_smi):
        print(f'{i}: {reactant_values[i]} x {reactants_keys[i]}, {smi}', flush=True)
        mol = Chem.MolFromSmiles(smi)
        Draw.MolToFile(mol, f'{render_dir}/reactant_{i}_{reactants_keys[i]}.png', size=size)

    print('Products:', flush=True)
    for i, smi in enumerate(products_smi):
        print(f'{i}: {products_values[i]} x {products_keys[i]}, {smi}', flush=True)
        mol = Chem.MolFromSmiles(smi)
        Draw.MolToFile(mol, f'{render_dir}/product_{i}_{products_keys[i]}.png', size=size)

    return reactants_keys, products_keys


def plot_reaction_id(id, data_r, data_c, render_dir='render', size=(600, 600)):
    """
    Plots the reaction corresponding to the given reaction ID.

    Parameters:
    id (int): The reaction ID to plot.
    data_r (pd.DataFrame): A pandas DataFrame containing reaction data with 'id' and 'reaction' columns.
    data_c (pd.DataFrame): A pandas DataFrame containing compound data with 'compound_id' and 'smiles' columns.
    render_dir (str, optional): The directory where the reaction images will be saved. Default is 'render'.
    size (tuple, optional): The size of the images to be rendered. Default is (600, 600).

    Returns:
    tuple: A tuple containing two lists:
           - The first list contains the keys of the reactants.
           - The second list contains the keys of the products.
    """
    # Get the reaction from the id
    eq = data_r.loc[data_r['id'] == id, 'reaction'].values[0]
    # Plot the eq
    return plot_eq_line(eq, data_c, render_dir=render_dir, size=size)


def to_smarts_rxn_line(eq, data_c):
    """
    Converts a chemical equation into a SMARTS reaction line.

    Parameters:
    -----------
    eq : str
        A string representing a chemical equation, with reactants and products separated by '<=>'.
    data_c : pandas.DataFrame
        A DataFrame containing compound data with 'compound_id' and 'smiles' columns.

    Returns:
    --------
    str
        A SMARTS reaction line in the format 'reactants>>products', where reactants and products are represented
        as concatenated SMARTS strings.

    Notes:
    ------
    - The function retrieves SMILES strings for reactants and products from the `data_c` DataFrame.
    - Reactants and products are combined into SMARTS strings using their stoichiometric coefficients.
    """
    reactants, products = eq_to_dict(eq)

    # Get SMILES for reactants and products
    reactants_smiles = {k: data_c.loc[data_c['compound_id'] == k, 'smiles'].values[0] for k in reactants}
    products_smiles = {k: data_c.loc[data_c['compound_id'] == k, 'smiles'].values[0] for k in products}

    # Combine SMILES into SMARTS strings
    reactants_smarts = ".".join(f"{reactants[k]}{v}" for k, v in reactants_smiles.items())
    products_smarts = ".".join(f"{products[k]}{v}" for k, v in products_smiles.items())

    return f"{reactants_smarts}>>{products_smarts}"
