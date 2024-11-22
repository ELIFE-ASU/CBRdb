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
    # formula_new = strip_ionic_states(formula_new)

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

    Parameters:
    dict1 (dict): The first dictionary.
    dict2 (dict): The second dictionary.

    Returns:
    tuple: A tuple containing two dictionaries:
           - The first dictionary contains key-value pairs from dict1 where the values differ from dict2.
           - The second dictionary contains key-value pairs from dict2 where the values differ from dict1.
    """
    diff_in_dict1 = {key: dict1.get(key) for key in set(dict1) | set(dict2) if dict1.get(key) != dict2.get(key)}
    diff_in_dict2 = {key: dict2.get(key) for key in set(dict1) | set(dict2) if dict1.get(key) != dict2.get(key)}
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
    return map(side_to_dict, eq.split('<=>'))


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
    ids = list(set(sorted(compound_dict.keys())))
    formulas = get_formulas_from_ids(ids, c_data)
    return {id: strip_plus_x(formula) for id, formula in zip(ids, formulas)}


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


def convert_form_dict_to_elements(form_dict):
    """
    Converts a dictionary of chemical formulas and their counts into a dictionary of elements and their counts.

    Parameters:
    form_dict (dict): A dictionary where keys are chemical formulas and values are their counts.

    Returns:
    dict: A dictionary where keys are element symbols and values are their counts in the formulas.
    """
    elements = {}
    for formula, count in form_dict.items():
        elements = add_dicts(elements, multiply_dict(convert_formula_to_dict(formula), count))
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


def get_formulas_from_eq(eq, c_data):
    """
    Converts a chemical equation string into dictionaries of reactants and products with their chemical formulas.

    Parameters:
    eq (str): A string representing a chemical equation, with reactants and products separated by '<=>'.
    c_data (DataFrame): A pandas DataFrame containing compound data with 'compound_id' and 'formula' columns.

    Returns:
    tuple: A tuple containing two dictionaries:
           - The first dictionary contains the converted reactants with chemical formulas as keys and their counts as values.
           - The second dictionary contains the converted products with chemical formulas as keys and their counts as values.
    """
    # Convert the Eq into the dicts
    reactants, products = eq_to_dict(eq)

    # Get the conversion of the ids to formulas
    react_id_form_key = get_ids_to_formulas(reactants, c_data)
    prod_id_form_key = get_ids_to_formulas(products, c_data)

    # Convert the reactants into formulas
    converted_reactants = convert_ids_to_formulas(reactants, react_id_form_key)
    converted_products = convert_ids_to_formulas(products, prod_id_form_key)

    return converted_reactants, converted_products


def get_elements_from_eq(eq, c_data):
    """
    Converts a chemical equation string into dictionaries of reactants and products,
    and then converts these dictionaries into dictionaries of elements and their counts.

    Parameters:
    eq (str): A string representing a chemical equation, with reactants and products separated by '<=>'.
    c_data (DataFrame): A pandas DataFrame containing compound data with 'compound_id' and 'formula' columns.

    Returns:
    tuple: A tuple containing four dictionaries:
           - The first dictionary contains the converted reactants with chemical formulas as keys and their counts as values.
           - The second dictionary contains the converted products with chemical formulas as keys and their counts as values.
           - The third dictionary contains the elements and their counts in the reactants.
           - The fourth dictionary contains the elements and their counts in the products.
    """
    # Convert the Eq into the formula dicts
    converted_reactants, converted_products = get_formulas_from_eq(eq, c_data)

    # Convert the formulas into reactants
    react_ele = convert_form_dict_to_elements(converted_reactants)
    prod_ele = convert_form_dict_to_elements(converted_products)

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
    return len(formulas) != len(ids)


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
    return dict(sorted(input_dict.items()))


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


def contains_n_m_x(reactants, products):
    """
    Checks if any of the reactant or product values contain the strings 'n', 'm', or 'x'.

    Parameters:
    reactants (dict): A dictionary of reactants with chemical formulas as keys and their counts as values.
    products (dict): A dictionary of products with chemical formulas as keys and their counts as values.

    Returns:
    bool: True if any reactant or product value contains 'n', 'm', or 'x', False otherwise.
    """
    for value in reactants.values():
        if any(char in value for char in ['n', 'm', 'x']):
            return True
    for value in products.values():
        if any(char in value for char in ['n', 'm', 'x']):
            return True
    return False


def solve_for_n(elements):
    """
    Solves for the smallest integer n that makes all elements in the list greater than 0.

    Parameters:
    elements (list): A list of strings representing expressions involving 'n'.

    Returns:
    int: The smallest integer n that makes all elements greater than 0.
    """
    min_n = 0
    for element in elements:
        # Replace 'n' with a symbolic variable
        expr = element.replace('n', 'n')
        # Solve for the smallest n that makes the expression greater than 0
        n_value = eval(expr.replace('n', '0'))
        if n_value <= 0:
            min_n = max(min_n, -n_value + 1)
    return min_n
