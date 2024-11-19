import re


def split_by_letters(input_string):
    return re.findall(r'[A-Z][a-z]*\d*', input_string)


def convert_formula_to_dict(input_string):
    parts = split_by_letters(input_string)
    result = {}
    for part in parts:
        element = re.match(r'[A-Za-z]+', part).group()
        count = int(re.search(r'\d+', part).group()) if re.search(r'\d+', part) else 1
        result[element] = count
    return result


def multiply_dict(input_dict, multiplier):
    return {key: value * multiplier for key, value in input_dict.items()}


def add_dicts(dict1, dict2):
    return {key: dict1.get(key, 0) + dict2.get(key, 0) for key in set(dict1) | set(dict2)}


def subtract_dicts(dict1, dict2):
    return {key: dict1.get(key, 0) - dict2.get(key, 0) for key in set(dict1) | set(dict2)}


def compare_dict_keys(dict1, dict2):
    missing_in_dict2 = set(dict1) - set(dict2)
    missing_in_dict1 = set(dict2) - set(dict1)
    return missing_in_dict2, missing_in_dict1


def compare_dicts(dict1, dict2):
    return dict1.keys() == dict2.keys() and all(dict1[key] == dict2[key] for key in dict1)


def compare_dict_values(dict1, dict2):
    diff_in_dict1 = {key: dict1.get(key) for key in set(dict1) | set(dict2) if dict1.get(key) != dict2.get(key)}
    diff_in_dict2 = {key: dict2.get(key) for key in set(dict1) | set(dict2) if dict1.get(key) != dict2.get(key)}
    return diff_in_dict1, diff_in_dict2


def get_formulas_from_ids(ids, data):
    return data.loc[data["compound_id"].isin(ids), "formula"].tolist()


def side_to_dict(side):
    result = {}
    for component in map(str.strip, side.split('+')):
        match = re.match(r'(\d*)\s*(C\d+)', component)
        if match:
            count = int(match.group(1)) if match.group(1) else 1
            molecule = match.group(2)
            result[molecule] = count
    return result


def eq_to_dict(eq):
    return map(side_to_dict, eq.split('<=>'))


def dict_to_side(d):
    return ' + '.join(f"{v} {k}" for k, v in d.items())


def dicts_to_eq(reactants, products):
    return f"{dict_to_side(reactants)} <=> {dict_to_side(products)}"


def strip_plus_x(input_string):
    return re.sub(r'\+\d+', '', input_string).replace('+', '') if '+' in input_string else input_string


def get_ids_to_formulas(compound_dict, data):
    ids = list(set(sorted(compound_dict.keys())))
    formulas = get_formulas_from_ids(ids, data)
    return {id: strip_plus_x(formula) for id, formula in zip(ids, formulas)}


def convert_ids_to_formulas(in_dict, react_id_form):
    return {react_id_form[id]: count for id, count in in_dict.items()}


def convert_formulas_to_ids(formulas_dict, react_id_form):
    inverse_react_id_form = {v: k for k, v in react_id_form.items()}
    return {inverse_react_id_form[formula]: count for formula, count in formulas_dict.items()}


def convert_form_dict_to_elements(form_dict):
    elements = {}
    for formula, count in form_dict.items():
        elements = add_dicts(elements, multiply_dict(convert_formula_to_dict(formula), count))
    return elements


def get_eq(old_eq, reactants, products, data):
    lhs, rhs = eq_to_dict(old_eq)
    l_key = get_ids_to_formulas(lhs, data)
    r_key = get_ids_to_formulas(rhs, data)

    reactants = convert_formulas_to_ids(reactants, l_key)
    products = convert_formulas_to_ids(products, r_key)

    eq_left = [f"{v} {k}" for k, v in reactants.items()]
    eq_right = [f"{v} {k}" for k, v in products.items()]
    return " + ".join(eq_left) + " <=> " + " + ".join(eq_right)


def get_elements_from_eq(eq, data):
    # Convert the Eq in to the dicts
    reactants, products = eq_to_dict(eq)
    # Get the conversion of the ids to formulas
    react_id_form_key = get_ids_to_formulas(reactants, data)
    prod_id_form_key = get_ids_to_formulas(products, data)

    # Convert the reactants into formulas
    converted_reactants = convert_ids_to_formulas(reactants, react_id_form_key)
    converted_products = convert_ids_to_formulas(products, prod_id_form_key)

    # Convert the formulas into reactants
    react_ele = convert_form_dict_to_elements(converted_reactants)
    prod_ele = convert_form_dict_to_elements(converted_products)

    return converted_reactants, converted_products, react_ele, prod_ele


def get_missing_elements(react_ele, prod_ele):
    missing_in_prod, missing_in_react = compare_dict_keys(react_ele, prod_ele)
    return list(missing_in_react), list(missing_in_prod)


def check_missing_elements(react_ele, prod_ele):
    missing_in_react, missing_in_prod = get_missing_elements(react_ele, prod_ele)
    if len(missing_in_react + missing_in_prod) > 0:
        return True
    else:
        return False


def check_missing_formulas(eq, data):
    reactants, products = eq_to_dict(eq)
    ids = list(reactants.keys()) + list(products.keys())
    formulas = get_formulas_from_ids(ids, data)
    return len(formulas) != len(ids)


def check_eq_unbalanced(react_ele, prod_ele):
    # Check if all values in the dictionaries are positive and non-zero
    all_positive_react = all(value > 0 for value in react_ele.values())
    all_positive_prod = all(value > 0 for value in prod_ele.values())

    # Check if the dictionaries are unbalanced
    diff_ele_react, diff_ele_prod = compare_dict_values(react_ele, prod_ele)
    unbalanced = len(diff_ele_react) + len(diff_ele_prod) > 0

    return not all_positive_react or not all_positive_prod or unbalanced
