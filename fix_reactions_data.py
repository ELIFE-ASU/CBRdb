import os
import re
import time

import pandas as pd
import requests
from chempy import balance_stoichiometry
from requests import Session
from requests.adapters import HTTPAdapter
from urllib3.util import Retry


def check_known_ids(id):
    if id == "C00138":  # Reduced ferredoxin
        return "Fe2S2"
    elif id == "C00139":  # Oxidized ferredoxin
        return "Fe2S2"
    elif id == "C05359":  # electron
        return " "
    else:
        return False


def prepare_session():
    # Make the session
    s = Session()
    # Add retries
    retries = Retry(
        total=5,
        backoff_factor=0.1,
        status_forcelist=[502, 503, 504],
        allowed_methods={'POST'},
    )
    # Mount the session
    s.mount('https://', HTTPAdapter(max_retries=retries))
    return s


def get_formula_from_web(id, kegg_website="https://rest.kegg.jp/get/", request_sleep=0.2):
    # Prepare the session
    s = prepare_session()
    # Get the data
    try:
        response = s.get(f"{kegg_website}{id}", timeout=10.0)
        # Limit the number of requests
        time.sleep(request_sleep)
    except requests.exceptions.RequestException as e:
        # Some error in the connection
        print(f"Error in ID {id}, connection exception {e}")
        return None
    # Check if the response is ok
    if response.ok:
        # Strip the last section of the file
        res = response.text.split("> <ENTRY>")[0]
        data = res.split("\n")
        # Get the line which contains the formula
        data = [line for line in data if "FORMULA" in line]
        if len(data) == 0:
            return None
        return data[0].split("FORMULA")[1].strip()
    else:
        # Some error in the response
        print(f"Error in ID {id}, response {response.status_code}")
        return None


def file_list(mypath=None):
    mypath = mypath or os.getcwd()  # If no path is provided, use the current working directory
    return [f for f in os.listdir(mypath) if
            os.path.isfile(os.path.join(mypath, f))]  # Return a list of all files in the directory


def filter_files(file_paths, substring):
    return list(filter(lambda file_path: substring in os.path.basename(file_path), file_paths))


def sub_file_list(mypath, substring):
    paths = file_list(mypath)
    return filter_files(paths, substring)


def file_list_all(mypath=None):
    mypath = mypath or os.getcwd()  # If no path is provided, use the current working directory
    files = []
    # os.walk generates the file names in a directory tree by walking the tree either top-down or bottom-up
    for dirpath, dirnames, filenames in os.walk(mypath):
        for filename in filenames:
            # os.path.join joins one or more path components intelligently
            files.append(os.path.join(dirpath, filename))
    return files


def select_numbers(arr):
    return [x for x in arr if x.isdigit()]


def split_by_letters(input_string):
    return re.findall(r'[A-Z][a-z]*\d*', input_string)


def convert_formula_to_dict(input_string):
    parts = split_by_letters(input_string)
    result = {}
    for part in parts:
        element = re.match(r'[A-Za-z]+', part).group()
        count = re.search(r'\d+', part)
        count = int(count.group()) if count else 1
        result[element] = count
    return result


def multiply_dict(input_dict, multiplier):
    return {key: value * multiplier for key, value in input_dict.items()}


def add_dicts(dict1, dict2):
    result = {}
    for key in set(dict1) | set(dict2):
        result[key] = dict1.get(key, 0) + dict2.get(key, 0)
    return result


def subtract_dicts(dict1, dict2):
    result = {}
    for key in set(dict1) | set(dict2):
        result[key] = dict1.get(key, 0) - dict2.get(key, 0)
    return result


def compare_dict_keys(dict1, dict2):
    """
    # Example usage
    dict1 = {'a': 1, 'b': 2, 'c': 3}
    dict2 = {'b': 3, 'c': 4, 'd': 5}
    missing_in_dict2, missing_in_dict1 = compare_dict_keys(dict1, dict2)
    print(f"Keys missing in dict2: {missing_in_dict2}")  # Output: {'a'}
    print(f"Keys missing in dict1: {missing_in_dict1}")  # Output: {'d'}
    Args:
        dict1:
        dict2:

    Returns:

    """
    keys1 = set(dict1.keys())
    keys2 = set(dict2.keys())
    missing_in_dict2 = keys1 - keys2
    missing_in_dict1 = keys2 - keys1
    return missing_in_dict2, missing_in_dict1


def compare_dicts(dict1, dict2):
    """
    # Example usage
    dict1 = {'a': 1, 'b': 2, 'c': 3}
    dict2 = {'a': 1, 'b': 2, 'c': 3}
    dict3 = {'a': 1, 'b': 2, 'c': 4}
    print(compare_dicts(dict1, dict2))  # Output: True
    print(compare_dicts(dict1, dict3))  # Output: False
    Args:
        dict1:
        dict2:

    Returns:

    """
    # Check if both dictionaries have the same keys
    if dict1.keys() != dict2.keys():
        return False

    # Check if the values for each key are the same
    for key in dict1:
        if dict1[key] != dict2[key]:
            return False

    return True


def compare_dict_values(dict1, dict2):
    """
    # Example usage
    dict1 = {'a': 1, 'b': 2, 'c': 3}
    dict2 = {'a': 1, 'b': 3, 'c': 4, 'd': 5}
    diff_in_dict1, diff_in_dict2 = compare_dict_values(dict1, dict2)
    print(f"Differences in dict1: {diff_in_dict1}")
    print(f"Differences in dict2: {diff_in_dict2}")
    Args:
        dict1:
        dict2:

    Returns:

    """
    diff_in_dict1 = {}
    diff_in_dict2 = {}

    for key in set(dict1) | set(dict2):
        value1 = dict1.get(key)
        value2 = dict2.get(key)
        if value1 != value2:
            diff_in_dict1[key] = value1
            diff_in_dict2[key] = value2

    return diff_in_dict1, diff_in_dict2


def sum_formulas(formulas):
    total = {}
    for formula in formulas:
        formula_dict = convert_formula_to_dict(formula)
        total = add_dicts(total, formula_dict)
    return total


def get_formulas_from_ids(ids, file_path='Data/kegg_data_C.csv.zip'):
    data = pd.read_csv(file_path)
    return data.loc[data["compound_id"].isin(ids), "formula"].tolist()


def side_to_dict(side):
    components = side.split('+')
    result = {}
    for component in components:
        component = component.strip()
        match = re.match(r'(\d*)\s*(C\d+)', component)
        if match:
            count = int(match.group(1)) if match.group(1) else 1
            molecule = match.group(2)
            result[molecule] = count
    return result


def eq_to_dict(eq):
    reactants_side, products_side = eq.split('<=>')
    reactants = side_to_dict(reactants_side)
    products = side_to_dict(products_side)
    return reactants, products


def dict_to_side(d):
    return ' + '.join([f"{v} {k}" for k, v in d.items()])


def dicts_to_eq(reactants, products):
    reactants_str = dict_to_side(reactants)
    products_str = dict_to_side(products)
    return f"{reactants_str} <=> {products_str}"


def strip_plus_x(input_string):
    # Check if the input string has a plus and attempt to remove it
    if "+" in input_string:
        return re.sub(r'\+\d+', '', input_string).replace('+', '')
    else:
        return input_string


def get_ids_to_formulas(compound_dict, web=False, file_path='Data/kegg_data_C.csv.zip'):
    ids = list(compound_dict.keys())
    # sort the ids
    ids.sort()
    if web:
        formulas = [get_formula_from_web(id) for id in ids]
    else:
        formulas = get_formulas_from_ids(ids, file_path)

    # Remake the dict and strip the plus x from the formulas
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


def get_eq(old_eq, reactants, products, web=False, file_path='Data/kegg_data_C.csv.zip'):
    # Convert the Eq in to the dicts
    lhs, rhs = eq_to_dict(old_eq)
    # Get the conversion of the ids to formulas
    l_key = get_ids_to_formulas(lhs, web=web, file_path=file_path)
    r_key = get_ids_to_formulas(rhs, web=web, file_path=file_path)

    # Convert the dict back into eq form
    reactants = convert_formulas_to_ids(reactants, l_key)
    products = convert_formulas_to_ids(products, r_key)

    react_keys = list(reactants.keys())
    react_vals = list(reactants.values())
    prod_keys = list(products.keys())
    prod_vals = list(products.values())

    eq_left = [f"{react_vals[i]} {react_keys[i]}" for i in range(len(react_keys))]
    eq_right = [f"{prod_vals[i]} {prod_keys[i]}" for i in range(len(prod_keys))]
    return " + ".join(eq_left) + " <=> " + " + ".join(eq_right)


def get_elements_from_eq(eq, verbose=False, web=False, file_path='Data/kegg_data_C.csv.zip'):
    # Convert the Eq in to the dicts
    reactants, products = eq_to_dict(eq)
    if verbose:
        print("Reactants:                   ", reactants)
        print("Products:                    ", products)
    # Get the conversion of the ids to formulas
    react_id_form_key = get_ids_to_formulas(reactants, web=web, file_path=file_path)
    prod_id_form_key = get_ids_to_formulas(products, web=web, file_path=file_path)
    if verbose:
        print("Reactant key id to formula:  ", react_id_form_key)
        print("Product key id to formula:   ", prod_id_form_key)

    # Convert the reactants into formulas
    converted_reactants = convert_ids_to_formulas(reactants, react_id_form_key)
    converted_products = convert_ids_to_formulas(products, prod_id_form_key)
    if verbose:
        print("Converted reactants:         ", converted_reactants)
        print("Converted products:          ", converted_products)
    # Convert the formulas into reactants
    react_ele = convert_form_dict_to_elements(converted_reactants)
    prod_ele = convert_form_dict_to_elements(converted_products)
    if verbose:
        print("reactant element dict:       ", react_ele)
        print("product element dict:        ", prod_ele)
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


def check_missing_formulas(eq, web=False, file_path='Data/kegg_data_C.csv.zip'):
    # Convert the Eq in to the dicts
    reactants, products = eq_to_dict(eq)
    ids = list(reactants.keys()) + list(products.keys())
    if web:
        formulas = [get_formula_from_web(id) for id in ids]
    else:
        formulas = get_formulas_from_ids(ids, file_path)

    if len(formulas) != len(ids):
        return True
    else:
        return False


def check_eq_unbalanced(react_ele, prod_ele):
    diff_ele_react, diff_ele_prod = compare_dict_values(react_ele, prod_ele)
    if len(diff_ele_react) + len(diff_ele_prod) > 0:
        return True
    else:
        return False


def inject_compounds(eq_line, missing_r, missing_p, missing_dict):
    # Get the left and right side of the equation
    eq_left, eq_right = map(str.strip, eq_line.split("<=>"))

    # Loop over the missing "elements"
    for item in missing_r:
        eq_left += f" + {missing_dict[item]}"
    for item in missing_p:
        eq_right += f" + {missing_dict[item]}"

    # Make the new equation
    return f"{eq_left} <=> {eq_right}"


def fix_imbalance_core(eq_line, diff_ele_react, diff_ele_prod, inject):
    # Get the reactant and product sides
    eq_left, eq_right = map(str.strip, eq_line.split("<=>"))
    # Find which side has the lowest
    if sum(diff_ele_react.values()) < sum(diff_ele_prod.values()):
        eq_left += f" + {inject}"
    else:
        eq_right += f" + {inject}"
    # Update eq_line with the new equation
    return f"{eq_left} <=> {eq_right}"


def fix_simple_imbalance(eq_line, diff_ele_react, diff_ele_prod):
    """

    """
    # Find the difference in elements
    diff_ele = set(diff_ele_react) | set(diff_ele_prod)
    # Find the difference in values
    diff_val = abs(sum(diff_ele_react.values()) - sum(diff_ele_prod.values()))

    # Attempt to fix issue with missing
    if diff_ele == {"H"} and diff_val % 2 != 0:
        print("Adding H", flush=True)
        return fix_imbalance_core(eq_line, diff_ele_react, diff_ele_prod, "C00080")
    elif diff_ele == {"H"} and diff_val % 2 == 0:
        print("Adding H2", flush=True)
        return fix_imbalance_core(eq_line, diff_ele_react, diff_ele_prod, "C00282")
    elif diff_ele == {"O", "H"}:
        print("Adding H2O", flush=True)
        return fix_imbalance_core(eq_line, diff_ele_react, diff_ele_prod, "C00001")
    elif diff_ele == {"O"}:
        print("Adding O2", flush=True)
        return fix_imbalance_core(eq_line, diff_ele_react, diff_ele_prod, "C00007")
    elif diff_ele == {"C", "O"}:
        print("Adding CO2", flush=True)
        return fix_imbalance_core(eq_line, diff_ele_react, diff_ele_prod, "C00011")
    elif diff_ele == {"N", "H"}:
        print("Adding NH3", flush=True)
        return fix_imbalance_core(eq_line, diff_ele_react, diff_ele_prod, "C00014")
    elif diff_ele == {"C", "H"}:
        print("Adding CH4", flush=True)
        return fix_imbalance_core(eq_line, diff_ele_react, diff_ele_prod, "C01438")
    elif diff_ele == {"N", "O"}:
        print("Adding NO2", flush=True)
        return fix_imbalance_core(eq_line, diff_ele_react, diff_ele_prod, "C00088")
    elif diff_ele == {"S", "O"}:
        print("Adding SO2", flush=True)
        return fix_imbalance_core(eq_line, diff_ele_react, diff_ele_prod, "C09306")
    else:
        print("Could not fix the imbalance", flush=True)
        return eq_line


def check_glycan(eq_line):
    if "G" in eq_line:
        return True
    else:
        return False


def main(eq_file="Data/kegg_data_R.csv.zip"):
    missing_dict = {"H2O": "C00001",
                    "H": "C00080",
                    "Fe": "C00023",  # C14819, C14818 https://www.kegg.jp/entry/C00023
                    "Na": "C01330",
                    "Ca": "C00076",
                    "Cu": "C00070",
                    "Al": "C06264",
                    "Cl": "C00698",
                    "Co": "C00175",
                    "Ni": "C19609",
                    "Mo": "C00150",
                    "O2": "C00007",
                    "H2Se": "C01528",
                    "SeO3": "C05684",
                    "H3PO4": "C00009",  # this is bad up anything smaller?
                    "CO2": "C00011",
                    "NH3": "C00014",
                    "H2": "C00282",
                    "Mn": "C19610",  # https://www.kegg.jp/entry/C00034
                    "Zn": "C00038",
                    "CH4": "C01438",
                    "NO2": "C00088",
                    "CH2O": "C00067",
                    "SO2": "C09306",
                    "S": "C00087",
                    "H2S": "C00283",
                    "I": "C00708",
                    "Br": "C01324",
                    "Pb": "C06696",
                    "K": "C00238",
                    "Mg": "C00305",
                    "F": "C00742",
                    "Te": "C99999",
                    }

    fix_comp_dict = {"H2SO4": "C00059",  # Sulfuric acid
                     "H4P2O7": "C00013",  # Pyrophosphate
                     "C3H7O6P": "C00111",  # Glycerone phosphate
                     "CH3NO": "C00488",  # Formamide
                     "C2H4": "C06547",  # Ethene
                     "H2O2": "C00027",  # Hydrogen peroxide
                     "HCN": "C01326",  # Hydrogen cyanide
                     "H2CO3": "C01353",  # Carbonic acid
                     "C5H9O9P": "C22278",  # 3-Oxoisoapionate 4-phosphate
                     "CH2O": "C00067",  # Formaldehyde
                     "HCl": "C01327",  # Hydrochloric acid
                     "H5P3O10": "C00536",  # Triphosphate
                     "CH5O3P": "C20396",  # Methylphosphonate
                     "SeO3": "C05684",  # Selenite
                     "WH2O4": "C20679",  # Tungstic acid
                     "H2MoO4": "C06232",  # Molybdate
                     "HCO3": "C00288",  # Hydrogencarbonate
                     "H2Se": "C01528",  # Selenous acid
                     }

    out_eq_file = f"{eq_file.split(".")[0]}_processed.csv.zip"
    bad_file = "Data/R_IDs_bad.dat"
    # read the bad file
    with open(bad_file, "r") as f:
        bad_data = f.read()
    bad_ids = bad_data.split("\n")[1:]

    # Load the processed data
    data = pd.read_csv(eq_file)

    # Filter out the bad ids
    print("Filtering out bad ids", flush=True)
    data = data.loc[~data["id"].isin(bad_ids)]

    # Get the data from the dataframe
    ids = data["id"].tolist()
    eq_lines = data["reaction"].tolist()
    ec = data["ec"].tolist()
    print("Data loaded", flush=True)
    print("Data columns", data.columns, flush=True)
    print("Data shape", data.shape, flush=True)
    print("Data head", data.head(4).values, flush=True)

    # Get the size of the data
    n_ids = len(ids)
    print(f"Total number of reactions {n_ids}", flush=True)
    # Init the lists
    bad_n = []
    bad_eq = []
    bad_missing_mol = []
    bad_missing_ele = []
    bad_no_balance = []
    missing_ele = []
    # Define the output file lists
    out_eq_lines = []
    out_ids = []
    out_ec = []

    # Loop over the reactions data
    for i, re_id in enumerate(ids):
        # two injections R04795, R04808
        # if re_id != "R05923":  # R00538, R00634, R00915, R00916, R01317, R01350, R01409
        #     continue
        eq_line = eq_lines[i]
        print(f"\nProcessing {i}/{n_ids} {re_id}", flush=True)
        print("Equation line:", eq_line, flush=True)

        # Check if the equation has n
        if "n" in eq_line:
            print(f"Warning! Equation has n. Skipping ID: {re_id}", flush=True)
            bad_n.append(re_id)
            continue

        # Check if the equation has reactant and product side
        if len(eq_line.split("<=>")) != 2:
            print(f"Warning! Equation does not have a reactant and product side. Skipping ID: {re_id}", flush=True)
            bad_eq.append(re_id)
            continue

        if check_glycan(eq_line):
            print("Warning! Glycan in the equation", flush=True)
            bad_eq.append(re_id)
            continue

        # Check if the compound has missing formulas and skip as it will break the rest of the code
        if check_missing_formulas(eq_line, web=False):
            print("Warning! No formula", flush=True)
            bad_missing_mol.append(re_id)
            continue

        reactants, products, react_ele, prod_ele = get_elements_from_eq(eq_line, verbose=False, web=False)

        if check_missing_elements(react_ele, prod_ele):
            print("Warning! Missing elements", flush=True)
            missing_in_react, missing_in_prod = get_missing_elements(react_ele, prod_ele)
            print("Missing in reactants:        ", missing_in_react, flush=True)
            print("Missing in products:         ", missing_in_prod, flush=True)
            print("Attempting to fix missing products!", flush=True)
            try:
                eq_line = inject_compounds(eq_line, missing_in_react, missing_in_prod, missing_dict)
            except KeyError as e:
                print("No item in the missing dict that could fix; ", e, flush=True)
            # With the new equation line lets try again
            reactants, products, react_ele, prod_ele = get_elements_from_eq(eq_line, verbose=False, web=False)
            if check_missing_elements(react_ele, prod_ele):
                missing_in_react, missing_in_prod = get_missing_elements(react_ele, prod_ele)
                print("Missing in reactants:        ", missing_in_react, flush=True)
                print("Missing in products:         ", missing_in_prod, flush=True)
                bad_missing_ele.append(re_id)
                for val in missing_in_react:
                    missing_ele.append(val)
                for val in missing_in_prod:
                    missing_ele.append(val)
                print("Missing ele set:             ", set(missing_ele), flush=True)
                continue
            else:
                print("Fix worked!", flush=True)

        if check_eq_unbalanced(react_ele, prod_ele):
            print("Warning! unbalanced equation", flush=True)
            # Get the difference in the elements
            diff_ele_react, diff_ele_prod = compare_dict_values(react_ele, prod_ele)
            print("Differences in reactants:    ", diff_ele_react, flush=True)
            print("Differences in products:     ", diff_ele_prod, flush=True)

            try:
                print("Attempt balancing eq x1", flush=True)
                reactants, products = balance_stoichiometry(set(reactants.keys()),
                                                            set(products.keys()),
                                                            underdetermined=None)
                reactants = dict(reactants)
                products = dict(products)
                # Convert the dict back into eq form
                eq_line = get_eq(eq_line, reactants, products)
                print("Rebalance success!", flush=True)

            except:
                print("Could not find stoichiometry on first attempt", flush=True)
                # Attempt a more manual injection to help balance, this simply looks at the pop in-balance
                eq_line = fix_simple_imbalance(eq_line, diff_ele_react, diff_ele_prod)
                print("New eq line:", eq_line, flush=True)
                # Update values
                reactants, products, react_ele, prod_ele = get_elements_from_eq(eq_line, verbose=False, web=False)
                diff_ele_react, diff_ele_prod = compare_dict_values(react_ele, prod_ele)
                if check_eq_unbalanced(react_ele, prod_ele):
                    try:
                        print("Attempt balancing eq x2", flush=True)
                        print("Differences in reactants:    ", diff_ele_react, flush=True)
                        print("Differences in products:     ", diff_ele_prod, flush=True)
                        reactants, products = balance_stoichiometry(set(reactants.keys()),
                                                                    set(products.keys()),
                                                                    underdetermined=None)
                        reactants = dict(reactants)
                        products = dict(products)
                        # Convert the dict back into eq form
                        eq_line = get_eq(eq_line, reactants, products)
                        print("Rebalance success!", flush=True)

                    except:
                        print("Could not find stoichiometry after injection!", flush=True)
                        bad_no_balance.append(re_id)
                        # Skip the iteration
                        continue
                else:
                    print("Rebalance success!", flush=True)

        # Allocate the result to the lists
        out_ids.append(re_id)
        out_eq_lines.append(eq_line)
        out_ec.append(ec[i])

    # Store the data in a dataframe
    df = pd.DataFrame({'id': out_ids, 'reaction': out_eq_lines, 'ec': out_ec})
    # Write the data to a file
    df.to_csv(out_eq_file, compression='zip', encoding='utf-8', index=False)

    # print out the bad files
    print(f"bad n: {bad_n}", flush=True)
    print(f"bad eq: {bad_eq}", flush=True)
    print(f"bad_missing_mol: {bad_missing_mol}", flush=True)
    print(f"bad_missing_ele: {bad_missing_ele}", flush=True)
    print(f"bad no balance: {bad_no_balance}", flush=True)
    # print out the length of each of the lists
    print(f"len bad n: {len(bad_n)}/{n_ids}", flush=True)
    print(f"len bad eq: {len(bad_eq)}/{n_ids}", flush=True)
    print(f"len bad_missing_mol: {len(bad_missing_mol)}/{n_ids}", flush=True)
    print(f"len bad_missing_ele: {len(bad_missing_ele)}/{n_ids}", flush=True)
    print(f"len bad no balance: {len(bad_no_balance)}/{n_ids}", flush=True)
    print(f"missing ele {set(missing_ele)}", flush=True)

    # # Save the data to file
    # bad_df = pd.DataFrame({'bad_n': bad_n,
    #                        'bad_eq': bad_eq,
    #                        'bad_missing_mol': bad_missing_mol,
    #                        'bad_missing_ele': bad_missing_ele,
    #                        'bad_no_balance': bad_no_balance})
    # bad_df.to_csv(bad_file, compression='zip', encoding='utf-8', index=False)


if __name__ == "__main__":
    print("Program started", flush=True)
    main(eq_file="Data/kegg_data_R.csv.zip")
    main(eq_file="Data/atlas_data_kegg_R.csv.zip")
    main(eq_file="Data/atlas_data_R.csv.zip")
    print("Program finished!", flush=True)
