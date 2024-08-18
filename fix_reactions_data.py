import os
import re
import time

import pandas as pd
import requests
from chempy import balance_stoichiometry
from requests import Session
from requests.adapters import HTTPAdapter
from urllib3.util import Retry

"""
things to do
- check if the equation has n
- check if the equation has a reactant and product side
- check if the equation balances

fucked buckets
1) Has an n
2) Missing mol file (no formula)
3) Equations do not balance
"""


def preprocess_kegg_r(target_dir, outfile):
    # Get a list of all files in the directory
    paths = file_list_all(target_dir)
    N = len(paths)
    id_list = []
    eq_list = []
    ec_list = []

    # Loop over the reactions data
    for i, path in enumerate(paths):
        # Get the ID
        re_id = os.path.basename(path).split(".")[0]
        if i % 100 == 0:
            print(f"Processing {i}/{N} {re_id}", flush=True)
        # Load the data
        with open(path, "r") as f:
            data = f.read()
            # Split the data by new lines
            data = data.split("\n")
        # Get the line which contains the equation
        eq_line = [d for d in data if "EQUATION" in d][0].split("EQUATION")[1].strip()
        # Get the line which contains the enzyme class
        try:
            ec_line = [d for d in data if "ENZYME" in d][0].split("ENZYME")[1].strip()
        except:
            ec_line = " "
        # Append the data to the lists
        id_list.append(re_id)
        eq_list.append(eq_line)
        ec_list.append(ec_line)
    # Make the dataframe for the id and the equation
    df = pd.DataFrame({'ID': id_list, 'Reaction': eq_list})
    # Write the data to a file
    df.to_csv(outfile, compression='zip', encoding='utf-8')
    return None


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
        return "fucked"
    # Check if the response is ok
    if response.ok:
        # Strip the last section of the file
        res = response.text.split("> <ENTRY>")[0]
        data = res.split("\n")
        # Get the line which contains the formula
        data = [line for line in data if "FORMULA" in line]
        if len(data) == 0:
            return "fucked"
        return data[0].split("FORMULA")[1].strip()
    else:
        # Some error in the response
        print(f"Error in ID {id}, response {response.status_code}")
        return "fucked"


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
    return re.findall(r'[A-Z][^A-Z]*', input_string)


def convert_formula_to_dict(input_string):
    parts = split_by_letters(input_string)
    result = {}
    for part in parts:
        result[part[0]] = int(part[1:] or '1')
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


def get_formulas_from_ids(ids, file_path):
    # "Data/kegg_data_C.csv.zip"
    data = pd.read_csv(file_path, index_col=0)
    formulas = data.loc[data["CID"].isin(ids), "Formula"]
    if formulas.empty:
        return None
    elif len(formulas) == 1:
        return formulas.tolist()[0]
    else:
        return formulas.tolist()


if __name__ == "__main__":
    print("Program started", flush=True)
    # in1 = "C10H16N5O13P3"
    # in2 = "C10H15N5O10P2"
    # out1 = convert_formula_to_dict(in1)
    # out2 = convert_formula_to_dict(in2)
    # print(out1,flush=True)
    # print(out2,flush=True)
    #
    # print(multiply_dict(out1, 2),flush=True)
    #
    # diff_in_dict1, diff_in_dict2 = compare_dict_values(out1, out2)
    # print(f"Differences in dict1: {diff_in_dict1}",flush=True)
    # print(f"Differences in dict2: {diff_in_dict2}",flush=True)
    #
    # exit()
    f_preprocess = False
    target_dir = r"..\data\kegg_data_R"

    eq_file = "Data/kegg_data_R_eq.csv.zip"
    reac_file = "Data/kegg_data_R_eq.csv.zip"
    if f_preprocess:
        preprocess_kegg_r(target_dir, reac_file)
        print("Preprocessing done", flush=True)

    # Load the processed data
    data = pd.read_csv(eq_file, index_col=0)
    # Get the IDs and the formulas
    ids = data["ID"].tolist()
    formulas = data["Reaction"].tolist()
    # Get the size of the data
    N = len(ids)
    print(f"Total number of reactions {N}", flush=True)
    fucked_n = []
    fucked_eq = []
    fucked_missing_mol = []
    fucked_no_balance = []
    # Loop over the reactions data
    for i, re_id in enumerate(ids):
        if i < 3:
            continue
        print(f"\nProcessing {i}/{N} {re_id}", flush=True)
        eq_line = formulas[i]
        print("Equation line:", eq_line, flush=True)

        # Check if the equation has n
        if "n" in eq_line:
            print(f"Warning! Equation has n. Skipping ID: {re_id}", flush=True)
            fucked_n.append(re_id)
            continue
        # Check if the equation has reactant and product side
        eq_sides = eq_line.split("<=>")
        if len(eq_sides) != 2:
            print(f"Warning! Equation does not have a reactant and product side. Skipping ID: {re_id}", flush=True)
            fucked_eq.append(re_id)
            continue

        tmp1 = eq_sides[0].split("+")
        print(tmp1, flush=True)

        # Get all the numbers in front of the equation line
        tmp1 = eq_sides[0].split()
        print(tmp1, flush=True)
        # count the number of molecules
        n_mols = len([i for i in tmp1 if str(i) == "+"]) + 1
        filtered_array = select_numbers(tmp1)
        print(filtered_array, flush=True)

        all_ids = [id for id in eq_sides[0].split(" ") if id.startswith("C")]
        eq_react = [get_formula_from_web(id) for id in all_ids]
        print("Reactant", eq_react, flush=True)
        if "fucked" in eq_react:
            print("Warning! No formula in reactant")
            fucked_missing_mol.append(re_id)
            continue

        all_ids = [id for id in eq_sides[1].split(" ") if id.startswith("C")]
        eq_prod = [get_formula_from_web(id) for id in all_ids]
        print("Product", eq_prod, flush=True)
        if "fucked" in eq_prod:
            print("Warning! No formula in product")
            fucked_missing_mol.append(re_id)
            continue
        # Trying to balance_stoichiometry
        try:
            reac, prod = balance_stoichiometry(eq_react, eq_prod, underdetermined=None)
            print(dict(reac))
            print(dict(prod))
        except:
            print("Could not find stoichiometry")
            fucked_no_balance.append(re_id)

    # print out the bad files
    print(f"fucked n: {fucked_n}")
    print(f"fucked eq: {fucked_eq}")
    print(f"fucked_missing_mol: {fucked_missing_mol}")
    print(f"fucked no balance: {fucked_no_balance}")
    # print out the length of each of the lists
    print(f"len fucked n: {len(fucked_n)}")
    print(f"len fucked eq: {len(fucked_eq)}")
    print(f"len fucked_missing_mol: {len(fucked_missing_mol)}")
    print(f"len fucked no balance: {len(fucked_no_balance)}")
