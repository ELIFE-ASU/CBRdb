import re

import pandas as pd
from chempy import balance_stoichiometry


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


def inject_compounds(eq_line, missing_r, missing_p, missing_dict):
    eq_left, eq_right = map(str.strip, eq_line.split("<=>"))
    eq_left += ''.join(f" + {missing_dict[item]}" for item in missing_r)
    eq_right += ''.join(f" + {missing_dict[item]}" for item in missing_p)
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


def main(r_file="Data/kegg_data_R.csv.zip",
         c_file="Data/kegg_data_C.csv.zip",
         bad_file="Data/R_IDs_bad.dat",
         f_fresh=True):
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

    out_eq_file = f"{r_file.split(".")[0]}_processed.csv.zip"

    # read the bad file
    with open(bad_file, "r") as f:
        bad_data = f.read()
    bad_ids = [line.split(',')[0].strip() for line in bad_data.split("\n")[1:]]
    # Load the C data
    data_c = pd.read_csv(c_file)

    # Load the processed data
    data_r = pd.read_csv(r_file)

    # Filter out the bad ids
    print("Filtering out bad ids", flush=True)
    data_r = data_r.loc[~data_r["id"].isin(bad_ids)]

    # Get the data from the dataframe
    ids = data_r["id"].tolist()
    eq_lines = data_r["reaction"].tolist()
    ec = data_r["ec"].tolist()
    print("Data loaded", flush=True)
    print("Data columns", data_r.columns, flush=True)
    print("Data shape", data_r.shape, flush=True)

    # Init the lists
    bad_n = []
    bad_eq = []
    bad_missing_mol = []
    bad_missing_ele = []
    bad_no_balance = []
    missing_ele = []

    # Checks if you want to start fresh or not
    if not f_fresh:
        # Load the data from the file
        df = pd.read_csv(out_eq_file)
        out_ids = df["id"].tolist()
        out_eq_lines = df["reaction"].tolist()
        out_ec = df["ec"].tolist()
        print("Data loaded", flush=True)
    else:
        # Define the output file lists
        out_eq_lines = []
        out_ids = []
        out_ec = []

    # Get the size of the data
    n_ids = len(ids)
    print(f"Total number of reactions {n_ids}", flush=True)

    # Loop over the reactions data
    for i, re_id in enumerate(ids):
        # Skip the reactions that have already been processed
        if re_id in out_ids:
            continue

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

        # Check if the compound has missing formulas and skip as it will break the rest of the code
        if check_missing_formulas(eq_line, data_c):
            print("Warning! No formula", flush=True)
            bad_missing_mol.append(re_id)
            continue

        reactants, products, react_ele, prod_ele = get_elements_from_eq(eq_line, data_c)
        print("Reactants: ", reactants, flush=True)
        print("Products:  ", products, flush=True)

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
            reactants, products, react_ele, prod_ele = get_elements_from_eq(eq_line, data_c)
            if check_missing_elements(react_ele, prod_ele):
                missing_in_react, missing_in_prod = get_missing_elements(react_ele, prod_ele)
                print("Missing in reactants:        ", missing_in_react, flush=True)
                print("Missing in products:         ", missing_in_prod, flush=True)
                bad_missing_ele.append(re_id)
                for val in missing_in_react:
                    missing_ele.append(val)
                for val in missing_in_prod:
                    missing_ele.append(val)
                continue
            else:
                print("Fix worked!", flush=True)

        if check_eq_unbalanced(react_ele, prod_ele):
            print("Warning! unbalanced equation", flush=True)
            # Get the difference in the elements
            diff_ele_react, diff_ele_prod = compare_dict_values(react_ele, prod_ele)
            print("Differences in reactants:    ", diff_ele_react, flush=True)
            print("Differences in products:     ", diff_ele_prod, flush=True)

            # Find the difference in elements
            diff_ele = set(diff_ele_react) | set(diff_ele_prod)
            # Find the difference in values
            diff_val = abs(sum(diff_ele_react.values()) - sum(diff_ele_prod.values()))
            # Try injecting a H to help balance
            if diff_ele == {"H"} and diff_val == 1:
                print("Adding H", flush=True)
                eq_line = fix_imbalance_core(eq_line, diff_ele_react, diff_ele_prod, "C00080")
                reactants, products, react_ele, prod_ele = get_elements_from_eq(eq_line, data_c)
                if check_eq_unbalanced(react_ele, prod_ele):
                    diff_ele_react, diff_ele_prod = compare_dict_values(react_ele, prod_ele)
                else:
                    # Allocate the result to the lists
                    out_ids.append(re_id)
                    out_eq_lines.append(eq_line)
                    out_ec.append(ec[i])
                    continue

            try:
                print("Attempt balancing eq x1", flush=True)
                reactants, products = balance_stoichiometry(set(reactants.keys()),
                                                            set(products.keys()),
                                                            underdetermined=None)
                reactants = dict(reactants)
                products = dict(products)
                # Convert the dict back into eq form
                eq_line = get_eq(eq_line, reactants, products, data_c)
                print("Rebalance success!", flush=True)

            except:
                print("Could not find stoichiometry on first attempt", flush=True)
                # Attempt a more manual injection to help balance, this simply looks at the pop in-balance
                eq_line = fix_simple_imbalance(eq_line, diff_ele_react, diff_ele_prod)
                print("New eq line:", eq_line, flush=True)
                # Update values
                reactants, products, react_ele, prod_ele = get_elements_from_eq(eq_line, data_c)
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
                        eq_line = get_eq(eq_line, reactants, products, data_c)
                        print("Rebalance success!", flush=True)
                        print("New eq line:", eq_line, flush=True)

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

    # check if the data is fresh
    if not f_fresh:
        # Load the data from the file
        df_old = pd.read_csv(out_eq_file)
        # Append the data
        df = pd.concat([df_old, df])
        # Drop the duplicates
        df = df.drop_duplicates(subset="id", keep="last")
        # # Reset the index
        # df = df.reset_index(drop=True)
        # sort the data
        df = df.sort_values(by="id")

    # Write the data to a file
    # get the shape of the data
    print("Data shape", df.shape, flush=True)
    df.to_csv(out_eq_file, compression='zip', encoding='utf-8', index=False)

    # print out the bad files
    print(f"bad n: {bad_n}", flush=True)
    print(f"bad eq: {bad_eq}", flush=True)
    print(f"bad_missing_mol: {bad_missing_mol}", flush=True)
    print(f"bad_missing_ele: {bad_missing_ele}", flush=True)
    print(f"bad no balance: {bad_no_balance}", flush=True)
    # print out the length of each of the lists
    n_fail = len(bad_n) + len(bad_eq) + len(bad_missing_mol) + len(bad_missing_ele) + len(bad_no_balance)
    print(f"Total failed reactions: {n_fail}/{n_ids}", flush=True)
    print(f"len bad n: {len(bad_n)}/{n_ids}", flush=True)
    print(f"len bad eq: {len(bad_eq)}/{n_ids}", flush=True)
    print(f"len bad_missing_mol: {len(bad_missing_mol)}/{n_ids}", flush=True)
    print(f"len bad_missing_ele: {len(bad_missing_ele)}/{n_ids}", flush=True)
    print(f"len bad no balance: {len(bad_no_balance)}/{n_ids}", flush=True)
    print(f"missing ele {set(missing_ele)}", flush=True)


if __name__ == "__main__":
    print("Program started", flush=True)
    main(r_file="Data/kegg_data_R.csv.zip")
    main(r_file="Data/atlas_data_kegg_R.csv.zip")
    main(r_file="Data/atlas_data_R.csv.zip")
    print("Program finished!", flush=True)
