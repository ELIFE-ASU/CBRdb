import os
import time

import pandas as pd
# import dask.dataframe as dd
import swifter
from chempy import balance_stoichiometry

from .tools_eq import (get_eq,
                       get_elements_from_eq,
                       compare_dict_values,
                       check_missing_elements,
                       check_eq_unbalanced,
                       get_missing_elements,
                       standardise_eq,
                       check_contains_var_list,
                       check_missing_formulas,
                       full_check_eq_unbalanced,
                       rebalance_eq,
                       fix_imbalance_core,
                       inject_compounds,
                       )
from .tools_mols import (get_small_compounds, get_small_compounds_all, get_compounds_with_matching_elements)


def fix_simple_imbalance(eq_line, diff_ele_react, diff_ele_prod):
    """
    Attempts to fix the imbalance in a reaction equation by injecting common compounds based on the element differences.

    Parameters:
    eq_line (str): The original reaction equation line.
    diff_ele_react (dict): A dictionary of element differences on the reactant side.
    diff_ele_prod (dict): A dictionary of element differences on the product side.

    Returns:
    str: The updated reaction equation with the injected compound if a fix is found, otherwise the original equation.
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


def kitchen_sink(eq, data_c, small_compounds):
    reactants, products, react_ele, prod_ele = get_elements_from_eq(eq, data_c)
    diff_ele_react, diff_ele_prod = compare_dict_values(react_ele, prod_ele)

    compounds = get_compounds_with_matching_elements(small_compounds, diff_ele_react, diff_ele_prod)
    print("Compounds that might match:  ", compounds, flush=True)

    eq_new = eq
    attempt = 1
    for compound in compounds:
        # Counter
        print(f"Attempt {attempt}", flush=True)
        attempt += 1

        # Get the elements
        _, _, react_ele, prod_ele = get_elements_from_eq(eq_new, data_c)
        # Check the difference
        diff_ele_react, diff_ele_prod = compare_dict_values(react_ele, prod_ele)
        print("Differences in reactants:    ", diff_ele_react, flush=True)
        print("Differences in products:     ", diff_ele_prod, flush=True)

        # Inject the compound
        eq_new = fix_imbalance_core(eq_new, diff_ele_react, diff_ele_prod, compound)
        # rebalance the equation
        eq_new = rebalance_eq(eq_new, data_c)
        if eq_new is False:
            eq_new = eq
        else:
            break

    return eq_new


def fix_reactions_data(r_file="../data/kegg_data_R.csv",
                       c_file="../data/kegg_data_C.csv",
                       bad_file="../data/R_IDs_bad.dat",
                       f_fresh=True,
                       f_assume_var=True,
                       f_save_intermediate=False,
                       f_parallel=True):
    """creates processed version of r_file csv from original r_file csv; returns processed r_file as pd.DataFrame"""
    # f_assume_var=True => assume that the equation contains a var list are correct

    if f_parallel:
        swifter.set_defaults(allow_dask_on_strings=True, force_parallel=True, progress_bar=False)

    # Get the absolute paths
    r_file = os.path.abspath(r_file)
    c_file = os.path.abspath(c_file)
    bad_file = os.path.abspath(bad_file)

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

    out_eq_file = f"{r_file.split('.')[0]}_processed.csv".replace('_deduped', '')

    # Read the bad reactions file
    with open(bad_file, "r") as f:
        bad_data = f.read()
    bad_ids = [line.split(',')[0].strip() for line in bad_data.split("\n")[1:]]

    # Load the processed compound data
    print("Loading the compound data...", flush=True)
    data_c = pd.read_csv(c_file)
    print(f"Compound data shape: {data_c.shape}", flush=True)

    # Load the small compounds
    data_c_1 = get_small_compounds(c_path=c_file, n=1)
    data_c_2 = get_small_compounds(c_path=c_file, n=2)
    data_c_3 = get_small_compounds(c_path=c_file, n=3)
    # print the size of the data
    print(f"Small compound size: 1:{data_c_1.shape}, 2:{data_c_2.shape}, 3:{data_c_3.shape}", flush=True)

    data_c_small = get_small_compounds_all(c_path=c_file, n=1)

    # Load the processed reaction data
    print("Loading the reaction data...", flush=True)
    data_r = pd.read_csv(r_file)
    print("data loaded", flush=True)
    print("data columns", data_r.columns, flush=True)
    print("data shape", data_r.shape, flush=True)

    # Sort by the index
    data_r = data_r.sort_values(by="id")

    # Filter out the bad ids, data that is shortcut/general/incomplete
    print("Filtering out bad ids", flush=True)
    bool_bad_ids = data_r["id"].isin(bad_ids)
    data_r = data_r.loc[~bool_bad_ids]
    print(f"Number of bad ids removed: {sum(bool_bad_ids)}", flush=True)

    # Filter out the data that has missing formulas
    print("Filtering out missing formulas", flush=True)
    if f_parallel:
        t0 = time.time()
        bool_missing_data = data_r['reaction'].swifter.force_parallel(enable=True).allow_dask_on_strings(
            enable=True).apply(check_missing_formulas, args=(data_c,))
        print("Time to check missing formulas: ", time.time() - t0)
    else:
        t0 = time.time()
        bool_missing_data = data_r['reaction'].apply(check_missing_formulas, args=(data_c,))
        print("Time to check missing formulas: ", time.time() - t0)

    data_r_missing_data = data_r[bool_missing_data]
    if f_save_intermediate:
        data_r_missing_data.to_csv(f"{r_file.split('.')[0]}_missing_data.csv",
                                   encoding='utf-8',
                                   index=False)
    data_r = data_r[~bool_missing_data]
    print(f"Number of missing formulas removed: {sum(bool_missing_data)}", flush=True)

    # Filter out the reactions that contain a var list
    print("Filtering out var list", flush=True)
    if f_parallel:
        t0 = time.time()
        bool_var_list = data_r['reaction'].swifter.force_parallel(enable=True).apply(check_contains_var_list,
                                                                                     args=(data_c,))
        print("Time to check missing formulas: ", time.time() - t0)
    else:

        t0 = time.time()
        bool_var_list = data_r['reaction'].apply(check_contains_var_list, args=(data_c,))
        print("Time to check missing formulas: ", time.time() - t0)

    data_r_var_list = data_r[bool_var_list]
    if f_save_intermediate:
        data_r_var_list.to_csv(f"{r_file.split('.')[0]}_var_list.csv",
                               encoding='utf-8',
                               index=False)
    data_r = data_r[~bool_var_list]
    print(f"Number of var list reactions removed: {sum(bool_var_list)}", flush=True)

    # Filter out the data that is not balanced
    print("Filtering out unbalanced reactions", flush=True)
    if f_parallel:
        t0 = time.time()
        bool_unbalanced = data_r['reaction'].swifter.force_parallel(enable=True).apply(full_check_eq_unbalanced,
                                                                                       args=(data_c,))
        print("Time to check missing formulas: ", time.time() - t0)
    else:
        t0 = time.time()
        bool_unbalanced = data_r['reaction'].apply(full_check_eq_unbalanced, args=(data_c,))
        print("Time to check missing formulas: ", time.time() - t0)

    # Get the data that is unbalanced
    data_r_unbalanced = data_r[bool_unbalanced]
    if f_save_intermediate:
        data_r_unbalanced.to_csv(f"{r_file.split('.')[0]}_unbalanced.csv",
                                 encoding='utf-8',
                                 index=False)
    # Get the data that is balanced
    data_r = data_r[~bool_unbalanced]
    # Determine the number of reactions that have been removed
    print(f"Number of unbalanced reactions: {sum(bool_unbalanced)}", flush=True)

    # # Filter out the var data that is not balanced
    # print("Filtering out unbalanced reactions", flush=True)
    # bool_unbalanced = data_r['reaction'].apply(full_check_eq_unbalanced, args=(data_c,))
    # # Get the data that is unbalanced
    # data_r_unbalanced = data_r[bool_unbalanced]
    # # Get the data that is balanced
    # data_r = data_r[~bool_unbalanced]
    # # Determine the number of reactions that have been removed
    # print(f"Number of unbalanced reactions: {sum(bool_unbalanced)}", flush=True)
    # print(data_r_unbalanced, flush=True)
    # print(data_r_unbalanced["id"].item, flush=True)

    # Get the data from the unbalanced dataframe
    ids = data_r_unbalanced["id"].tolist()
    eq_lines = data_r_unbalanced["reaction"].tolist()
    n_ids = len(ids)

    ids_out = []
    eq_lines_out = []

    # Loop over the unbalanced reactions data
    for i in range(n_ids):
        # Standardise the equation
        eq_line = standardise_eq(eq_lines[i])
        id = ids[i]
        print(f"\nProcessing {i}/{n_ids} {id}", flush=True)
        print("Equation line:", eq_line, flush=True)
        reactants, products, react_ele, prod_ele = get_elements_from_eq(eq_line, data_c)
        print("Reactants: ", reactants, flush=True)
        print("Products:  ", products, flush=True)
        eq_line_new = rebalance_eq(eq_line, data_c)
        if eq_line_new is False:
            # Need to get dirty and try to fix the imbalance
            eq_line_new = kitchen_sink(eq_line, data_c, data_c_1)
            if eq_line_new is False:
                print("Could not fix the imbalance", flush=True)
        else:
            ids_out.append(id)
            eq_lines_out.append(eq_line_new)

    data_r_rebalanced = data_r_unbalanced

    print("Combining the data", flush=True)
    # Combine the data
    if f_assume_var:
        # Here we have assumed that the data_r_var_list reactions data is correct
        # This is questionable as the data may be incorrect
        print("Assuming equations with a var list data is correct...", flush=True)
        df_final = pd.concat([data_r, data_r_var_list, data_r_rebalanced])
    else:
        df_final = pd.concat([data_r, data_r_rebalanced])

    # Get the finial length of the data
    print(f"Final data shape: {df_final.shape}", flush=True)

    # Sort by the index
    df_final = df_final.sort_values(by="id")
    # Write the data to a file
    df_final.to_csv(out_eq_file, encoding='utf-8', index=False)

    exit()

    # Get the data from the dataframe
    ids = data_r["id"].tolist()
    eq_lines = data_r["reaction"].tolist()
    ec = data_r["ec"].tolist()

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
        out_ids = df["index"].tolist()
        out_eq_lines = df["reaction"].tolist()
        out_ec = df["ec"].tolist()
        print("data loaded", flush=True)
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
                print("Reactants: ", reactants, flush=True)
                print("Products:  ", products, flush=True)
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

        # Standardise the equation
        eq_line = standardise_eq(eq_line)

        # Allocate the result to the lists
        out_ids.append(re_id)
        out_eq_lines.append(eq_line)
        out_ec.append(ec[i])

    # Store the data in a dataframe
    df = pd.DataFrame({'index': out_ids, 'reaction': out_eq_lines, 'ec': out_ec})

    # check if the data is fresh
    if not f_fresh:
        # Load the data from the file
        df_old = pd.read_csv(out_eq_file)
        # Append the data
        df = pd.concat([df_old, df])
        # Drop the duplicates
        df = df.drop_duplicates(subset="index", keep="last")
        # # Reset the index
        # df = df.reset_index(drop=True)
        # sort the data
        df = df.sort_values(by="index")

    # Write the data to a file
    # get the shape of the data
    print("data shape", df.shape, flush=True)
    df.to_csv(out_eq_file, encoding='utf-8', index=False)

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

    return df
