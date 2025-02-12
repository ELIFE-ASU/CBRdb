import os
import time

import pandas as pd
import swifter

from .tools_eq import (get_elements_from_eq,
                       compare_dict_values,
                       standardise_eq,
                       check_contains_var_list,
                       check_missing_formulas,
                       full_check_eq_unbalanced,
                       rebalance_eq,
                       fix_imbalance_core,
                       )
from .tools_mols import (get_small_compounds, get_compounds_with_matching_elements)


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
        print(f"Attempt {attempt}, {compound}", flush=True)
        attempt += 1

        # Get the elements
        _, _, react_ele, prod_ele = get_elements_from_eq(eq_new, data_c)
        # Check the difference
        diff_ele_react, diff_ele_prod = compare_dict_values(react_ele, prod_ele)
        print("Differences in reactants:    ", diff_ele_react, flush=True)
        print("Differences in products:     ", diff_ele_prod, flush=True)

        # Inject the compound
        eq_new = fix_imbalance_core(eq_new, diff_ele_react, diff_ele_prod, compound)
        print("New equation:                ", eq_new, flush=True)

        # Get the elements
        _, _, react_ele, prod_ele = get_elements_from_eq(eq_new, data_c)
        # Check the difference
        diff_ele_react, diff_ele_prod = compare_dict_values(react_ele, prod_ele)
        print("Differences in reactants:    ", diff_ele_react, flush=True)
        print("Differences in products:     ", diff_ele_prod, flush=True)

        # Check if the equation is balanced by checking the differences
        if not diff_ele_react and not diff_ele_prod:
            print("Balanced equation found", flush=True)
            return eq_new

        # rebalance the equation
        eq_new = rebalance_eq(eq_new, data_c)
        if eq_new is False:
            eq_new = eq
        else:
            break

    return eq_new


def dict_ele_contains_star(react_ele, prod_ele):
    """
    Checks if there is a '*' in either reactants or products.

    Parameters:
    react_ele (dict): A dictionary of reactant elements.
    prod_ele (dict): A dictionary of product elements.

    Returns:
    bool: True if there is a '*' in either reactants or products, False otherwise.
    """
    return '*' in react_ele or '*' in prod_ele


def fix_reactions_data(r_file="../data/kegg_data_R.csv",
                       c_file="../data/kegg_data_C.csv",
                       bad_file="../data/R_IDs_bad.dat",
                       f_fresh=True,
                       f_assume_var=True,
                       f_assume_star=True,
                       f_save_intermediate=False,
                       f_parallel=True):
    """creates processed version of r_file csv from original r_file csv; returns processed r_file as pd.DataFrame"""
    # f_assume_var=True => assume that the equation contains a var list are correct
    # f_assume_star=True => assume that the equation contains a star are correct

    if f_parallel:
        swifter.set_defaults(allow_dask_on_strings=True, force_parallel=True, progress_bar=False)

    # Get the absolute paths
    r_file = os.path.abspath(r_file)
    c_file = os.path.abspath(c_file)
    bad_file = os.path.abspath(bad_file)

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
    print(f"Small compound size: 1:{data_c_1.shape}, 2:{data_c_2.shape}, 3:{data_c_3.shape}", flush=True)

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

        # Check if the equation contains a '*' in either reactants or products
        if dict_ele_contains_star(react_ele, prod_ele):
            # We assume that the equation is correct and add it to the output
            if f_assume_star:
                print(f"Assuming {id} is correct", flush=True)
                ids_out.append(id)
                eq_lines_out.append(eq_line)
                continue
            # We assume that the equation is incorrect and skip it as we cannot fix it
            else:
                print(f"Skipping {id}", flush=True)
                continue

        if id == 'R00263':
            print(f"Skipping {id}", flush=True)
            continue

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
