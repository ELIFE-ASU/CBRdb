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


def print_and_log(statement, file):
    """
    Prints and logs a statement.

    Parameters:
    statement (str): The statement to be printed and logged.
    file (file object): The file object to log the statement.

    Returns:
    None
    """
    print(statement, flush=True)
    file.write(statement + "\n")
    return None


def kitchen_sink(eq, data_c, small_compounds):
    """
    Attempts to balance a chemical equation by injecting small compounds based on element differences.

    Parameters:
    eq (str): The original chemical equation.
    data_c (pd.DataFrame): The DataFrame containing compound data.
    small_compounds (list): A list of small compounds to be used for balancing.

    Returns:
    str: The balanced chemical equation if a solution is found, otherwise the original equation.
    """
    reactants, products, react_ele, prod_ele = get_elements_from_eq(eq, data_c)
    diff_ele_react, diff_ele_prod = compare_dict_values(react_ele, prod_ele)

    compounds = get_compounds_with_matching_elements(small_compounds, diff_ele_react, diff_ele_prod)
    print(f"Compounds that might match:  {compounds}", flush=True)

    if not compounds:
        print("No compounds found that match the element differences.", flush=True)
        return False

    eq_new = eq
    attempt = 1
    for compound in compounds:
        # Counter
        print(f"Attempt {attempt}, {compound}", flush=True)

        # Get the elements
        _, _, react_ele, prod_ele = get_elements_from_eq(eq_new, data_c)
        # Check the difference
        diff_ele_react, diff_ele_prod = compare_dict_values(react_ele, prod_ele)
        print(f"Differences in reactants:    {diff_ele_react}", flush=True)
        print(f"Differences in products:     {diff_ele_prod}", flush=True)

        # Inject the compound
        eq_new = fix_imbalance_core(eq_new, diff_ele_react, diff_ele_prod, compound)
        print(f"New equation:                {eq_new}", flush=True)

        # Get the elements
        _, _, react_ele, prod_ele = get_elements_from_eq(eq_new, data_c)

        # Check the difference
        diff_ele_react, diff_ele_prod = compare_dict_values(react_ele, prod_ele)
        print(f"Differences in reactants:    {diff_ele_react}", flush=True)
        print(f"Differences in products:     {diff_ele_prod}", flush=True)

        # Check if the equation is balanced by checking the differences
        if not diff_ele_react and not diff_ele_prod:
            print("Differences null, balanced equation found!", flush=True)
            # Return the new equation line
            print(f"Final equation:              {eq_new}", flush=True)
            return eq_new

        # Try to rebalance the equation
        eq_new = rebalance_eq(eq_new, data_c)
        if eq_new is False:
            # Revert to the original equation and try the next compound
            eq_new = eq
            print("Rebalance failed, trying next compound...", flush=True)
        else:
            # Return the new equation line
            print(f"Final equation:              {eq_new}", flush=True)
            return eq_new
        # Increment the attempt
        attempt += 1
    return False


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
                       log_file="fix_reactions_data_log.dat",
                       rebalance_file="R_IDs_bad_rebalance.dat",
                       f_assume_var=True,
                       f_assume_star=True,
                       f_save_intermediate=False,
                       f_parallel=True,
                       rebalance_depth=1):
    if f_parallel:
        swifter.set_defaults(allow_dask_on_strings=True, force_parallel=True, progress_bar=False)

    # Get the absolute paths
    r_file = os.path.abspath(r_file)
    c_file = os.path.abspath(c_file)
    bad_file = os.path.abspath(bad_file)

    # Get the path from r_file
    tmp_path = os.path.dirname(r_file)

    # Get the file names
    log_file = os.path.basename(r_file).split('.')[0] + f'_{log_file}'
    rebalance_file = os.path.basename(r_file).split('.')[0] + f'_{rebalance_file}'

    # Get the absolute paths
    log_file = os.path.abspath(os.path.join(tmp_path, log_file))
    rebalance_file = os.path.abspath(os.path.join(tmp_path, rebalance_file))

    # Open the log file and rebalance file
    f_log = open(log_file, "w")
    f_rebalance = open(rebalance_file, "w")

    # Write the header
    f_log.write("# This file contains a log of information on the rebalancer run\n")
    f_rebalance.write("# Bad IDs, reason\n")

    # Get the output file name
    out_eq_file = f"{r_file.split('.')[0]}_processed.csv".replace('_deduped', '')

    # Read the bad reactions file
    with open(bad_file, "r") as f:
        bad_data = f.read()
    bad_ids = [line.split(',')[0].strip() for line in bad_data.split("\n")[1:]]

    # Load the processed compound data
    print_and_log("Loading the compound data...", f_log)
    data_c = pd.read_csv(c_file)
    print_and_log(f"Compound data shape: {data_c.shape}", f_log)

    # Load the small compounds
    data_c_1 = get_small_compounds(c_path=c_file, n=1)
    data_c_2 = get_small_compounds(c_path=c_file, n=2)
    data_c_3 = get_small_compounds(c_path=c_file, n=3)
    print_and_log(f"Small compound size: 1:{data_c_1.shape}, 2:{data_c_2.shape}, 3:{data_c_3.shape}", f_log)

    # Load the processed reaction data
    print("Loading the reaction data...", flush=True)
    data_r = pd.read_csv(r_file)
    print("data loaded", flush=True)
    print(f"data columns: {data_r.columns}", flush=True)
    print(f"data shape: {data_r.shape}", flush=True)

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
        print(f"Time to check missing formulas: {time.time() - t0}", flush=True)
    else:
        t0 = time.time()
        bool_missing_data = data_r['reaction'].apply(check_missing_formulas, args=(data_c,))
        print(f"Time to check missing formulas: {time.time() - t0}", flush=True)

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
        print(f"Time to check missing formulas: {time.time() - t0}", flush=True)
    else:
        t0 = time.time()
        bool_var_list = data_r['reaction'].apply(check_contains_var_list, args=(data_c,))
        print(f"Time to check missing formulas: {time.time() - t0}", flush=True)

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
        print(f"Time to check missing formulas: {time.time() - t0}", flush=True)
    else:
        t0 = time.time()
        bool_unbalanced = data_r['reaction'].apply(full_check_eq_unbalanced, args=(data_c,))
        print(f"Time to check missing formulas: {time.time() - t0}", flush=True)

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

    # Get the data from the unbalanced dataframe
    ids = data_r_unbalanced["id"].tolist()
    eq_lines = data_r_unbalanced["reaction"].tolist()
    n_ids = len(ids)

    # Initialise the output lists
    ids_out = []
    eq_lines_out = []
    ids_failed = []

    # Loop over the unbalanced reactions data
    for i in range(n_ids):
        # Enforce that the equation is standardised
        eq_line = standardise_eq(eq_lines[i])
        # Get the id
        id = ids[i]

        print(f"\nProcessing {i}/{n_ids} {id}", flush=True)
        print(f"Original equation line: {eq_line}", flush=True)
        reactants, products, react_ele, prod_ele = get_elements_from_eq(eq_line, data_c)
        print(f"Reactants: {reactants}", flush=True)
        print(f"Products:  {products}", flush=True)

        # Check if the equation contains a '*' in either reactants or products
        if dict_ele_contains_star(react_ele, prod_ele):
            # We assume that the equation is correct and add it to the output
            if f_assume_star:
                print(f"Assuming eq {id} is correct", flush=True)
                ids_out.append(id)
                eq_lines_out.append(eq_line)
                continue
            # We assume that the equation is incorrect and skip it as we cannot fix it
            else:
                print(f"Assume that the eq {id} is incorrect and skip it as we cannot fix it..", flush=True)
                continue

        if id == 'R00263':
            print(f"Skipping {id} due to it being malformed", flush=True)
            continue

        # Check if the equation is balanced
        eq_line_new = rebalance_eq(eq_line, data_c)
        # Equation is not balance...
        if eq_line_new is False:
            # Attempt injecting methods with the c1 compounds
            eq_line_new = kitchen_sink(eq_line, data_c, data_c_1)
            if eq_line_new is False:
                print(f"Could not fix the imbalance for eq {id} using c1 list!", flush=True)
                if rebalance_depth > 1:
                    # Attempt with the c2 compounds
                    eq_line_new = kitchen_sink(eq_line, data_c, data_c_2)
                    if eq_line_new is False:
                        print(f"Could not fix the imbalance for eq {id} using c2 list!", flush=True)
                        if rebalance_depth > 2:
                            # Attempt with the c3 compounds
                            eq_line_new = kitchen_sink(eq_line, data_c, data_c_3)
                            if eq_line_new is False:
                                print(f"Could not fix the imbalance for eq {id} using c3 list!", flush=True)
                                ids_failed.append(id)
                                continue
                        else:
                            ids_failed.append(id)
                            continue
                else:
                    ids_failed.append(id)
                    continue

        print(f"rebalanced {id} with new equation line: {eq_line_new}", flush=True)
        ids_out.append(id)
        eq_lines_out.append(eq_line_new)

    print(f"Number of reactions rebalanced: {len(ids_out)}", flush=True)
    print(f"Number of reactions failed: {len(ids_failed)}", flush=True)
    print(f"Reactions failed: {ids_failed}", flush=True)

    # Update the "reaction" column in data_r_unbalanced using eq_lines_out and ids_out
    data_r_unbalanced.loc[data_r_unbalanced["id"].isin(ids_out), "reaction"] = eq_lines_out

    data_r_rebalanced = data_r_unbalanced

    print("Combining the data!", flush=True)
    print(f"data_r shape: {data_r.shape}", flush=True)
    print(f"data_r_var_list shape: {data_r_var_list.shape}", flush=True)
    print(f"data_r_rebalanced shape: {data_r_rebalanced.shape}", flush=True)

    # Combine the data
    if f_assume_var:
        # Here we have assumed that the data_r_var_list reactions data is correct
        # This is questionable as the data may be incorrect
        print("Merging data assuming equations with a var list data are correct...", flush=True)
        df_final = pd.concat([data_r, data_r_var_list, data_r_rebalanced])
    else:
        print("Merging data assuming equations with a var list data are incorrect...", flush=True)
        df_final = pd.concat([data_r, data_r_rebalanced])

    # Get the final length of the data
    print(f"Final data shape: {df_final.shape}", flush=True)

    # Sort by the index
    df_final = df_final.sort_values(by="id")
    # Write the data to a file
    df_final.to_csv(out_eq_file, encoding='utf-8', index=False)

    # Close the log files
    f_log.close()
    f_rebalance.close()

    return df_final
