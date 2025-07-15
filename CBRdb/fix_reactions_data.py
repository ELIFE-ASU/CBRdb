import os
import time

import pandas as pd
import swifter
from rdkit import Chem as Chem

from .tools_eq import (convert_formula_to_dict,
                       side_to_dict,
                       get_elements_from_eq,
                       compare_dict_values,
                       standardise_eq,
                       check_contains_var_list,
                       check_missing_formulas,
                       full_check_eq_unbalanced,
                       rebalance_eq,
                       fix_imbalance_core,
                       )
from .tools_mols import (get_small_compounds, get_compounds_with_matching_elements, standardize_mol)


def print_and_log(statement, file=None):
    """
    Prints and logs a statement.

    Parameters:
    statement (str): The statement to be printed and logged.
    file (file object): The file object to log the statement.

    Returns:
    None
    """
    print(statement, flush=True)
    if file:
        file.write(statement + "\n")
    return None


def kitchen_sink(eq, data_c, small_compounds, f_log):
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
    print_and_log(f"Compounds that might match:  {compounds}", f_log)

    if not compounds:
        print_and_log("No compounds found that match the element differences.", f_log)
        return False

    eq_new = eq
    attempt = 1
    for compound in compounds:
        # Counter
        print_and_log(f"Attempt {attempt}, {compound}", f_log)

        # Get the elements
        _, _, react_ele, prod_ele = get_elements_from_eq(eq_new, data_c)
        # Check the difference
        diff_ele_react, diff_ele_prod = compare_dict_values(react_ele, prod_ele)
        print_and_log(f"Differences in reactants:    {diff_ele_react}", f_log)
        print_and_log(f"Differences in products:     {diff_ele_prod}", f_log)

        # Inject the compound
        eq_new = fix_imbalance_core(eq_new, diff_ele_react, diff_ele_prod, compound)
        print_and_log(f"New equation:                {eq_new}", f_log)

        # Get the elements
        _, _, react_ele, prod_ele = get_elements_from_eq(eq_new, data_c)

        # Check the difference
        diff_ele_react, diff_ele_prod = compare_dict_values(react_ele, prod_ele)
        print_and_log(f"Differences in reactants:    {diff_ele_react}", f_log)
        print_and_log(f"Differences in products:     {diff_ele_prod}", f_log)

        # Check if the equation is balanced by checking the differences
        if not diff_ele_react and not diff_ele_prod:
            print_and_log("Differences null, balanced equation found!", f_log)
            # Return the new equation line
            print_and_log(f"Final equation:              {eq_new}", f_log)
            return eq_new

        # Try to rebalance the equation
        eq_new = rebalance_eq(eq_new, data_c)
        if eq_new is False:
            # Revert to the original equation and try the next compound
            eq_new = eq
            print_and_log("Rebalance failed, trying next compound...", f_log)
        else:
            # Return the new equation line
            print_and_log(f"Final equation:              {eq_new}", f_log)
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
                       rebalance_depth=1,
                       bad_criterion='shortcut|structure_missing'):
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
    f_rebalance.write("id, reason\n")

    # Get the output file name
    out_eq_file = f"{r_file.split('.')[0]}_processed.csv".replace('_deduped', '')

    # Read the bad reactions file
    bad_ids = pd.read_csv(bad_file, index_col=0).query('reason.str.contains(@bad_criterion)').index.tolist()

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
    print_and_log("Loading the reaction data...", f_log)
    data_r = pd.read_csv(r_file)
    print_and_log("data loaded", f_log)
    print_and_log(f"data columns: {data_r.columns}", f_log)
    print_and_log(f"data shape: {data_r.shape}", f_log)

    # Sort by the index
    data_r = data_r.sort_values(by="id")

    # Filter out the bad ids, data that is shortcut/general/incomplete
    print_and_log("Filtering out bad ids", f_log)
    bool_bad_ids = data_r["id"].isin(bad_ids)
    data_r = data_r.loc[~bool_bad_ids]
    print_and_log(f"Number of bad ids removed: {sum(bool_bad_ids)}", f_log)

    # Filter out the data that has missing formulas
    print_and_log("Filtering out missing formulas", f_log)
    if f_parallel:
        t0 = time.time()
        bool_missing_data = data_r['reaction'].swifter.force_parallel(enable=True).allow_dask_on_strings(
            enable=True).apply(check_missing_formulas, args=(data_c,))
        print_and_log(f"Time to check missing formulas: {time.time() - t0}", f_log)
    else:
        t0 = time.time()
        bool_missing_data = data_r['reaction'].apply(check_missing_formulas, args=(data_c,))
        print_and_log(f"Time to check missing formulas: {time.time() - t0}", f_log)

    data_r_missing_data = data_r[bool_missing_data]
    if f_save_intermediate:
        data_r_missing_data.to_csv(f"{r_file.split('.')[0]}_missing_data.csv",
                                   encoding='utf-8',
                                   index=False)

    # Save the missing data to the rebalance file
    for id in data_r_missing_data["id"].tolist():
        f_rebalance.write(f"{id}, missing data\n")

    # Remove the missing data
    data_r = data_r[~bool_missing_data]
    print_and_log(f"Number of missing formulas removed: {sum(bool_missing_data)}", f_log)

    # Filter out the reactions that contain a var list
    print_and_log("Filtering out var list", f_log)
    if f_parallel:
        t0 = time.time()
        bool_var_list = data_r['reaction'].swifter.force_parallel(enable=True).apply(check_contains_var_list,
                                                                                     args=(data_c,))
        print_and_log(f"Time to check var list: {time.time() - t0}", f_log)
    else:
        t0 = time.time()
        bool_var_list = data_r['reaction'].apply(check_contains_var_list, args=(data_c,))
        print_and_log(f"Time to check var list: {time.time() - t0}", f_log)

    data_r_var_list = data_r[bool_var_list]
    if f_save_intermediate:
        data_r_var_list.to_csv(f"{r_file.split('.')[0]}_var_list.csv",
                               encoding='utf-8',
                               index=False)

    # Save the var list data to the rebalance file
    if not f_assume_var:
        for id in data_r_var_list["id"].tolist():
            f_rebalance.write(f"{id}, var list\n")

    # Remove the var list data
    data_r = data_r[~bool_var_list]
    print_and_log(f"Number of var list reactions removed: {sum(bool_var_list)}", f_log)

    # Filter out the data that is not balanced
    print_and_log("Filtering out unbalanced reactions", f_log)
    if f_parallel:
        t0 = time.time()
        bool_unbalanced = data_r['reaction'].swifter.force_parallel(enable=True).apply(full_check_eq_unbalanced,
                                                                                       args=(data_c,))
        print_and_log(f"Time to check if unbalanced: {time.time() - t0}", f_log)
    else:
        t0 = time.time()
        bool_unbalanced = data_r['reaction'].apply(full_check_eq_unbalanced, args=(data_c,))
        print_and_log(f"Time to check if unbalanced: {time.time() - t0}", f_log)

    # Get the data that is unbalanced
    data_r_unbalanced = data_r[bool_unbalanced]
    if f_save_intermediate:
        data_r_unbalanced.to_csv(f"{r_file.split('.')[0]}_unbalanced.csv",
                                 encoding='utf-8',
                                 index=False)
    # Get the data that is balanced
    data_r = data_r[~bool_unbalanced]
    # Determine the number of reactions that have been removed
    print_and_log(f"Number of unbalanced reactions: {sum(bool_unbalanced)}", f_log)

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

        print_and_log(f"\nProcessing {i}/{n_ids} {id}", f_log)
        print_and_log(f"Original equation line: {eq_line}", f_log)
        reactants, products, react_ele, prod_ele = get_elements_from_eq(eq_line, data_c)
        print_and_log(f"Reactants: {reactants}", f_log)
        print_and_log(f"Products:  {products}", f_log)

        # Check if the equation contains a '*' in either reactants or products
        if dict_ele_contains_star(react_ele, prod_ele):
            # We assume that the equation is correct and add it to the output
            if f_assume_star:
                print_and_log(f"Assuming eq {id} is correct", f_log)
                ids_out.append(id)
                eq_lines_out.append(eq_line)
                continue
            # We assume that the equation is incorrect and skip it as we cannot fix it
            else:
                print_and_log(f"Assume that the eq {id} is incorrect and skip it as we cannot fix it..", f_log)
                continue

        if id == 'R00263':
            print_and_log(f"Skipping {id} due to it being malformed", f_log)
            continue

        # Check if the equation is balanced
        eq_line_new = rebalance_eq(eq_line, data_c)
        # Equation is not balance...
        if eq_line_new is False:
            # Attempt injecting methods with the c1 compounds
            eq_line_new = kitchen_sink(eq_line, data_c, data_c_1, f_log)
            if eq_line_new is False:
                print_and_log(f"Could not fix the imbalance for eq {id} using c1 list!", f_log)
                if rebalance_depth > 1:
                    # Attempt with the c2 compounds
                    eq_line_new = kitchen_sink(eq_line, data_c, data_c_2, f_log)
                    if eq_line_new is False:
                        print_and_log(f"Could not fix the imbalance for eq {id} using c2 list!", f_log)
                        if rebalance_depth > 2:
                            # Attempt with the c3 compounds
                            eq_line_new = kitchen_sink(eq_line, data_c, data_c_3, f_log)
                            if eq_line_new is False:
                                print_and_log(f"Could not fix the imbalance for eq {id} using c3 list!", f_log)
                                ids_failed.append(id)
                                continue
                        else:
                            ids_failed.append(id)
                            continue
                else:
                    ids_failed.append(id)
                    continue

        print_and_log(f"rebalanced {id} with new equation line: {eq_line_new}", f_log)
        ids_out.append(id)
        eq_lines_out.append(eq_line_new)

    print_and_log(f"Number of reactions rebalanced: {len(ids_out)}", f_log)
    print_and_log(f"Number of reactions failed: {len(ids_failed)}", f_log)

    # Write the failed ids to the rebalance file
    for id in ids_failed:
        f_rebalance.write(f"{id}, failed to balance\n")

    # Update the "reaction" column in data_r_unbalanced using eq_lines_out and ids_out
    data_r_unbalanced.loc[data_r_unbalanced["id"].isin(ids_out), "reaction"] = eq_lines_out

    data_r_rebalanced = data_r_unbalanced

    print_and_log("Combining the data!", f_log)
    print_and_log(f"data_r shape: {data_r.shape}", f_log)
    print_and_log(f"data_r_var_list shape: {data_r_var_list.shape}", f_log)
    print_and_log(f"data_r_rebalanced shape: {data_r_rebalanced.shape}", f_log)

    # Combine the data
    if f_assume_var:
        # Here we have assumed that the data_r_var_list reactions data is correct
        # This is questionable as the data may be incorrect
        print_and_log("Merging data assuming equations with a var list data are correct...", f_log)
        to_concat = [data_r, data_r_var_list, data_r_rebalanced]
    else:
        print_and_log("Merging data assuming equations with a var list data are incorrect...", f_log)
        to_concat = [data_r, data_r_rebalanced]

    df_final = pd.concat(to_concat).set_index('id').drop(ids_failed, errors='ignore')
    if f_assume_var and len(data_r_var_list.index)>0:
            df_final['var_coeff'] = df_final.index.isin(data_r_var_list['id'])

    # Get the final length of the data
    print_and_log(f"Final data shape: {df_final.shape}", f_log)

    # Sort by the index
    df_final = df_final.sort_index().reset_index()
    # Write the data to a file
    df_final.to_csv(out_eq_file, encoding='utf-8', index=False)

    # Close the log files
    f_log.close()
    f_rebalance.close()

    return df_final


def compound_lookup_tables(data_c, f_log=None):

    # DataFrame of compound attributes relevant for balancing reactions
    print_and_log('making "cpd_data": DataFrame of compound attributes relevant for balancing reactions', f_log)
    cpd_data = data_c[['formula','smiles']].copy(deep=True).dropna().assign(
        formula_dict = lambda x: x.formula.map(convert_formula_to_dict),
        starred = lambda x: x.formula_dict.map(lambda y: '*' in y),
        n_elements = lambda x: x.formula_dict.map(len),
        formal_charge = lambda x: (x.smiles.map(Chem.MolFromSmiles) # returns None for 7 CBRdb_C IDs
                                   .map(standardize_mol, na_action='ignore')
                                   .map(Chem.GetFormalCharge, na_action='ignore')
                                   .fillna(0).astype(int))) # when standardized from file, all 7 have charge=0

    # DataTable indicating, for each compound (column), the count (value) of each element (row)
    print_and_log('making "formula_table": matrix-like DataFrame of element counts for each compound', f_log)
    formula_table = cpd_data['formula_dict'].apply(pd.Series).fillna(0).astype(int).T
    
    return cpd_data, formula_table


def filter_reactions_pandas(data_r, data_c, formula_table=None, f_log=None):
    
    # DataFrame of compound attributes relevant for balancing reactions
    if set(['formula_dict', 'starred', 'formal_charge']).issubset(set(data_c.columns)) and formula_table is not None:
        cpd_data = data_c
    else:
        cpd_data, formula_table = compound_lookup_tables(data_c, f_log=f_log)

    # DataFrame of dicts of reactants and products
    print_and_log('making "sides": DataFrame of reactant and product dicts', f_log)
    sides = data_r['reaction'].str.split('<=>', expand=True).map(side_to_dict)
    # Function to get a GroupBy object returning whether reaction participants are in a listlike object
    cpd_group_attrs = lambda lstlike: sides.map(lambda x: x.keys()).stack().explode().isin(lstlike).groupby(level=0)

    # DataFrame of reaction attributes
    print_and_log('making "rns": DataFrame of reaction attributes e.g. missing data, starred compounds, vars-as-coefficients', f_log)
    rn_attrs = pd.DataFrame({'bool_missing_data': ~ cpd_group_attrs(cpd_data.index).all(),
                            'cpd_starred': cpd_group_attrs(cpd_data.query('starred').index).any(),
                            'bool_var_list': ~ sides.map(lambda x: all([str(i).isnumeric() for i in x.values()])).all(axis=1)})
    rn_attrs['rebalanceable'] = ~ rn_attrs.any(axis=1)

    # DataFrame for each side, indicating the numeric coefficient (value) for each compound (row) in each reaction (column)
    print_and_log('making "reactant_cps", "product_cps": shows L+R sides as matrix-like DataFrame of coefficients for each compound in each reaction', f_log)
    reactant_cps, product_cps = [i.drop(rn_attrs.query('bool_missing_data or bool_var_list').index)
                                 .apply(pd.Series).T for _, i in sides.items()]

    # calculate stoichiometry of each reaction as currently written
    print_and_log('making "reactant_els", "product_els": matrix-like DataFrames of current element counts for each reaction', f_log)
    reactant_els, product_els = dict(), dict()
    for id, reactants in reactant_cps.items():  # for each reaction
        coeffs = reactants.dropna()  # access reactant IDs and coefficients
        reactant_els[id] = (coeffs * formula_table[coeffs.index]).sum(axis=1)  # sum elements
    for id, products in product_cps.items():
        coeffs = products.dropna()
        product_els[id] = (coeffs * formula_table[coeffs.index]).sum(axis=1)
    reactant_els = pd.DataFrame(reactant_els).T.astype(int)
    product_els = pd.DataFrame(product_els).T.astype(int)

    # ID whether reaction is balanced, balanceable, or (to inform treatment of starred reactions) balanced except for *
    print_and_log('using "reactant_els and "product_els" to add columns to DataFrame "rns":', f_log)
    print_and_log('  * is_balanced: whether the reaction is balanced as written', f_log)
    rn_attrs['is_balanced'] = (product_els == reactant_els).all(axis=1)
    print_and_log('  * to_rebalance: unbalanced reactions that lack missing data, starred compounds, or vars-as-coefficients', f_log)
    rn_attrs['to_rebalance'] = rn_attrs['is_balanced'].eq(False) & rn_attrs['rebalanceable'].eq(True)
    print_and_log('  * is_balanced_except_star: reactions for which "*" is the only imbalanced element', f_log)
    rn_attrs['is_balanced_except_star'] = (product_els == reactant_els).drop(columns='*', errors='ignore').all(axis=1) & \
                                          rn_attrs['is_balanced'].eq(False)
    print_and_log('  * charge_R-L: charge (im)balance across the reaction (charge_R - charge_L). If positive, right needs (-) or left needs (+)', f_log)
    rn_attrs['charge_R-L'] = ((product_cps.T * cpd_data['formal_charge'].loc[product_cps.index]).sum(axis=1) -
                            (reactant_cps.T * cpd_data['formal_charge'].loc[reactant_cps.index]).sum(axis=1)).astype(int)

    # pd.Series listing sets of formulas for each side of the reaction; can use Series.apply chempy.balance_stoichiometry to check in bulk.
    print_and_log('making "formula_sides": pd.Series listing sets of formulas for to-rebalance reactions above.', f_log)
    formula_sides = pd.Series(
        sides.loc[rn_attrs['to_rebalance']].map(lambda x: set(data_c.formula.loc[x.keys()])).T.to_dict(orient='list'))

    # DataFrame of element count diffs (R - L); same count diff might indicate same set of compound-injector solutions, saving iterations
    print_and_log('making "el_diff_groups": DataFrame assigning group numbers to reactions based on their element-count diff (R - L).', f_log)
    observed_el_diffs = (product_els - reactant_els).loc[rn_attrs.query('to_rebalance').index]
    observed_el_diffs['group_num'] = (observed_el_diffs.astype(str)+' ').groupby(by=list(observed_el_diffs.columns)).ngroup()
    observed_el_diffs['group_size'] = observed_el_diffs['group_num'].map(observed_el_diffs['group_num'].value_counts())
    observed_el_diffs = observed_el_diffs.sort_values(by=['group_size', 'group_num'], ascending=[False, True]).drop(
        'group_size', axis=1).set_index('group_num', append=True)

    print_and_log('Combining all dataframes into a dict.', f_log)
    dfs = {'reactant_cps': reactant_cps, 'product_cps': product_cps,
           'reactant_els': reactant_els, 'product_els': product_els,
           'rns': rn_attrs, 'formula_table': formula_table,
           'sides': sides, 'formula_sides': formula_sides,
           'cpd_data': cpd_data, 
           'el_diff_groups': observed_el_diffs}

    print_and_log(f'Entries: {list(dfs.keys())}', f_log)
    return dfs


def get_charge_balanced_injections_1el(data_c, dfs):

    cpd_pool = data_c.query('n_elements==1 & starred==False').copy(deep=True).assign(
        el_sym = lambda x: x.formula_dict.explode(),
        el_num = lambda x: x.formula_dict.map(lambda y: y.values()).explode()).reset_index()
    cpd_pool = cpd_pool.groupby(by=['el_sym', 'formal_charge', 'el_num'], as_index=False)['compound_id'].agg(lambda x: x)

    diffs = dfs['el_diff_groups'][dfs['el_diff_groups']['*'].eq(0)]     # consider only reactions with no starred compounds
    sum_all_atoms = diffs.abs().sum(axis=1)
    single_el_diffs = (diffs[(diffs.abs().T == sum_all_atoms).any()] # el diff = total diff necessarily means a one-element mismatch
                    .reset_index(level=0).melt(id_vars='id', value_name='diff', var_name='el_sym') # so only one diff per reaction
                    .query('diff!=0').set_index('id').assign(formal_charge=dfs['rns']['charge_R-L'])) # diffs are RE: element and charge.

    cpd_pool = cpd_pool.query('el_sym in @single_el_diffs.el_sym')
    single_el_diffs = single_el_diffs.assign(compound_id=None, left_side=None, count=None)
    queries = {'formal_charge == 0': '(formal_charge == 0)',
            'formal_charge == 1': '(formal_charge == diff)',
            'formal_charge == -1': '(formal_charge == -1 * diff)'}
    
    for charge_specs, solution_query in queries.items():
        for _, soln in cpd_pool.query(charge_specs).iterrows():
            solvable = single_el_diffs.query(f'(el_sym == @soln.el_sym) & (diff.abs() % @soln.el_num == 0) & {solution_query}')
            if len(solvable.index) > 0:
                single_el_diffs.loc[solvable.index, 'compound_id'] = soln['compound_id']
                single_el_diffs.loc[solvable.index, 'count'] = (solvable['diff'].abs() / soln['el_num']).astype(int)
                single_el_diffs.loc[solvable.index, 'left_side'] = solvable['diff'] > 0

    return single_el_diffs.dropna()


def get_charge_balanced_injections_OH(data_c, dfs):
    solutions = pd.DataFrame(columns=['el_sym', 'diff', 'formal_charge', 'compound_id', 'left_side', 'count'])
    diffs = dfs['el_diff_groups'][dfs['el_diff_groups']['*'].eq(0)]
    sum_all_atoms = diffs.abs().sum(axis=1)
    diffs = diffs.loc[diffs.query('H.ne(0) & O.ne(0)')[['H','O']].sum(axis=1).abs().eq(sum_all_atoms)][['H','O']].reset_index(level=1, drop=True)
    diffs['formal_charge'] = dfs['rns']['charge_R-L'].loc[diffs.index.get_level_values(0)].values
    H2O_solves = diffs.query('formal_charge.eq(0) & H == 2*O')
    H2O2_solves = diffs.query('H==O & formal_charge.eq(0) & H%2==0')
    solutions = pd.DataFrame()
    if len(H2O_solves.index)>0:
        H2O_solves = H2O_solves.assign( el_sym = 'H2O', diff = lambda x: x.O, compound_id = 'C00001', 
                                       left_side = lambda x: x.H.gt(0), count = lambda x: x.O.abs()).drop(['H','O'], axis=1)
        solutions = pd.concat([solutions, H2O_solves])
    if len(H2O2_solves.index)>0:
        H2O2_solves = H2O2_solves.assign(el_sym = 'H2O2', diff = lambda x: (x.O/2).astype(int), compound_id = 'C00027', 
                                         left_side= lambda x: x['diff'].gt(0), count = lambda x: x['diff'].abs()).drop(['H','O'], axis=1)
        solutions = pd.concat([solutions, H2O2_solves])
    return solutions
