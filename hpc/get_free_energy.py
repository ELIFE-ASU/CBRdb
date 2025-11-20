from functools import partial

import numpy as np
import pandas as pd
from equilibrator_api import ComponentContribution, Q_

import CBRdb


def err_from_sig(sigma_fin, sigma_inf, rmse_inf=1e-5):
    """
    Calculate the error from sigma values.

    This function computes the root mean square error (RMSE) based on the provided
    sigma values and a specified RMSE for infinite sigma.

    Parameters:
    -----------
    sigma_fin : array-like
        The finite sigma values as an array or list.
    sigma_inf : array-like
        The infinite sigma values as an array or list. Only the first element is used.
    rmse_inf : float, optional
        The RMSE value for infinite sigma. Default is 1e-5.

    Returns:
    --------
    float
        The calculated root mean square error.
    """
    term = np.sum(np.asarray(sigma_fin, float) ** 2)
    sigma_inf = np.asarray(sigma_inf[0], float)
    term += np.sum((rmse_inf * sigma_inf) ** 2)
    return np.sqrt(term)


def get_working_cids(data_c):
    """
    Process a DataFrame of compounds to retrieve valid ComponentContribution objects.

    This function adds KEGG-formatted IDs to the input DataFrame, retrieves the corresponding
    ComponentContribution objects, and filters out compounds that are not present in the
    ComponentContribution database.

    Parameters:
    -----------
    data_c : pandas.DataFrame
        A DataFrame containing a 'compound_id' column with compound identifiers.

    Returns:
    --------
    pandas.DataFrame
        A filtered DataFrame with an additional 'kegg_id' column (KEGG-formatted IDs)
        and 'cc_comp' column (ComponentContribution objects).
        Only compounds present in the ComponentContribution database are retained.
    """
    # Get the all compound ID lists and add KEGG-formatted IDs as a new column.
    data_c['kegg_id'] = data_c['compound_id'].apply(lambda cid: f'kegg:{cid}')

    cc = ComponentContribution()
    data_c['cc_comp'] = data_c['kegg_id'].apply(lambda cid: cc.get_compound(cid))
    # Filter to only those compounds that are in the ComponentContribution database
    data_c = data_c.dropna(subset=['cc_comp'])
    data_c = data_c.reset_index(drop=True)
    return data_c


def get_energy_of_formation(in_file='../CBRdb_C.csv',
                            out_file='CBRdb_C_formation_energies.csv.gz',
                            compact_output=True,
                            run_physiological=True):
    """
    Calculate and save the Gibbs free energy of formation for compounds.

    This function reads compound data from a specified input file, calculates the
    standard Gibbs free energy of formation using the `equilibrator-api`, and
    optionally calculates it under physiological conditions. The results, including
    energies and errors, are saved to a gzipped CSV file.

    Parameters
    ----------
    in_file : str, optional
        Path to the input CSV file containing compound data (e.g., from CBRdb).
        Default is '../CBRdb_C.csv'.
    out_file : str, optional
        Path to the output gzipped CSV file where formation energies will be saved.
        Default is 'CBRdb_C_formation_energies.csv.gz'.
    compact_output : bool, optional
        If True, the output file will contain only the essential columns
        (ID, energy, and error). Default is True.
    run_physiological : bool, optional
        If True, calculations will also be performed under physiological conditions
        (pH 7.5, pMg 3.0) in addition to standard conditions. Default is True.

    Returns
    -------
    None
    """
    print('Reading CBRdb compound data...', flush=True)
    data_c = pd.read_csv(in_file, low_memory=False)

    # Make a new dataframe with only the compound_id and name columns
    data_c = pd.DataFrame(data_c['compound_id'])
    # Get the working compound list
    data_c = get_working_cids(data_c)
    print(f'Calculating formation energies for {len(data_c)} compounds...', flush=True)
    cc = ComponentContribution()
    cc.p_h = Q_(7.0)
    cc.p_mg = Q_(10.0)
    cc.ionic_strength = Q_(0.25, "M")
    cc.temperature = Q_(298.15, "K")
    # Extract lists from dataframe and calculate standard Gibbs free energy and sigmas
    data_c[['std_dgf', 'sigma_fin', 'sigma_inf']] = pd.DataFrame(
        map(cc.standard_dg_formation, data_c['cc_comp'].tolist())
    )
    # Drop None values
    data_c = data_c.dropna(subset=['std_dgf'])
    data_c = data_c.reset_index(drop=True)

    # Calculate the standard error from the sigmas
    data_c['std_dgf_error'] = data_c.apply(lambda row: err_from_sig(row['sigma_fin'], row['sigma_inf']), axis=1)
    print('First 10 formation energies:', flush=True)
    for i in range(10):
        print(
            f"{data_c['compound_id'][i]}: {data_c['std_dgf'][i]:.2f} ± {data_c['std_dgf_error'][i]:.2f} kJ/mol")

    if run_physiological:
        print('Calculating formation energies under physiological conditions...', flush=True)
        cc = ComponentContribution()
        cc.p_h = Q_(7.5)
        cc.p_mg = Q_(3.0)
        cc.ionic_strength = Q_(0.25, "M")
        cc.temperature = Q_(298.15, "K")
        # Extract lists from dataframe and calculate standard Gibbs free energy and sigmas
        data_c[['std_dgf_p', 'sigma_fin_p', 'sigma_inf_p']] = pd.DataFrame(
            map(cc.standard_dg_formation, data_c['cc_comp'].tolist())
        )
        # Drop None values
        data_c = data_c.dropna(subset=['std_dgf_p'])
        data_c = data_c.reset_index(drop=True)
        # Calculate the standard error from the sigmas
        data_c['std_dgf_p_error'] = data_c.apply(lambda row: err_from_sig(row['sigma_fin_p'], row['sigma_inf_p']),
                                                 axis=1)
        print('First 10 formation energies:', flush=True)
        for i in range(10):
            print(
                f"{data_c['compound_id'][i]}: {data_c['std_dgf_p'][i]:.2f} ± {data_c['std_dgf_p_error'][i]:.2f} kJ/mol")

    if compact_output:
        # Keep only relevant columns
        cols_to_keep = ['compound_id', 'std_dgf', 'std_dgf_error']
        if run_physiological:
            cols_to_keep += ['std_dgf_p', 'std_dgf_p_error']
        data_c = data_c[cols_to_keep]

    # Save the results to a CSV file
    print(f'Saving results to {out_file}', flush=True)
    print(f'Calculated formation energies for {len(data_c)} compounds.', flush=True)
    data_c.to_csv(out_file, index=False, compression='gzip')
    print('Done.', flush=True)
    print(flush=True)
    return None


def get_rev_index(cc, r):
    """
    Calculate the natural logarithm of the reversibility index for a reaction.

    This function computes the natural logarithm of the reversibility index for a given
    reaction using the ComponentContribution model. If an error occurs during the calculation,
    it logs the error and returns None.

    Parameters
    ----------
    cc : ComponentContribution
        An instance of the ComponentContribution class used for the calculation.
    r : equilibrator_api.Reaction
        The reaction object for which the reversibility index is calculated.

    Returns
    -------
    float or None
        The natural logarithm of the reversibility index as a dimensionless value, or None
        if an error occurs during the calculation.
    """
    try:
        return cc.ln_reversibility_index(r).value.m_as("dimensionless")
    except Exception as e:
        print(f"Error calculating reversibility index for reaction {r}: {e}", flush=True)
        return None


def get_std_dg(cc, data_r):
    """
    Calculate the standard Gibbs free energy and its uncertainty for reactions.

    This function computes the standard Gibbs free energy (ΔG°') and its associated
    uncertainty for a set of reactions using the ComponentContribution model.

    Parameters
    ----------
    cc : ComponentContribution
        An instance of the ComponentContribution class used for the calculation.
    data_r : pandas.DataFrame
        A DataFrame containing a 'cc_reac' column with parsed reaction objects.

    Returns
    -------
    tuple
        A tuple containing:
        - dg (float): The calculated standard Gibbs free energy (ΔG°') in kJ/mol.
        - uncertainty (float): The uncertainty of the calculated ΔG°' in kJ/mol.
        If an error occurs, both values will be None.
    """
    try:
        dg, uc = cc.standard_dg_prime_multi(data_r['cc_reac'].tolist(), uncertainty_representation="cov")
        return dg.magnitude, np.sqrt(np.diag(uc.magnitude))
    except Exception as e:
        print(f"Error calculating standard Gibbs free energy for reactions: {e}", flush=True)
        return None, None


def prepare_reaction_data(data_r, data_c):
    """
    Prepare and clean reaction data for further processing.

    This function processes a DataFrame of reactions by:
    - Filtering out reactions with compounds not in the working compound list.
    - Removing reactions with unknown stoichiometry or invalid formatting.
    - Ensuring reactions contain at least two compounds.
    - Formatting reaction equations to KEGG-compatible format.

    Parameters
    ----------
    data_r : pandas.DataFrame
        A DataFrame containing reaction data with at least 'id' and 'reaction' columns.
    data_c : pandas.DataFrame
        A DataFrame containing compound data with a 'compound_id' column.

    Returns
    -------
    pandas.DataFrame
        A cleaned and formatted DataFrame of reactions with an additional 'kegg_reaction' column.
    """
    # Trim them down to only the necessary columns
    data_r = pd.DataFrame(data_r[['id', 'reaction']])

    # Create a set of working compound IDs
    working_compounds = set(data_c['compound_id'].tolist())
    print(f'We have {len(working_compounds)} working compounds.', flush=True)

    # Get the compound IDs for each reaction if not already present
    if 'CBRdb_C_ids' not in data_r.columns:
        data_r['CBRdb_C_ids'] = data_r['reaction'].apply(lambda i: CBRdb.get_eq_all_cids(i))

    # Prune the reactions that have compounds not in the working compound list
    def is_reaction_prunable(cids):
        return any(cid not in working_compounds for cid in cids)

    print(f'Initial number of reactions: {len(data_r)}', flush=True)
    data_r['prunable'] = data_r['CBRdb_C_ids'].apply(is_reaction_prunable)
    print(f'Pruning {data_r["prunable"].sum()} reactions that contain unknown compounds.', flush=True)
    data_r = data_r[~data_r['prunable']]
    data_r = data_r.reset_index(drop=True)
    print(f'After pruning, {len(data_r)} reactions remain.', flush=True)

    # Drop rows where the equation has an n, m or x (unknown stoichiometry)
    data_r = data_r[~data_r['reaction'].str.contains(r'[nmx]')]
    data_r = data_r.reset_index(drop=True)
    # Drop rows where the equation is NaN or empty
    data_r = data_r.dropna(subset=['reaction'])
    # Drop rows that do not contain a valid reaction arrow
    data_r = data_r[data_r['reaction'].str.contains(r'<=>')]
    data_r = data_r.reset_index(drop=True)
    # Count the number of 'C' in each reaction and drop rows with less than 2
    data_r = data_r[data_r['reaction'].str.count('C') >= 2]
    data_r = data_r.reset_index(drop=True)

    print(f'After removing unworkable or incomplete. {len(data_r)} reactions remain.', flush=True)

    # Fix the formatting
    data_r['kegg_reaction'] = data_r['reaction'].apply(
        lambda eq: eq.replace('C', 'kegg:C').replace('<=>', '=')
    )
    return data_r


def _get_energy_of_reaction(data_r, run_physiological=True):
    """
    Calculate the Gibbs free energy and reversibility index for reactions.

    This function computes the standard Gibbs free energy (ΔG°') and reversibility index
    for a set of reactions under standard conditions. Optionally, it also calculates these
    values under physiological conditions.

    Parameters
    ----------
    data_r : pandas.DataFrame
        A DataFrame containing reaction data with a 'kegg_reaction' column.
    run_physiological : bool, optional
        If True, the calculations are also performed under physiological conditions
        (pH 7.5, pMg 3.0). Default is True.

    Returns
    -------
    pandas.DataFrame
        The input DataFrame with additional columns:
        - 'cc_reac': Parsed reaction objects.
        - 'rev_index': Reversibility index under standard conditions.
        - 'std_dg': Standard Gibbs free energy under standard conditions.
        - 'std_dg_error': Uncertainty of the standard Gibbs free energy.
        - 'rev_index_p' (if run_physiological=True): Reversibility index under physiological conditions.
        - 'std_dg_p' (if run_physiological=True): Gibbs free energy under physiological conditions.
        - 'std_dg_p_error' (if run_physiological=True): Uncertainty of the Gibbs free energy under physiological conditions.
    """
    # Apply to get the eq cc format
    cc = ComponentContribution()
    cc.p_h = Q_(7.0)
    cc.p_mg = Q_(10.0)
    cc.ionic_strength = Q_(0.25, "M")
    cc.temperature = Q_(298.15, "K")
    print('Calculating reversibility for standard conditions...', flush=True)
    func_rev_idx = partial(get_rev_index, cc)
    data_r['cc_reac'] = data_r['kegg_reaction'].apply(lambda r: cc.parse_reaction_formula(r))
    data_r['rev_index'] = data_r['cc_reac'].apply(lambda r: func_rev_idx(r))

    print(f'Calculating reaction energies for {len(data_r)} reactions...', flush=True)
    data_r['std_dg'], data_r['std_dg_error'] = get_std_dg(cc, data_r)

    if run_physiological:
        print('Calculating reaction energies under physiological conditions...', flush=True)
        cc = ComponentContribution()
        cc.p_h = Q_(7.5)
        cc.p_mg = Q_(3.0)
        cc.ionic_strength = Q_(0.25, "M")
        cc.temperature = Q_(298.15, "K")
        print('Calculating reversibility indices under physiological conditions...', flush=True)
        func_rev_idx = partial(get_rev_index, cc)
        data_r['cc_reac'] = data_r['kegg_reaction'].apply(lambda r: cc.parse_reaction_formula(r))
        data_r['rev_index_p'] = data_r['cc_reac'].apply(lambda r: func_rev_idx(r))

        print(f'Calculating reaction energies for {len(data_r)} reactions...', flush=True)
        data_r['std_dg_p'], data_r['std_dg_p_error'] = get_std_dg(cc, data_r)
    return data_r


def get_energy_of_reaction(in_file_c='CBRdb_C_formation_energies.csv.gz',
                           in_file_r='../CBRdb_R.csv',
                           out_file='CBRdb_R_reaction_energies.csv.gz',
                           compact_output=True,
                           run_physiological=True):
    """
    Calculate and save the Gibbs free energy and reversibility index for reactions.

    This function reads compound and reaction data from specified input files, processes
    the reactions in segments to avoid memory issues, and calculates the Gibbs free energy
    and reversibility index for each reaction. The results are saved to a gzipped CSV file.

    Parameters
    ----------
    in_file_c : str, optional
        Path to the input gzipped CSV file containing compound formation energies.
        Default is 'CBRdb_C_formation_energies.csv.gz'.
    in_file_r : str, optional
        Path to the input CSV file containing reaction data.
        Default is '../CBRdb_R.csv'.
    out_file : str, optional
        Path to the output gzipped CSV file where reaction energies will be saved.
        Default is 'CBRdb_R_reaction_energies.csv.gz'.
    compact_output : bool, optional
        If True, the output file will contain only the essential columns
        (ID, energy, error, and reversibility index). Default is True.
    run_physiological : bool, optional
        If True, calculations will also be performed under physiological conditions
        (pH 7.5, pMg 3.0) in addition to standard conditions. Default is True.

    Returns
    -------
    None
    """
    print('Reading CBRdb reaction data...', flush=True)

    # Load the compound and reaction data
    data_c = pd.read_csv(in_file_c, low_memory=False)
    data_r = pd.read_csv(in_file_r, low_memory=False)

    # Prepare and clean the reaction data
    data_r_full = prepare_reaction_data(data_r, data_c)
    total_reactions = len(data_r_full)
    total_segments = 10  # Number of segments to divide the data into
    segment_size = total_reactions // total_segments  # Size of each segment
    print(segment_size)

    # Process each segment sequentially to avoid memory issues
    all_results = []
    for i in range(total_segments):
        print(f"Processing segment {i + 1}/{total_segments}", flush=True)
        # Select the segment
        start_idx = i * segment_size
        end_idx = (i + 1) * segment_size if i < total_segments - 1 else total_reactions
        seg_r = data_r_full.iloc[start_idx:end_idx].reset_index(drop=True)

        # Process the segment
        segment_results = _get_energy_of_reaction(seg_r, run_physiological=run_physiological)
        all_results.append(segment_results)
        print(flush=True)

    # Combine the results from all segments
    data_r = pd.concat(all_results, ignore_index=True)

    if compact_output:
        # Keep only relevant columns
        cols_to_keep = ['id', 'std_dg', 'std_dg_error', 'rev_index']
        if run_physiological:
            cols_to_keep += ['std_dg_p', 'std_dg_p_error', 'rev_index_p']
        data_r = data_r[cols_to_keep]

    # Save the results to a gzipped CSV file
    print(f'Saving results to {out_file}', flush=True)
    print(f'Calculated reaction energies for {len(data_r)} reactions.', flush=True)
    data_r.to_csv(out_file, index=False, compression='gzip')
    print('Done.', flush=True)
    print(flush=True)
    return None


if __name__ == "__main__":
    print(flush=True)
    get_energy_of_formation()
    get_energy_of_reaction()
