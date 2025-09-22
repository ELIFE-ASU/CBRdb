from functools import partial

import numpy as np
import pandas as pd
from equilibrator_api import ComponentContribution, Q_

import CBRdb


def err_from_sig(sigma_fin, sigma_inf, rmse_inf=1e-5):
    term = np.sum(np.asarray(sigma_fin, float) ** 2)
    sigma_inf = np.asarray(sigma_inf[0], float)
    term += np.sum((rmse_inf * sigma_inf) ** 2)
    return np.sqrt(term)


def get_working_cids(data_c):
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
    try:
        return cc.ln_reversibility_index(r).value.m_as("dimensionless")
    except Exception as e:
        print(f"Error calculating reversibility index for reaction {r}: {e}", flush=True)
        return None


def get_std_dg(cc, data_r):
    try:
        dg, uc = cc.standard_dg_prime_multi(data_r['cc_reac'].tolist(), uncertainty_representation="cov")
        return dg.magnitude, np.sqrt(np.diag(uc.magnitude))
    except Exception as e:
        print(f"Error calculating standard Gibbs free energy for reactions: {e}", flush=True)
        return None, None


def prepare_reaction_data(data_r, data_c):
    # Trim them down to only the necessary columns
    data_r = pd.DataFrame(data_r[['id', 'reaction']])

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
    # count the number of 'C' in each reaction and drop rows with less than 2
    data_r = data_r[data_r['reaction'].str.count('C') >= 2]
    data_r = data_r.reset_index(drop=True)

    print(f'After removing unworkable or incomplete. {len(data_r)} reactions remain.', flush=True)

    # Fix the formatting
    data_r['kegg_reaction'] = data_r['reaction'].apply(
        lambda eq: eq.replace('C', 'kegg:C').replace('<=>', '=')
    )
    return data_r


def _get_energy_of_reaction(data_r, run_physiological=True):
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
    print('Reading CBRdb reaction data...', flush=True)

    # Load the compound and reaction data
    data_c = pd.read_csv(in_file_c, low_memory=False)
    data_r = pd.read_csv(in_file_r, low_memory=False)

    data_r_full = prepare_reaction_data(data_r, data_c)
    total_reactions = len(data_r_full)
    total_segments = 10
    segment_size = total_reactions // total_segments
    print(segment_size)

    # Process each segment sequentially otherwise we run out of memory
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

    # Join the data
    data_r = pd.concat(all_results, ignore_index=True)

    if compact_output:
        # Keep only relevant columns
        cols_to_keep = ['id', 'std_dg', 'std_dg_error', 'rev_index']
        if run_physiological:
            cols_to_keep += ['std_dg_p', 'std_dg_p_error', 'rev_index_p']
        data_r = data_r[cols_to_keep]

    # Save the results to a CSV file
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
