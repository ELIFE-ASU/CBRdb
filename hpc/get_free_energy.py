import numpy as np
from equilibrator_api import ComponentContribution
import CBRdb
import pandas as pd
import warnings
from functools import partial


# warnings.filterwarnings("error")


def get_dg_prime(eq_line, conditions='standard'):
    cc = ComponentContribution()
    # replace 'C' with 'kegg:C' to use KEGG IDs
    eq_line = eq_line.replace('C', 'kegg:C')
    # replace '<=>' with '=' for equilibrium
    eq_line = eq_line.replace('<=>', '=')

    try:
        react_line = cc.parse_reaction_formula(eq_line)
        rev = cc.ln_reversibility_index(react_line).value.m_as("dimensionless")
        if conditions == 'standard':
            dg = cc.standard_dg_prime(react_line)
        else:
            # assuming physiological conditions
            dg = cc.physiological_dg_prime(react_line)
        val = dg.value.m_as("kJ/mol")
        error = dg.error.m_as("kJ/mol")
    except Exception as e:
        print(f"Error: {e}")
        return np.nan, np.nan, np.nan

    if error == 100000.0:
        error = np.nan
    return val, error, rev


def stderr_from_sigmas(sigma_fin, sigma_inf, rmse_inf=1e-5):
    term = np.sum(np.asarray(sigma_fin, float) ** 2)
    sigma_inf = np.asarray(sigma_inf[0], float)
    term += np.sum((rmse_inf * sigma_inf) ** 2)
    return np.sqrt(term)


def get_energy_of_formation(in_file='../CBRdb_C.csv',
                            out_file='CBRdb_C_formation_energies.csv',
                            run_physiological_transform=True):
    data_c = pd.read_csv(in_file, low_memory=False)

    # Make a new dataframe with only the compound_id and name columns
    data_c = pd.DataFrame(data_c['compound_id'])
    data_c.head(100)
    # Make sure there are no duplicates or NaNs
    data_c = data_c.drop_duplicates()
    data_c = data_c.dropna()
    data_c = data_c.reset_index(drop=True)

    # Get the all compound ID lists and add KEGG-formatted IDs as a new column.
    data_c['kegg_id'] = data_c['compound_id'].apply(lambda cid: f'kegg:{cid}')

    cc = ComponentContribution()
    data_c['cc_comp'] = data_c['kegg_id'].apply(lambda cid: cc.get_compound(cid))
    # Filter to only those compounds that are in the ComponentContribution database
    data_c = data_c.dropna(subset=['cc_comp'])
    data_c = data_c.reset_index(drop=True)

    # Extract lists from dataframe and calculate standard Gibbs free energy and sigmas
    data_c[['standard_dgf_mu', 'sigma_fin', 'sigma_inf']] = pd.DataFrame(
        map(cc.standard_dg_formation, data_c['cc_comp'].tolist())
    )
    # Drop None values
    data_c = data_c.dropna(subset=['standard_dgf_mu'])
    data_c = data_c.reset_index(drop=True)

    # Apply stderr_from_sigmas to each row to get the standard error
    data_c['standard_error'] = data_c.apply(lambda row: stderr_from_sigmas(row['sigma_fin'], row['sigma_inf']), axis=1)

    # Loop over the first 10 and print the compound_id, standard_dgf_mu, and standard_error
    for i in range(10):
        print(
            f"{data_c['compound_id'][i]}: {data_c['standard_dgf_mu'][i]:.2f} ± {data_c['standard_error'][i]:.2f} kJ/mol")

    if run_physiological_transform:
        # Apply the Legendre transform to convert from the standard ΔGf to the standard ΔG'f
        data_c['delta_dgf'] = data_c['cc_comp'].apply(
            lambda cpd: cpd.transform(cc.p_h, cc.ionic_strength, cc.temperature, cc.p_mg).m_as("kJ/mol")
        )
        data_c['standard_dgf_prime_mu'] = data_c['standard_dgf_mu'] + data_c['delta_dgf']
        # calculate the covariance matrix
        sigmas_fin = np.array(data_c['sigma_fin'].tolist()).T
        sigmas_inf = np.array(data_c['sigma_inf'].tolist()).T
        standard_dgf_cov = sigmas_fin @ sigmas_fin.T + 1e6 * sigmas_inf @ sigmas_inf.T
        # convert the covariance matrix to a single error value (1-σ) for each compound
        data_c['standard_dgf_prime_error'] = np.sqrt(np.diag(standard_dgf_cov))

        # loop over the first 10 and print the compound_id, standard_dgf_prime_mu, and standard_dgf_prime_error
        for i in range(10):
            print(
                f"{data_c['compound_id'][i]}: {data_c['standard_dgf_prime_mu'][i]:.2f} ± {data_c['standard_dgf_prime_error'][i]:.2f} kJ/mol")

    # Save the results to a CSV file
    data_c.to_csv(out_file, index=False)
    return data_c


if __name__ == "__main__":
    print(flush=True)
    get_energy_of_formation()

    exit()

    # https://equilibrator.readthedocs.io/en/latest/equilibrator_examples.html#Using-formation-energies-to-calculate-reaction-energies

    data = pd.read_csv('../CBRdb_R.csv', low_memory=False)
    # get the reaction equations
    ids = data['id'].tolist()
    eq_lines = data['reaction'].tolist()[:1000]
    eq_lines = [l.replace('C', 'kegg:C') for l in eq_lines]
    eq_lines = [l.replace('<=>', '=') for l in eq_lines]
    cc = ComponentContribution()
    tmp = cc.get_compound('kegg:C98976')  # kegg:C98976
    print(tmp, flush=True)
    exit()

    reactions = [cc.parse_reaction_formula(l) for l in eq_lines]

    dg, uc = cc.standard_dg_prime_multi(reactions, uncertainty_representation="cov")
    print(uc, flush=True)
    # get the diagonal of the covariance matrix
    uc = np.sqrt(np.diag(uc.magnitude))

    print(dg, flush=True)
    print(uc, flush=True)

    # broken 1.41421356e+05

    # for i in range(len(eq_lines)):
    #     print(get_dg_prime(eq_lines[i]))

    # tmp = partial(get_dg_prime, conditions='physiological')
    # # use multiprocessing to speed up the calculation
    # results = CBRdb.mp_calc(tmp, eq_lines)
