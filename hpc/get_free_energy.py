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


# Combine finite + (optional) infinite-uncertainty parts into a scalar 1-σ (kJ/mol)
def stderr_from_sigmas(sigma_fin, sigma_inf):
    # finite part
    term = np.sum(np.asarray(sigma_fin, float)**2)
    # some versions also return sigma_inf; scale by rmse_inf (default 1e-5 kJ/mol)
    sigma_inf = np.asarray(sigma_inf[0], float)
    rmse_inf = 1e-5  # same default as the constructor
    term += np.sum((rmse_inf * sigma_inf)**2)
    return np.sqrt(term)




if __name__ == "__main__":
    print(flush=True)
    data_c = pd.read_csv('../CBRdb_C.csv', low_memory=False)
    print(data_c)
    compound_ids = data_c['compound_id'].tolist()[:2]
    print(f"Total compounds: {len(compound_ids)}", flush=True)
    compound_ids = [f'kegg:{cid}' for cid in compound_ids]

    # Get the all compound ID lists.

    # Check we can poll the compound database
    cc = ComponentContribution()
    compound_list = [cc.get_compound(i) for i in compound_ids]
    # Get the index of the bad compounds
    bad_idx = [i for i in range(len(compound_list)) if compound_list[i] is None]
    bad_compounds = [compound_ids[i] for i in bad_idx]
    print(f"Bad compounds: {len(bad_compounds)}", flush=True)

    # Remove the bad compounds from the list
    compound_list = [compound_list[i] for i in range(len(compound_list)) if i not in bad_idx]
    compound_ids = [compound_ids[i] for i in range(len(compound_ids)) if i not in bad_idx]
    print(f"Good compounds: {len(compound_list)}", flush=True)

    # # Store the bad compound IDs
    # with open('bad_compounds.txt', 'w') as f:
    #     for cid in bad_compounds:
    #         f.write(f"{cid}\n")

    # Get the standard dg formation for each compound
    print("Calculating standard formation energies...", flush=True)
    standard_dgf_mu, sigmas_fin, sigmas_inf = zip(*map(cc.standard_dg_formation, compound_list))

    # Loop over the items and remove any that have NaN values


    

    standard_dgf_mu = np.array(standard_dgf_mu)
    sigmas_fin = np.array(sigmas_fin)
    sigmas_inf = np.array(sigmas_inf)

    for i in range(len(standard_dgf_mu)):
        print(compound_ids[i], standard_dgf_mu[i])

    err_atp = stderr_from_sigmas(sigmas_fin[0], sigmas_inf[0])
    err_adp = stderr_from_sigmas(sigmas_fin[1], sigmas_inf[0])

    print(f"ATP ΔfG° = {standard_dgf_mu[0]:.2f} ± {err_atp:.2f} kJ/mol (1σ)")
    print(f"ADP ΔfG° = {standard_dgf_mu[1]:.2f} ± {err_adp:.2f} kJ/mol (1σ)")



    # we now apply the Legendre transform to convert from the standard ΔGf to the standard ΔG'f
    delta_dgf_list = np.array([
        cpd.transform(cc.p_h, cc.ionic_strength, cc.temperature, cc.p_mg).m_as("kJ/mol")
        for cpd in compound_list
    ])
    standard_dgf_prime_mu = standard_dgf_mu + delta_dgf_list

    # to create the formation energy covariance matrix, we need to combine the two outputs
    # sigma_fin and sigma_inf
    standard_dgf_cov = sigmas_fin @ sigmas_fin.T + 1e6 * sigmas_inf @ sigmas_inf.T

    print(f"μ(ΔGf'0) in kJ / mol: {standard_dgf_prime_mu}")
    print(f"Σ(ΔGf'0) in kJ^2 / mol^2: {standard_dgf_cov}")
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
