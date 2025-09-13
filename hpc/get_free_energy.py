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


if __name__ == "__main__":
    print(flush=True)

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
