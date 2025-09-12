import numpy as np
import pubchempy as pcp
from equilibrator_api import ComponentContribution, Q_

import warnings
warnings.filterwarnings("error")


def get_dg_prime(eq_line, conditions='nstandard'):
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
        error= np.nan
    return val, error, rev

if __name__ == "__main__":
    print(flush=True)
    eq_line = '2 C00084 <=> 1 C00466'
    print(get_dg_prime(eq_line))