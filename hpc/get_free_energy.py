import os
import sys

import pandas as pd
from rdkit import Chem as Chem
import pubchempy as pcp
from equilibrator_api import ComponentContribution, Q_
# import CBRdb

# conda install -c conda-forge equilibrator-api
# conda install -c conda-forge pubchempy

if __name__ == "__main__":
    print(flush=True)
    smiles = "C1=CC=C(C=C1)C=O"  # example
    hits = pcp.get_compounds(smiles, namespace="smiles")
    print(hits[0].cid)

    cc = ComponentContribution()
    cc.p_h = Q_(7.4)
    cc.p_mg = Q_(3.0)
    cc.ionic_strength = Q_("0.25M")
    cc.temperature = Q_("298.15K")
    atpase_reaction = cc.parse_reaction_formula(
        "bigg.metabolite:atp + bigg.metabolite:h2o = "
        "bigg.metabolite:adp + bigg.metabolite:pi"
    )
    dG_prime = cc.standard_dg_prime(atpase_reaction)
    print(f"ΔG'° = {dG_prime}")

    dGm_prime = cc.physiological_dg_prime(atpase_reaction)
    print(f"ΔG'm = {dGm_prime}")

    dG_prime_value_in_kj_per_mol = dG_prime.value.m_as("kJ/mol")
    dG_prime_error_in_kj_per_mol = dG_prime.error.m_as("kJ/mol")
    print(
        f"ΔG'° = {dG_prime_value_in_kj_per_mol:.1f} +/- "
        f"{dG_prime_error_in_kj_per_mol:.1f} kJ/mol"
    )
    print(f"ln(Reversibility Index) = {cc.ln_reversibility_index(atpase_reaction)}")