import os
import sys

import pandas as pd
from rdkit import Chem as Chem
import CBRdb

if __name__ == "__main__":
    print(flush=True)
    # Grab the command line arguments
    i = int(sys.argv[1])
    file = str(sys.argv[2])

    base_dir = os.getcwd()
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)

    data_dir = os.path.join(base_dir, f"{file}")
    out_file = os.path.join(base_dir, f"{file}_out.dat")

    # Load the data
    df = pd.read_csv(data_dir)

    # Select data that does not contain a star in the smiles
    df = df[~df['smiles'].str.contains(r'\*')]

    # Sort by n_heavy_atoms
    df = df.sort_values(by='n_heavy_atoms', ascending=True)
    n_heavy_max = 100
    n_heavy_min = 1
    df = df[df['n_heavy_atoms'] <= n_heavy_max]
    df = df[df['n_heavy_atoms'] >= n_heavy_min]

    id_list = df['compound_id'].to_list()
    smi_list = df['smiles'].to_list()

    mol = Chem.MolFromSmiles(smi_list[i])
    (d_free,
     d_enthalpy,
     d_entropy,
     free,
     enthalpy,
     entropy,
     vib_energies) = CBRdb.calculate_free_energy_formation_mace(mol)

    # Write the data to the shared file
    out_str = f"{id_list[i]}; {smi_list[i]}; {d_free}; {d_enthalpy}; {d_entropy}; {free}; {enthalpy}; {entropy}; {vib_energies}  \n"
    CBRdb.write_to_shared_file(out_str, out_file)
