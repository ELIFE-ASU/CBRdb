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

    temperatures_k = [
        223.15, 273.15, 298.15, 323.15, 348.15, 373.15, 398.15, 423.15,
        448.15, 473.15, 498.15, 523.15, 548.15, 573.15, 598.15, 623.15,
        223.15, 273.15, 298.15, 323.15, 348.15, 373.15, 398.15, 423.15,
        448.15, 473.15, 498.15, 523.15, 548.15, 573.15, 598.15, 623.15
    ]
    pressures_pa = [
        100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 101322.0, 232014.4, 475716.9,
        891804.9, 1553649.9, 2547860.3, 3973649.3, 5943125.1, 8583784.3, 12045757.2, 16521128.9,
        100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 100000.0,
        100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 100000.0, 100000.0
    ]

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
    print(f'{smi_list[i]}', flush=True)
    mol = Chem.MolFromSmiles(smi_list[i])
    free_energy, _ = CBRdb.calculate_free_energy_formation_mace(mol,
                                                                temperature=temperatures_k,
                                                                pressure=pressures_pa)

    # Write the data to the shared file
    out_str = f"{id_list[i]}; {smi_list[i]}; {free_energy}\n"
    CBRdb.write_to_shared_file(out_str, out_file)
