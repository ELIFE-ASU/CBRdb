import os
import sys

import numpy as np
from rdkit import Chem as Chem

import CBRdb

if __name__ == "__main__":
    print(flush=True)
    timeout = 3.0 * 60.0 * 60.0  # 3 hours
    # Grab the command line arguments
    i = int(sys.argv[1])
    file = str(sys.argv[2])
    calc_type = str(sys.argv[3])

    base_dir = os.path.abspath(os.path.expanduser(f"../data/"))
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)

    data_dir = os.path.join(base_dir, f"{file}.dat")
    out_file = os.path.join(base_dir, f"{file}_out_{calc_type}.dat")

    # Load the data
    smi = np.loadtxt(data_dir, dtype=str, comments=None)[i]

    # Get the smiles
    mol = Chem.AddHs(Chem.MolFromSmiles(smi, sanitize=True))

    # Calculate the IR spectrum
    data = CBRdb.calculate_vib_spectrum(mol, calc_type=calc_type, n_procs=1)

    # Add them to the list
    freq = data['Frequency (cm^-1)'].to_list()
    intensity = data['Intensity (km/mol)'].to_list()

    smi = Chem.MolToSmiles(mol, isomericSmiles=True, kekuleSmiles=True, canonical=True)
    smi = CBRdb.standardise_smiles(smi)
    # Write the data to the shared file
    CBRdb.write_to_shared_file(f"{i}; {smi}; {freq}; {intensity}\n", out_file)
