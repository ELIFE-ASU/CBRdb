import os

import numpy as np
from rdkit import Chem as Chem

import CBRdb

if __name__ == "__main__":
    print(flush=True)

    file = "CBRdb_C"

    base_dir = os.path.abspath(os.path.expanduser(f"../data/"))
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)

    data_dir = os.path.join(base_dir, f"{file}.dat")
    out_file = os.path.join(base_dir, f"{file}_out.dat")

    # Load the data
    smi_list = np.loadtxt(data_dir, dtype=str, comments=None)

    print(f"Number of SMILES: {len(smi_list)}", flush=True)

    # Get the smiles
    mols = [Chem.AddHs(Chem.MolFromSmiles(smi, sanitize=True)) for smi in smi_list]

    for i in range(len(smi_list)):
        print(f"Processing molecule: {i + 1}/{len(smi_list)}, smi: {smi_list[i]}", flush=True)
        mol = mols[i]

        # Calculate the IR spectrum
        data = CBRdb.calculate_vib_spectrum(mol)

        # Add them to the list
        freq = data['Frequency (cm^-1)'].to_list()
        intensity = data['Intensity (km/mol)'].to_list()

        smi = Chem.MolToSmiles(mol, isomericSmiles=True, kekuleSmiles=True, canonical=True)
        # Write the data to the shared file
        CBRdb.write_to_shared_file(f"{i}; {smi}; {freq}; {intensity}\n", out_file)
