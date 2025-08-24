import os
import sys

import pandas as pd

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
    n_heavy_max = 50
    n_heavy_min = 2
    df = df[df['n_heavy_atoms'] <= n_heavy_max]
    df = df[df['n_heavy_atoms'] >= n_heavy_min]

    id_list = df['compound_id'].to_list()
    smi_list = df['smiles'].to_list()

    atoms, charge, multiplicity = CBRdb.smi_to_atoms(smi_list[i])
    # Calculate the IR spectrum
    data_ir, data_raman, _ = CBRdb.calculate_vib_spectrum(atoms,
                                                          charge=charge,
                                                          multiplicity=multiplicity,
                                                          n_procs=1)

    if data_ir is None or data_raman is None:
        print(f"Skipping molecule {i} due to missing data.", flush=True)
    else:
        # Add them to the list
        freq = data_ir['Frequency (cm^-1)'].to_list()
        eps_ir = data_ir['Epsilon'].to_list()
        int_ir = data_ir['Intensity (km/mol)'].to_list()
        int_ram = data_raman['Activity (A^4 amu^-1)'].to_list()
        dep_ram = data_raman['Depolarization'].to_list()
        # Write the data to the shared file
        out_str = f"{id_list[i]}; {smi_list[i]}; {freq}; {eps_ir}; {int_ir}; {int_ram}; {dep_ram} \n"
        CBRdb.write_to_shared_file(out_str, out_file)
