import ast

import pandas as pd


def load_mace_spectrum_data(file_path):
    """
    Load MACE spectrum data from a semicolon-separated file.

    Parameters:
    -----------
    file_path : str
        Path to the data file

    Returns:
    --------
    pandas.DataFrame
        DataFrame containing the parsed data
    """
    # Initialize lists for each column
    compound_ids = []
    smiles_list = []
    d_frees = []
    d_enthalpies = []
    d_entropies = []
    frees = []
    enthalpies = []
    entropies = []
    vib_energies_list = []

    with open(file_path, 'r') as f:
        for line in f:
            # Split by semicolon and strip whitespace
            parts = [part.strip() for part in line.split(';')]

            if len(parts) < 9:  # Skip incomplete lines
                continue

            # Parse data from each part
            compound_ids.append(parts[0])
            smiles_list.append(parts[1])

            try:
                d_frees.append(float(parts[2]))
                d_enthalpies.append(float(parts[3]))
                d_entropies.append(float(parts[4]))
                frees.append(float(parts[5]))
                enthalpies.append(float(parts[6]))
                entropies.append(float(parts[7]))

                # Parse vibrational energies safely
                try:
                    vib_energies = ast.literal_eval(parts[8])
                    vib_energies_list.append(vib_energies)
                except (SyntaxError, ValueError):
                    vib_energies_list.append(parts[8])
            except ValueError:
                continue

    # Create DataFrame
    df = pd.DataFrame({
        'compound_id': compound_ids,
        'smiles': smiles_list,
        'd_free': d_frees,
        'd_enthalpy': d_enthalpies,
        'd_entropy': d_entropies,
        'free': frees,
        'enthalpy': enthalpies,
        'entropy': entropies,
        'vib_energies': vib_energies_list
    })

    return df


if __name__ == "__main__":
    print(flush=True)
    data_path = r"../../data/CBRdb_C.dat"
    df = load_mace_spectrum_data(data_path)
    # Drop duplicate compound_id entries, keeping the first occurrence
    df = df.drop_duplicates(subset='compound_id', keep='first')
    out_path = r"CBRdb_C_mace_spectrum.csv.gz"
    df.to_csv(out_path, index=False)
