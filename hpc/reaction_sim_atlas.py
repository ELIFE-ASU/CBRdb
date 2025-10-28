from functools import partial

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import CBRdb

plt.rcParams['axes.linewidth'] = 2.0

if __name__ == "__main__":
    print(flush=True)
    method = 'drfp'
    assert method in ['rdkit', 'drfp'], "Method must be 'rdkit' or 'drfp'"

    data_c = pd.read_csv('../CBRdb_C.csv', low_memory=False)
    data_r = pd.read_csv('../CBRdb_R.csv', low_memory=False)

    # Trim them down to only the necessary columns
    data_r = pd.DataFrame(data_r[['id', 'reaction']])

    # randomly select X reactions from data_r for testing
    # data_r = data_r.sample(n=10_000, random_state=1).reset_index(drop=True)

    # if the smarts column does not exist, create it
    if 'smarts' not in data_r.columns:
        print("Generating SMARTS for reactions...", flush=True)
        func_smarts = partial(CBRdb.to_smarts_rxn_line, data_c=data_c, add_stoich=False)
        data_r['smarts'] = CBRdb.mp_calc(func_smarts, data_r['reaction'].tolist())

    # Only select the reactions where the id starts with 'R'
    sel = data_r['id'].str.startswith('R')
    # Split the data into KEGG and ATLAS reactions
    data_r_atlas = data_r[~sel].copy()
    data_r_kegg = data_r[sel].copy()

    # Print the number of reactions
    print(f"Number of reactions (KEGG) : {len(data_r_kegg)}", flush=True)
    print(f"Number of reactions (ATLAS): {len(data_r_atlas)}", flush=True)

    print("Generating fingerprints for reactions...", flush=True)
    if method == 'rdkit':
        fp_kegg = CBRdb.mp_calc(CBRdb.get_rxn_fingerprint, data_r_kegg['smarts'].tolist())
        fp_atlas = CBRdb.mp_calc(CBRdb.get_rxn_fingerprint, data_r_atlas['smarts'].tolist())
    else:
        fp_kegg = CBRdb.mp_calc(CBRdb.get_rxn_fingerprint_drfp, data_r_kegg['smarts'].tolist())
        fp_atlas = CBRdb.mp_calc(CBRdb.get_rxn_fingerprint_drfp, data_r_atlas['smarts'].tolist())

        # reshape the fingerprints into a 2D array
        fp_kegg = np.array(fp_kegg)
        fp_atlas = np.array(fp_atlas)
        fp_kegg = fp_kegg.reshape((len(fp_kegg), -1))
        fp_atlas = fp_atlas.reshape((len(fp_atlas), -1))

    print("KEGG fingerprints shape:", np.shape(fp_kegg), flush=True)
    print("ATLAS fingerprints shape:", np.shape(fp_atlas), flush=True)

    print("Calculating similarities...", flush=True)
    if method == 'rdkit':
        func_sim = partial(CBRdb.find_max_similar_rxn, rxn_db=fp_kegg)
        sim_max, sim_max_idx = zip(*CBRdb.mp_calc(func_sim, fp_atlas))
    else:
        func_sim_drfp = partial(CBRdb.find_max_similar_rxn_drfp, rxn_db=fp_kegg, replace_unity=True)
        sim_max, sim_max_idx = zip(*CBRdb.mp_calc(func_sim_drfp, fp_atlas))
    print("Similarities calculated.", flush=True)

    # convert the idx into the corresponding reaction id
    sim_max_id = [data_r_kegg.iloc[idx]['id'] for idx in sim_max_idx]

    # add the results to the atlas dataframe
    data_r_atlas['sim_max'] = sim_max
    data_r_atlas['sim_max_id'] = sim_max_id

    # set nan and inf values to 0.0
    # data_r_atlas['sim_max'] = data_r_atlas['sim_max'].replace([np.nan, np.inf, None], 0.0)

    # Print the lowest 10 similarity scores
    print("Lowest 10 similarity scores:", flush=True)
    print(data_r_atlas.nsmallest(10, 'sim_max')[['id', 'sim_max', 'sim_max_id']], flush=True)

    # print cases where smin_max is 0.0, nan, or inf
    print("Cases where similarity is 0.0, nan, or inf:", flush=True)
    print(data_r_atlas[(data_r_atlas['sim_max'] == 0.0) | (data_r_atlas['sim_max'].isna()) | (
            data_r_atlas['sim_max'] == np.inf)][['id', 'sim_max', 'sim_max_id']], flush=True)
    print(flush=True)

    # print cases where sim_max is 1.0
    print("Cases where similarity is 1.0:", flush=True)
    print(data_r_atlas[data_r_atlas['sim_max'] == 1.0][['id', 'sim_max', 'sim_max_id']], flush=True)
    print(flush=True)

    CBRdb.plot_histogram(sim_max, xlab='Best Similarity score')
    plt.savefig('ATLAS_R_sim_hist.png', dpi=300)
    plt.savefig('ATLAS_R_sim_hist.pdf', dpi=300)
    plt.show()

    CBRdb.plot_histogram(sim_max, xlab='Best Similarity score')
    plt.yscale('log')
    plt.savefig('ATLAS_R_sim_log_hist.png', dpi=300)
    plt.savefig('ATLAS_R_sim_log_hist.pdf', dpi=300)
    plt.show()

    data_r_atlas.to_csv('ATLAS_R_sim.csv.gz', index=False, compression='gzip')
    print('Done.', flush=True)
    print(flush=True)
