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
    sel = data_r['id'].str.startswith('R')
    data_r = data_r[sel]

    # randomly select X reactions from data_r for testing
    # data_r = data_r.sample(n=10, random_state=1).reset_index(drop=True)

    # if the smarts column does not exist, create it
    if 'smarts' not in data_r.columns:
        print("Generating SMARTS for reactions...", flush=True)
        func_smarts = partial(CBRdb.to_smarts_rxn_line, data_c=data_c, add_stoich=False)
        data_r['smarts'] = CBRdb.mp_calc(func_smarts, data_r['reaction'].tolist())

    # Print the number of reactions
    print(f"Number of reactions : {len(data_r)}", flush=True)
    data_r_atlas = data_r.copy()
    data_r_kegg = data_r.copy()

    print("Generating fingerprints for reactions...", flush=True)
    fp_kegg = CBRdb.mp_calc(CBRdb.get_rxn_fingerprint_drfp, data_r_kegg['smarts'].tolist())
    fp_atlas = CBRdb.mp_calc(CBRdb.get_rxn_fingerprint_drfp, data_r_atlas['smarts'].tolist())

    # reshape the fingerprints into a 2D array
    fp_kegg = np.array(fp_kegg)
    fp_atlas = np.array(fp_atlas)
    fp_kegg = fp_kegg.reshape((len(fp_kegg), -1))
    fp_atlas = fp_atlas.reshape((len(fp_atlas), -1))

    # # set the diagnal to zero
    # np.fill_diagonal(fp_atlas, 0)

    print("KEGG fingerprints shape:", np.shape(fp_kegg), flush=True)
    print("ATLAS fingerprints shape:", np.shape(fp_atlas), flush=True)

    print("Calculating similarities...", flush=True)

    func_sim_drfp = partial(CBRdb.tanimoto_batch_drfp, db_bits=fp_kegg)
    sim = CBRdb.mp_calc(func_sim_drfp, fp_atlas)
    sim = np.array(sim)
    # set the diagnal to 0
    np.fill_diagonal(sim, np.nan)
    print("Similarities array shape:", np.shape(sim), flush=True)
    # flatten the sim array
    sim = sim.flatten()

    # remove nan values
    sim = sim[~np.isnan(sim)]

    CBRdb.plot_histogram(sim, xlab='Similarity score')
    plt.savefig('KEGG_R_sim_hist.png', dpi=300)
    plt.savefig('KEGG_R_sim_hist.pdf', dpi=300)
    plt.show()

    CBRdb.plot_histogram(sim, xlab='Similarity score')
    plt.yscale('log')
    plt.savefig('KEGG_R_sim_log_hist.png', dpi=300)
    plt.savefig('KEGG_R_sim_log_hist.pdf', dpi=300)
    plt.show()
