from functools import partial

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import CBRdb

plt.rcParams['axes.linewidth'] = 2.0

if __name__ == "__main__":
    method = 'drfp'
    sample = False
    assert method in ['rdkit', 'drfp'], "Method must be 'rdkit' or 'drfp'"

    data_c = pd.read_csv('../CBRdb_C.csv', low_memory=False)
    data_r = pd.read_csv('../CBRdb_R.csv', low_memory=False)[['id', 'reaction', 'smarts']]

    if sample:
        data_r = data_r.sample(n=10_000, random_state=1).reset_index(drop=True)

    fp_data = np.array(CBRdb.mp_calc(CBRdb.get_rxn_fingerprint_drfp, data_r['smarts'].tolist()))
    fp_data = fp_data.reshape((len(fp_data), -1))
    func_sim_drfp = partial(CBRdb.tanimoto_batch_drfp, db_bits=fp_data)
    sim = np.array(CBRdb.mp_calc(func_sim_drfp, fp_data.copy()))
    np.fill_diagonal(sim, np.nan)
    sim = sim.flatten()
    sim = sim[~np.isnan(sim)]

    CBRdb.plot_histogram(sim, xlab='Similarity score')
    plt.savefig('CBRdb_R_sim_hist.png', dpi=300)
    plt.savefig('CBRdb_R_sim_hist.pdf', dpi=300)
    plt.show()

    CBRdb.plot_histogram(sim, xlab='Similarity score')
    plt.yscale('log')
    plt.savefig('CBRdb_R_sim_log_hist.png', dpi=300)
    plt.savefig('CBRdb_R_sim_log_hist.pdf', dpi=300)
    plt.show()
