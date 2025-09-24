from functools import partial

import matplotlib.pyplot as plt
import pandas as pd

import CBRdb

plt.rcParams['axes.linewidth'] = 2.0

if __name__ == "__main__":
    print(flush=True)
    data_c = pd.read_csv('../CBRdb_C.csv', low_memory=False)
    data_r = pd.read_csv('../CBRdb_R.csv', low_memory=False)

    # Trim them down to only the necessary columns
    data_r = pd.DataFrame(data_r[['id', 'reaction']])

    # if the smarts column does not exist, create it
    if 'smarts' not in data_r.columns:
        print("Generating SMARTS for reactions...", flush=True)
        func_smarts = partial(CBRdb.to_smarts_rxn_line, data_c=data_c, add_stoich=False)
        data_r['smarts'] = CBRdb.mp_calc(func_smarts, data_r['reaction'].tolist())

    print("Generating fingerprints for reactions...", flush=True)
    data_r['fp_struct'] = CBRdb.mp_calc(CBRdb.get_rxn_fingerprint, data_r['smarts'].tolist())

    # Only select the reactions where the id starts with 'R'
    sel = data_r['id'].str.startswith('R')
    # Split the data into KEGG and ATLAS reactions
    data_r_atlas = data_r[~sel]
    data_r_kegg = data_r[sel]

    # print the number of reactions
    print(f"Number of reactions (KEGG) : {len(data_r_kegg)}")
    print(f"Number of reactions (ATLAS): {len(data_r_atlas)}")

    fp_kegg = data_r_kegg['fp_struct'].tolist()
    fp_atlas = data_r_atlas['fp_struct'].tolist()

    func_sim = partial(CBRdb.find_max_similar_rxn, rxn_db=fp_kegg)
    sim_max, sim_max_idx = zip(*CBRdb.mp_calc(func_sim, fp_atlas))

    # convert the idx into the corresponding reaction id
    sim_max_id = [data_r_kegg.iloc[idx]['id'] for idx in sim_max_idx]

    # add the results to the atlas dataframe
    data_r_atlas['sim_max'] = sim_max
    data_r_atlas['sim_max_id'] = sim_max_id

    # Print the lowest 10 similarity scores
    print("Lowest 10 similarity scores:", flush=True)
    print(data_r_atlas.nsmallest(10, 'sim_max')[['id', 'sim_max', 'sim_max_id']], flush=True)

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
