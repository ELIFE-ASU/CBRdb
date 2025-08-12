try:
    import CBRdb
except ImportError:
    import sys
    import os

    sys.path.append(os.path.abspath('../'))
    import CBRdb

import pandas as pd

if __name__ == "__main__":
    print("Program started", flush=True)
    f_man = False
    f_download = False
    kegg_reactions_data = "../data/kegg_data_R"
    atlas_reactions_data = "../data/atlas_data_R"
    out_reactions_data = "../data/CBRdb_R"

    if f_download:
        print("Downloading all the data, the mol files and the (full) web pages", flush=True)
        C_attrs = CBRdb.download_data(target="C_full")
        C_attrs.query('~ATOM|~BOND')[[]].rename_axis('id').assign(reason='molless').to_csv('../data/C_ids_bad.dat',
                                                                                           mode='a')
        pd.read_csv('../data/C_ids_bad.dat').drop_duplicates().sort_values(by='id').to_csv('../data/C_ids_bad.dat',
                                                                                           index=False)
        C_mols = CBRdb.download_data(target="C", skip=C_attrs.query('~ATOM|~BOND').index)
        R_attrs = CBRdb.download_data(target="R")
    else:
        C_attrs = CBRdb.lets_get_kegg.log_attributes_found('../../data/kegg_data_C_full')
        C_mols = pd.DataFrame(
            {f.split('/')[-2]: {'mol_file': f} for f in CBRdb.file_list_all('../../data/kegg_data_C') if
             f.endswith('.mol')}).T
        R_attrs = CBRdb.lets_get_kegg.log_attributes_found('../../data/kegg_data_R')

    print("Preprocessing compounds (metadata and structure)...", flush=True)
    C_main = CBRdb.preprocess(target="C")  # generates kegg_data_C.csv

    print("Importing and standardizing ATLAS reactions...", flush=True)
    atlas_data_R = CBRdb.clean_atlas(f_exclude_kegg=False)  # generates atlas_data_R.csv

    print("Preprocessing KEGG reactions...", flush=True)
    kegg_data_R = CBRdb.preprocess(target="R")  # generates kegg_data_R.csv

    print("Instantiating specific halogen compounds from generic ones...", flush=True)
    cids_dict, smis_dict, specific_halogens = CBRdb.fix_halogen_compounds()  # turn generic halogens in the C data into specific halogens
    C_main = CBRdb.merge_halogen_compounds(cids_dict,
                                           smis_dict)  # merge these into the C data, modifies kegg_data_C.csv

    print("Fixing halogen reactions in KEGG and ATLAS...", flush=True)
    kegg_data_R = CBRdb.fix_halogen_reactions_without_existing_halogens(kegg_data_R, C_main, specific_halogens)
    atlas_data_R = CBRdb.fix_halogen_reactions_without_existing_halogens(atlas_data_R, C_main, specific_halogens)
    CBRdb.reaction_csv(kegg_data_R, '../data/kegg_data_R.csv')
    CBRdb.reaction_csv(atlas_data_R, '../data/atlas_data_R.csv')

    print("Finding and removing suspect KEGG reactions + de-duping compound dataset for compounds and reactions...",
          flush=True)
    dbs = CBRdb.iteratively_prune_entries(kegg_data_R, atlas_data_R, C_main,
                                          to_quarantine="shortcut|structure_missing")

    print("Fixing the reactions data...", flush=True)
    dbs['kegg_data_R_processed'] = CBRdb.fix_reactions_data(r_file=kegg_reactions_data + "_dedupedCs.csv",
                                                            c_file="../CBRdb_C.csv")
    dbs['atlas_data_R_processed'] = CBRdb.fix_reactions_data(r_file=atlas_reactions_data + "_dedupedCs.csv",
                                                             c_file="../CBRdb_C.csv")

    print("Combining and deduplicating the reactions data...", flush=True)
    dbs['reactions_joined'] = pd.concat(objs=[dbs['kegg_data_R_processed'], dbs['atlas_data_R_processed']],
                                        ignore_index=True)
    dbs['r_dupemap'] = CBRdb.tools_eq.generate_reaction_dupemap(dbs['reactions_joined'], prefix='T')
    dbs['CBRdb_R'] = CBRdb.merge_duplicate_reactions(dbs['reactions_joined'], dbs['r_dupemap'])
    CBRdb.reaction_csv(dbs['CBRdb_R'], "../CBRdb_R.csv")
    CBRdb.add_R_col_to_C_file(final_output_Cs_fp='../CBRdb_C.csv', final_output_Rs_fp='../CBRdb_R.csv')

    print("Program finished", flush=True)
