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
        CBRdb.download_data(target="C")
        CBRdb.download_data(target="C_full")
        CBRdb.download_data(target="R")

    print("Preprocessing compounds (metadata and structure)...", flush=True)
    C_main = CBRdb.preprocess(target="C")  # generates kegg_data_C.csv

    print("Importing and standardizing ATLAS reactions...", flush=True)
    atlas_data_R = CBRdb.clean_atlas(f_exclude_kegg=True)  # generates atlas_data_R.csv

    print("Preprocessing KEGG reactions...", flush=True)
    kegg_data_R = CBRdb.preprocess(target="R")  # generates kegg_data_R.csv

    print("Instantiating specific halogen compounds from generic ones...", flush=True)
    cids_dict, smis_dict = CBRdb.fix_halogen_compounds()    # turn generic halogens in the C data into specific halogens
    C_main = CBRdb.merge_halogen_compounds(cids_dict, smis_dict)    # merge these into the C data, modifies kegg_data_C.csv

    print("Fixing halogen reactions in KEGG and ATLAS...", flush=True)
    kegg_data_R = CBRdb.fix_halogen_reactions(cids_dict, r_id_file=kegg_reactions_data + ".csv")  # modifies kegg_data_R.csv
    atlas_data_R = CBRdb.fix_halogen_reactions(cids_dict, r_id_file=atlas_reactions_data + ".csv")  # modifies atlas_data_R.csv

    print("Finding and removing suspect KEGG reactions + de-duping compound dataset for compounds and reactions...", flush=True)
    dbs = CBRdb.iteratively_prune_entries(kegg_data_R, atlas_data_R, C_main, to_quarantine="shortcut|structure_missing") #unclear|incomplete|general

    print("Fixing the reactions data...", flush=True) 
    dbs['kegg_data_R_processed'] = CBRdb.fix_reactions_data(r_file=kegg_reactions_data + "_dedupedCs.csv",
                                                                c_file="../CBRdb_C.csv")
    dbs['atlas_data_R_processed'] = CBRdb.fix_reactions_data(r_file=atlas_reactions_data + "_dedupedCs.csv",
                                                                 c_file="../CBRdb_C.csv")

    print("Combining and deduplicating the reactions data...", flush=True)
    dbs['reactions_joined'] = pd.concat(objs=[dbs['kegg_data_R_processed'], dbs['atlas_data_R_processed']], ignore_index=True)
    dbs['r_dupemap'] = CBRdb.tools_eq.generate_reaction_dupemap(dbs['reactions_joined'], prefix='T')
    dbs['reactions_merged'] = CBRdb.merge_duplicate_reactions(dbs['reactions_joined'], dbs['r_dupemap'])
    dbs['reactions_merged'].to_csv("../CBRdb_R.csv", encoding='utf-8', index=False)

    print("Program finished", flush=True)
