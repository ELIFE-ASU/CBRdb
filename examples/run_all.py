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

    if f_man:     # Only run this file to construct a list of compounds to manually add to the database
        print("Manually adding compounds...", flush=True)
        CBRdb.compounds_manual_add()  # generates C_IDs_good.dat
        # At this point you need to manually add the reactions in C_IDs_good.dat and the compounds in C_IDs_manual.dat
    
    print("Preprocessing compounds (metadata and structure)...", flush=True)
    C_meta, C_main = CBRdb.preprocess(target="C")  # generates kegg_data_C.csv kegg_data_C_metadata.csv + kegg_data_C_dupemap.csv

    print("Importing and standardizing ATLAS reactions...", flush=True)
    atlas_data_R = CBRdb.clean_atlas(f_exclude_kegg=True)  # generates atlas_data_R.csv

    print("Preprocessing KEGG reactions...", flush=True)
    kegg_data_R = CBRdb.preprocess(target="R")  # generates kegg_data_R.csv

    print("Finding and removing suspect KEGG reactions...", flush=True)
    kegg_data_R = CBRdb.remove_suspect_reactions()  # adds to R_IDs_bad.dat, modifies kegg_data_R.csv

    print("Instantiating specific halogen compounds from generic ones...", flush=True)
    cids_dict, smis_dict = CBRdb.fix_halogen_compounds()    # turn generic halogens in the C data into specific halogens
    C_main = CBRdb.merge_halogen_compounds(cids_dict, smis_dict)    # merge these into the C data, modifies kegg_data_C.csv
    
    print("Fixing halogen reactions in KEGG and ATLAS...", flush=True)
    kegg_data_R = CBRdb.fix_halogen_reactions(cids_dict, r_id_file=kegg_reactions_data + ".csv")  # modifies kegg_data_R.csv
    atlas_data_R = CBRdb.fix_halogen_reactions(cids_dict, r_id_file=atlas_reactions_data + ".csv")  # modifies atlas_data_R.csv

    print("De-duplicating compounds by structure, then propagating to reactions...", flush=True)
    datasets = CBRdb.merge_data_sets.dedupe_compounds()  # modifies all compound and reaction data files above; returns de-duped versions of them
    
    print("Fixing the reactions data...", flush=True) 
    datasets['kegg_data_R_processed'] = CBRdb.fix_reactions_data(r_file=kegg_reactions_data + "_dedupedCs.csv")  
    datasets['atlas_data_R_processed'] = CBRdb.fix_reactions_data(r_file=atlas_reactions_data + "_dedupedCs.csv")

    print("Combining and deduplicating the reactions data...", flush=True)
    datasets['reactions_joined'] = pd.concat(objs=[datasets['kegg_data_R_processed'], datasets['atlas_data_R_processed']], ignore_index=True)
    datasets['r_dupemap'] = CBRdb.tools_eq.generate_reaction_dupemap(datasets['reactions_joined'], prefix='T')
    datasets['reactions_merged'] = CBRdb.merge_duplicate_reactions(datasets['reactions_joined'], datasets['r_dupemap'])
    datasets['reactions_merged'].to_csv(out_reactions_data + ".csv", encoding='utf-8', index=False)

    print("Program finished", flush=True)
