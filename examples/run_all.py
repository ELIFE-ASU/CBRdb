import CBRdb

if __name__ == "__main__":
    print("Program started", flush=True)
    f_man = False

    # Downloading all the data, the mol files and the (full) web pages
    print("Downloading all the data, the mol files and the (full) web pages", flush=True)
    CBRdb.download_data(target="C")
    CBRdb.download_data(target="C_full")
    CBRdb.download_data(target="R")
    # Provides the raw data files
    print("Done! \n", flush=True)

    # Manually adding the compounds
    if f_man:
        # Only run this file to construct a list of compounds to manually add to the database
        print("Manually adding the compounds", flush=True)
        CBRdb.compounds_manual_add()  # provides the C_IDs_good.dat file
        # Add this point you need to manually add the reactions in C_IDs_good.dat and the compounds in C_IDs_manual.dat
        print("Done! \n", flush=True)

    # Converting mol into smiles and cleaning up the data
    print("Converting mol into smiles and cleaning up the data", flush=True)
    CBRdb.preprocess(target="C")
    CBRdb.preprocess(target="R")
    # Provides the kegg_data_C.csv.zip and kegg_data_R.csv.zip files
    print("Done! \n", flush=True)

    # Convert the atlas files into a more readable format
    print("Converting atlas files", flush=True)
    CBRdb.clean_kegg_atlas()
    CBRdb.clean_atlas(f_exclude_kegg=True)
    # Provides the atlas_data_kegg_R.csv.zip and atlas_data_R.csv.zip files
    print("Done! \n", flush=True)

    # Clean up the data and remove the shortcut reactions
    print("Getting shortcut reactions", flush=True)
    CBRdb.clean_reaction_shortcuts()
    # Provides R_IDs_bad.dat
    print("Done! \n", flush=True)

    # Fix the problems with the reactions with the halogens
    print("Fixing the halogen compounds and reactions", flush=True)
    # Fix the halogen compounds in the C data
    cids_dict, smis_dict = CBRdb.fix_halogen_compounds()
    # Merge the halogen compounds into the C data
    CBRdb.merge_halogen_compounds(cids_dict, smis_dict)
    # Fix the halogen reactions in the R data
    CBRdb.fix_halogen_reactions(cids_dict, r_id_file="../data/atlas_data_R.csv.zip")
    CBRdb.fix_halogen_reactions(cids_dict, r_id_file="../data/atlas_data_kegg_R.csv.zip")
    CBRdb.fix_halogen_reactions(cids_dict, r_id_file="../data/kegg_data_R.csv.zip")
    print("Done! \n", flush=True)

    # Fix the reactions data
    print("Fixing the reactions data", flush=True)
    CBRdb.fix_reactions_data(r_file="../data/kegg_data_R.csv.zip")
    CBRdb.fix_reactions_data(r_file="../data/atlas_data_kegg_R.csv.zip")
    CBRdb.fix_reactions_data(r_file="../data/atlas_data_R.csv.zip")
    # Provides the processed R files
    print("Done! \n", flush=True)

    # Fix the reactions data with a n

    # Merge the data from the atlas and the kegg data
    CBRdb.merge_data(merge_col='reaction',
                     f_keep='first',
                     kegg_file="../data/atlas_data_kegg_R_processed.csv.zip",
                     atlas_file="../data/atlas_data_R_processed.csv.zip",
                     out_file="../data/atlas_kegg_processed_merged.csv.zip")

    CBRdb.merge_data(merge_col='reaction',
                     f_keep='first',
                     kegg_file="../data/kegg_data_R_processed.csv.zip",
                     atlas_file="../data/atlas_data_R_processed.csv.zip",
                     out_file="../data/kegg_atlas_processed_merged.csv.zip")
    # Fix the EC ids
    CBRdb.fix_ec_ids(input_file="../data/atlas_kegg_processed_merged.csv.zip")
    CBRdb.fix_ec_ids(input_file="../data/kegg_atlas_processed_merged.csv.zip")
    # Provides the full_processed_merged.csv.zip file
    print("Done! \n", flush=True)

    print("Program finished", flush=True)
