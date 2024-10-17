import atlas_converter
import clean_reaction_shortcuts
import compounds_manual_add
import fix_reactions_data
import lets_get_kegg
import merge_data_sets
import preprocessor

if __name__ == "__main__":
    print("Program started", flush=True)
    f_man = False
    # Downloading all the data, the mol files and the (full) web pages
    print("Downloading all the data, the mol files and the (full) web pages", flush=True)
    lets_get_kegg.main(target="C")
    lets_get_kegg.main(target="C_full")
    lets_get_kegg.main(target="R")
    # Provides the raw data files
    print("Done! \n", flush=True)

    # Manually adding the compounds
    if f_man:
        print("Manually adding the compounds", flush=True)
        compounds_manual_add.main()  # provides the C_IDs_good.dat file
        # Add this point you need to manually add the reactions in C_IDs_good.dat and the compounds in C_IDs_manual.dat
        print("Done! \n", flush=True)

    # Converting mol into smiles and cleaning up the data
    print("Converting mol into smiles and cleaning up the data", flush=True)
    preprocessor.main(target="C")
    preprocessor.main(target="R")
    # Provides the kegg_data_C.csv.zip and kegg_data_R.csv.zip files
    print("Done! \n", flush=True)

    # Convert the atlas files into a more readable format
    print("Converting atlas files", flush=True)
    atlas_converter.main()
    # Provides the atlas_data_kegg_R.csv.zip and atlas_data_R.csv.zip files
    print("Done! \n", flush=True)

    # Clean up the data and remove the shortcut reactions
    print("Getting shortcut and glycan reactions", flush=True)
    clean_reaction_shortcuts.main()
    # Provides R_IDs_bad.dat
    print("Done! \n", flush=True)

    # Fix the problems with the reactions with the halogens
    # clean_reaction_shortcuts.fix_halogen_reactions()

    # Fix the reactions data
    print("Fixing the reactions data", flush=True)
    fix_reactions_data.main(r_file="Data/kegg_data_R.csv.zip")
    fix_reactions_data.main(r_file="Data/atlas_data_kegg_R.csv.zip")
    fix_reactions_data.main(r_file="Data/atlas_data_R.csv.zip")
    # Provides the processed R files
    print("Done! \n", flush=True)

    # Merge the data from the atlas and the kegg data
    merge_data_sets.merge_data(merge_col='reaction',
                               f_keep='first',
                               kegg_file="Data/atlas_data_kegg_R_processed.csv.zip",
                               atlas_file="Data/atlas_data_R_processed.csv.zip",
                               out_file="Data/atlas_kegg_processed_merged.csv.zip")

    merge_data_sets.merge_data(merge_col='reaction',
                               f_keep='first',
                               kegg_file="Data/kegg_data_R_processed.csv.zip",
                               atlas_file="Data/atlas_data_R_processed.csv.zip",
                               out_file="Data/kegg_atlas_processed_merged.csv.zip")
    # Fix the EC ids
    merge_data_sets.fix_ec_ids(input_file="Data/atlas_kegg_processed_merged.csv.zip")
    merge_data_sets.fix_ec_ids(input_file="Data/kegg_atlas_processed_merged.csv.zip")
    # Provides the full_processed_merged.csv.zip file
    print("Done! \n", flush=True)

    print("Program finished", flush=True)
