import atlas_converter
import compounds_manual_add
import fix_reactions_data
import lets_get_kegg
import preprocessor
import clean_reaction_shortcuts

if __name__ == "__main__":
    print("Program started", flush=True)
    f_man = False
    # Downloading all the data, the mol files and the (full) web pages
    print("Downloading all the data, the mol files and the (full) web pages", flush=True)
    lets_get_kegg.main(target="C")
    lets_get_kegg.main(target="C_full")
    lets_get_kegg.main(target="R")
    # Manually adding the compounds
    if f_man:
        print("Manually adding the compounds", flush=True)
        compounds_manual_add.main() # provides the good_cids.dat file
        # Add this point you need to manually add the reactions in good_cids.dat and the compounds in manual_cids.dat

    # Converting mol into smiles and cleaning up the data
    print("Converting mol into smiles and cleaning up the data", flush=True)
    preprocessor.main(target="C")
    preprocessor.main(target="R")
    # provides the kegg_data_C.csv.zip and kegg_data_R.csv.zip files

    # Convert the atlas files into a more readable format
    atlas_converter.main()

    # Fix the reactions data
    fix_reactions_data.main(eq_file="Data/kegg_data_R.csv.zip")
    fix_reactions_data.main(eq_file="Data/atlas_data_kegg_R.csv.zip")
    fix_reactions_data.main(eq_file="Data/atlas_data_R.csv.zip")

    # Clean up the data and remove the shortcut reactions
    clean_reaction_shortcuts.main()

    # Merge the data from the atlas and the kegg data


    # To do
    # Log the data which was removed in the cleaning process of the reactions
    # Address problem with halogens
    # Get the clean reaction shortcuts working


    print("Program finished", flush=True)
