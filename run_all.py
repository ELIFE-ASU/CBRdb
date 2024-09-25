import atlas_converter
import compounds_manual_add
import fix_reactions_data
import lets_get_kegg
import preprocessor
import clean_reaction_shortcuts

if __name__ == "__main__":
    print("Program started", flush=True)
    # f_man = False
    # # Downloading all the data, the mol files and the (full) web pages
    # lets_get_kegg.main(target="C")
    # lets_get_kegg.main(target="C_full")
    # lets_get_kegg.main(target="R")
    # # Manually adding the compounds
    # if f_man:
    #     compounds_manual_add.main()
    #     # Add this point you need to manually add the reactions in good_cids.dat
    #
    # # Converting mol into smiles and cleaning up the data
    # preprocessor.main(target="C")
    # preprocessor.main(target="R")
    # # Convert the atlas files into a more readable format
    # # atlas_converter.main()
    # # Fix the reactions data
    # fix_reactions_data.main(eq_file="Data/kegg_data_R.csv.zip")
    fix_reactions_data.main(eq_file="Data/atlas_data_kegg_R.csv.zip")
    fix_reactions_data.main(eq_file="Data/atlas_data_R.csv.zip")

    # Clean up the data and remove the shortcut reactions
    clean_reaction_shortcuts.main()

    # To do
    # Address problem with halogens
    #

    print("Program finished", flush=True)
