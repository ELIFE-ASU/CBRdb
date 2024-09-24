import atlas_converter
import fix_reactions_data
import lets_get_kegg
import preprocessor

if __name__ == "__main__":
    print("Program started", flush=True)
    # Downloading all the data
    lets_get_kegg.main(target="C")
    lets_get_kegg.main(target="C_full")
    lets_get_kegg.main(target="R")
    # Converting all the files into a more readable format
    preprocessor.main(target="C")
    preprocessor.main(target="R")
    # Convert the atlas files into a more readable format
    atlas_converter.main()
    # Fix the reactions data
    fix_reactions_data.main()

    print("Program finished", flush=True)
