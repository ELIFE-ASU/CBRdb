import lets_get_kegg
import preprocessor
import atlas_converter
import fix_reactions_data

if __name__ == "__main__":
    print("Program started", flush=True)
    lets_get_kegg.main(target="C")
    lets_get_kegg.main(target="R")

    preprocessor.main(target="C")
    preprocessor.main(target="R")

    atlas_converter.main()

    fix_reactions_data.main()

    print("Program finished", flush=True)