import pandas as pd

if __name__ == "__main__":
    print("Program started", flush=True)
    # load the data
    data = pd.read_csv("../data/atlas_kegg_processed_merged_deduped.csv.zip", compression='zip') # , index_col=0
    print("data loaded", flush=True)
    print("data shape", data.shape, flush=True)
    print("data head", data.head(), flush=True)
    # print the data columns
    print("data columns", data.columns, flush=True)

    # get the entry with the id R04254
    print(data[data["id"] == "R04254"], flush=True)


    # load to compound data
    data_c = pd.read_csv("../data/kegg_data_C.csv.zip")
    print("data columns", data_c.columns, flush=True)
    print(data_c[data_c["compound_id"] == "C06143"].values, flush=True)
    print("Program end", flush=True)