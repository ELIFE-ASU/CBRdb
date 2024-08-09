import os
import time
import numpy as np
import re

import pandas as pd

file_path = r"C:\Users\louie\skunkworks\data\atlas_reactions.csv"

# load the data using pandas
data = pd.read_csv(file_path)
# get the set of the ids
atlas_ids = sorted(list(set(data["rn"])))
atlas_list = []
eq_list = []
# loop over the atlas ids
for i, atlas_id in enumerate(atlas_ids):
    # print(f"atlas_ids {atlas_id}")
    all_reactions = data[data['rn'] == atlas_id]
    # Get only the forward
    f_react = all_reactions[all_reactions['direction'] == "forward"]
    # Get the reverse
    r_react = all_reactions[all_reactions['direction'] == "reverse"]
    if len(f_react) != 0:
        item = f_react.values.tolist()
    elif len(r_react) != 0:
        item = r_react.values.tolist()
    else:
        print(f"Skipping {atlas_id} due to no forward or reverse")
        continue
    # Get the reactants and products
    lhs = []
    rhs = []
    for j, entry in enumerate(item):
        if entry[4] < 0:
            lhs.append(f"{int(abs(entry[4]))} {entry[3]}")
        else:
            rhs.append(f"{int(abs(entry[4]))} {entry[3]}")
    lhs = " + ".join(lhs)
    rhs = " + ".join(rhs)
    eq = f"{lhs} <=> {rhs}"
    eq_list.append(eq)
    atlas_list.append(atlas_id)
    # if i > 100:
    #     break

# Write the data to a file
outfile = r"C:\Users\louie\skunkworks\data\atlas_reactions.csv.zip"
# make the dataframe from the eq_list and the atlas_ids
df = pd.DataFrame({'ID': atlas_list, 'Reaction': eq_list})
df.to_csv(outfile, compression='zip', encoding='utf-8')
