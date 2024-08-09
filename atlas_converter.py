import os
import time
import numpy as np
import re

import pandas as pd

file_path = r"C:\Users\louie\skunkworks\data\atlas_reactions.csv"

# load the data using pandas
data = pd.read_csv(file_path)
print(data.head())
# get the 
atlas_ids = list(set(data["rn"]))

print(atlas_ids)