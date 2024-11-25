import os

import pandas as pd

from .tools_files import file_list, clean_empty_folders


def compounds_manual_add(molless_path='../data/C_IDs_bad.dat',
                         target_dir='../../data/kegg_data_C_full',
                         good_file="../data/C_IDs_good.dat"
                         ):
    # Set the absolute path
    molless_path = os.path.abspath(molless_path)
    target_dir = os.path.abspath(target_dir)
    good_file = os.path.abspath(good_file)

    data = pd.read_csv(molless_path)
    data = data["# Bad IDs"].tolist()
    print("data loaded", flush=True)
    print("data head", data[:4], flush=True)

    # Remove empty folders
    clean_empty_folders(target_dir)
    # List all files in the target directory
    files = file_list(target_dir)
    print(f"Files in the target directory: {files}", flush=True)
    r_list = []
    good_list = []
    # loop over the data
    for i, id in enumerate(data):
        # load the data
        f_path = os.path.abspath(f'{target_dir}/{id}/{id}.data')
        try:
            with open(f_path, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    line = line.strip()
                    # Check that the cid has an associated reaction
                    if line.startswith('REACTION'):
                        r_list.append(id)
        except FileNotFoundError:
            print(f"ID: {id} not found", flush=True)

    # Let us refine the list by removing the generic compounds
    for i, id in enumerate(r_list):
        good_flag = True
        # load the data
        f_path = os.path.abspath(f'{target_dir}/{id}/{id}.data')
        with open(f_path, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip().lower()
                # Check if there is a comment line
                if line.startswith('comment'):
                    # if "Generic compound" in line:
                    #     print(f"ID: {id} is a generic compound", flush=True)
                    #     good_flag = False
                    if "protein" in line:
                        print(f"ID: {id} is a protein", flush=True)
                        good_flag = False
                if line.startswith('name'):
                    # search if substring is in the line
                    if "glycan" in line:
                        print(f"ID: {id} is a glycan", flush=True)
                        good_flag = False
                    if "protein" in line:
                        print(f"ID: {id} is a protein", flush=True)
                        good_flag = False
                    if "peptide" in line:
                        print(f"ID: {id} is a peptide", flush=True)
                        good_flag = False
                    if "rna" in line:
                        print(f"ID: {id} is a RNA", flush=True)
                        good_flag = False
                    if "dna" in line:
                        print(f"ID: {id} is a DNA", flush=True)
                        good_flag = False
                    if "lase" in line:
                        print(f"ID: {id} is a enzyme", flush=True)
                        good_flag = False
                    if "steroid" in line:
                        print(f"ID: {id} is a steroid", flush=True)
                        good_flag = False
                    if "lipid" in line:
                        print(f"ID: {id} is a lipid", flush=True)
                        good_flag = False
                    if "lignin" in line:
                        print(f"ID: {id} is a lignin", flush=True)
                        good_flag = False

                if line.startswith("sequence"):
                    print(f"ID: {id} is a sequence", flush=True)
                    good_flag = False
            if good_flag:
                # append the good list
                good_list.append(id)

    # Now we have a list of all the reactions
    good_list = list(set(good_list))
    good_list.sort()
    print(f"Good list: {good_list}", flush=True)
    print(f"Good list length: {len(good_list)}", flush=True)
    # Save the good list
    with open(good_file, "w") as f:
        f.write("# Good compound IDs to follow up on\n")
        for id in good_list:
            f.write(f"{id}\n")
