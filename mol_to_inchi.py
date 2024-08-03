import os
import os.path

from rdkit.Chem import AllChem as Chem
import rdkit.Chem

def file_list(mypath=None):
    """
    This function generates a list of all files in a specified directory.
    If no directory is specified, it defaults to the current working directory.

    Parameters:
    mypath (str, optional): The path to the directory. Defaults to None, which means the current working directory.

    Returns:
    list: A list of all files in the specified directory.
    """
    mypath = mypath or os.getcwd()  # If no path is provided, use the current working directory
    return [f for f in os.listdir(mypath) if
            os.path.isfile(os.path.join(mypath, f))]  # Return a list of all files in the directory


def file_list_all(mypath=None):
    """
    This function generates a list of all files in a specified directory and its subdirectories.
    If no directory is specified, it defaults to the current working directory.

    Parameters:
    mypath (str, optional): The path to the directory. Defaults to None, which means the current working directory.

    Returns:
    list: A list of all files in the specified directory and its subdirectories.
    """
    mypath = mypath or os.getcwd()  # If no path is provided, use the current working directory
    files = []
    # os.walk generates the file names in a directory tree by walking the tree either top-down or bottom-up
    for dirpath, dirnames, filenames in os.walk(mypath):
        for filename in filenames:
            # os.path.join joins one or more path components intelligently
            files.append(os.path.join(dirpath, filename))
    return files

def replace_r_group(target_file, new_file, replace_atom="C"):
    """
    This function replaces all instances of "R#" in a .mol file with a specified atom.

    Parameters:
    target_file (str): The path to the .mol file to be read.
    new_file (str): The path to the new .mol file to be written.
    replace_atom (str, optional): The atom to replace the R groups with. Defaults to "H".

    Returns:
    None
    """
    with open(target_file, "r") as f, open(new_file, "w") as nf:
        f = [line for line in f if "M  " not in line]
        f = [l.replace("R#", replace_atom + " ")
             .replace("R", replace_atom)
             .replace("*", replace_atom) for l in f]
        f = f[:-1]
        for line in f:
            nf.write(line)
        nf.write("M  END\n\n")


def check_for_X_group(target_file):
    with open(target_file, "r") as f:
        for line in f:
            if "X" in line:
                return True

def embed_mol(mol):
    # Sanitise the molecule
    Chem.SanitizeMol(mol, catchErrors=True)
    # Update the properties
    mol.UpdatePropertyCache()
    # kekulize the molecule
    Chem.Kekulize(mol)
    # Add hydrogens
    mol = rdkit.Chem.AddHs(mol)
    # Embed the molecule
    Chem.EmbedMolecule(mol, maxAttempts=100000)
    return mol


if __name__ == "__main__":
    print("Program started", flush=True)
    target = "C"
    target_dir = os.path.abspath(r"C:/Users/louie/skunkworks/data/kegg_data_" + target)
    files = file_list_all(target_dir)
    # filter the files to only include .mol files
    files = [f for f in files if "mod" not in f]

    for i, file in enumerate(files):
        if i < 406:
            continue
        print(f"Processing file {i + 1}/{len(files)}: {os.path.basename(file)}", flush=True)
        # Load the molecule using RDKit
        # mol = Chem.MolFromMolFile(file)

        # Check if the molecule has an X group
        if check_for_X_group(file):
            print("Skipping file for X", flush=True)
            continue

        # Make a temporary new file
        f_load_new = file.split(".")[0] + "_mod.mol"
        # Replace the R groups with Hydrogens
        replace_r_group(file, f_load_new)
        # Load the new molecule
        mol = Chem.MolFromMolFile(f_load_new, removeHs=False, sanitize=False)
        # Remove the temporary file
        os.remove(f_load_new)

        # mol = Chem.MolFromMolFile(file, removeHs=False, sanitize=False)

        # Embed the molecule in 3D
        mol = embed_mol(mol)
        # convert to inchi
        inchi = Chem.MolToInchi(mol)
        print(inchi)



    print("Program finished", flush=True)




