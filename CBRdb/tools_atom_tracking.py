import os
import re

import matplotlib.pyplot as plt
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rxnmapper import RXNMapper

script_dir = os.path.dirname(os.path.abspath(__file__))


def canonicalize_smiles(smiles):
    """
    Convert a SMILES string to its canonical form, removing all atom map numbers.
    Returns None if input is invalid or cannot be parsed.
    """
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def smiles_to_mols(smiles):
    """
    Convert a dot-separated SMILES string into a list of RDKit Mol objects.
    Only valid fragments are returned as Mol objects.
    """
    mols = [Chem.MolFromSmiles(smi) for smi in smiles.strip().split('.') if smi]
    return [mol for mol in mols if mol is not None]


def label_atom(mol, atom_idx, map_num):
    """
    Label a specific atom in a molecule with a given atom map number.
    Returns a new RDKit Mol object with the label applied.
    """
    mol = Chem.RWMol(mol)
    if 0 <= atom_idx < mol.GetNumAtoms():
        mol.GetAtomWithIdx(atom_idx).SetAtomMapNum(map_num)
    return mol


def trace_atom(atom_map_num, mols):
    """
    Search a list of molecules for an atom with a specific atom map number.
    Returns the first molecule containing the atom, or None if not found.
    """
    for mol in mols:
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() == atom_map_num:
                return mol
    return None


def remove_atom_map_numbers(mol):
    """
    Remove all atom map numbers from a molecule and return its canonical SMILES.
    Returns None if input is None.
    """
    if mol is None:
        return None
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def visualize_mapped_rxn(mapped_rxn_smiles, reaction_id, save_images=False):
    """
    Visualize the mapped reaction SMILES. If save_images is True, saves the image
    in a folder named 'Reactions Visualizations' in the same directory as the output CSV,
    using the reaction ID as the filename.
    """
    try:
        rxn = Chem.rdChemReactions.ReactionFromSmarts(mapped_rxn_smiles, useSmiles=True)
        img = Draw.ReactionToImage(rxn, subImgSize=(600, 600))
        vis_dir = os.path.join(os.path.dirname(output_file), "Reactions Visualizations")
        os.makedirs(vis_dir, exist_ok=True)
        img_path = os.path.join(vis_dir, f"{reaction_id}.png")
        plt.imshow(img)
        plt.axis('off')
        plt.title(f"Reaction ID: {reaction_id}")
        plt.tight_layout()
        plt.savefig(img_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"[INFO] Saved reaction visualization to {img_path}")
    except Exception as e:
        print(f"[ERROR] Could not visualize reaction {reaction_id}: {e}")


def find_file(filename):
    """
    Search for a file in the current script directory and up to two parent directories.
    Returns the absolute path if found, otherwise raises FileNotFoundError.
    """
    for base in [script_dir, os.path.join(script_dir, ".."), os.path.join(script_dir, "..", "..")]:
        candidate = os.path.abspath(os.path.join(base, filename))
        if os.path.exists(candidate):
            return candidate
    raise FileNotFoundError(f"{filename} not found in current or parent directories.")


def get_or_build_canonical_index(df_molecules, output_file):
    """
    Load a canonical index from CSV if it exists, otherwise build it from the molecule DataFrame.
    Each fragment of each molecule is canonicalized and indexed.
    Returns a DataFrame with columns: compound_id, original_smiles, fragment_smiles, canonical_smiles.
    """
    if os.path.exists(output_file):
        print(f"[INFO] Canonical index found at {output_file}, loading it.")
        return pd.read_csv(output_file)
    else:
        print(f"[INFO] Canonical index not found. Building it...")
        index = []
        for _, row in df_molecules.iterrows():
            cid = row['compound_id']
            smiles = row['smiles']
            if pd.isna(smiles):
                continue
            fragments = smiles.strip().split('.')
            for frag in fragments:
                mol = Chem.MolFromSmiles(frag)
                if mol:
                    canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
                    index.append({
                        'compound_id': cid,
                        'original_smiles': smiles,
                        'fragment_smiles': frag,
                        'canonical_smiles': canonical
                    })
                else:
                    print(f"[WARNING] Invalid fragment: {frag} for {cid}")
        df_index = pd.DataFrame(index)
        df_index.to_csv(output_file, index=False)
        print(f"[INFO] Canonical index saved to {output_file}")
        return df_index


def lookup_compound_id(canonical_smiles):
    """
    Look up the compound ID for a given canonical SMILES string in the canonical index DataFrame.
    Returns the compound ID if found, otherwise 'Unknown'.
    """
    match = df_canonical_index[df_canonical_index['canonical_smiles'] == canonical_smiles]
    if not match.empty:
        return match['compound_id'].values[0]
    return 'Unknown'


def omit_shared_compounds(sub_ids, prod_ids):
    shared = set(sub_ids) & set(prod_ids)
    if shared:
        print(f"[INFO] Omitting compounds present on both sides: {shared}")
    return [cid for cid in sub_ids if cid not in shared], [cid for cid in prod_ids if cid not in shared]


def remove_stoichiometry(reaction_str):
    # Remove numbers at the start of each compound (e.g., '2 H2O' -> 'H2O')
    def repl(match):
        return match.group(2)

    # Remove numbers and spaces before compound IDs
    return re.sub(r'(\b|\s)(\d+\s+)([A-Za-z0-9_]+)', r'\1\3', reaction_str)


def get_cofactor_filtered_reaction(substrate_ids, product_ids, cofactor_ids, remove_cofactors):
    # Remove cofactors from both sides if remove_cofactors is True
    if not remove_cofactors:
        return ''
    filtered_subs = [cid for cid in substrate_ids if cid not in cofactor_ids]
    filtered_prods = [cid for cid in product_ids if cid not in cofactor_ids]
    # If no change, return empty string
    if filtered_subs == substrate_ids and filtered_prods == product_ids:
        return 'No cofactor found'
    sub_names = [cid for cid in filtered_subs]
    prod_names = [cid for cid in filtered_prods]
    return ' + '.join(sub_names) + ' <=> ' + ' + '.join(prod_names)


def count_chiral_centers(mol):
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    return len(chiral_centers)


def selection_key_highest(i):
    ccc = count_chiral_centers(substrate_mols[i])
    heavy = substrate_mols[i].GetNumHeavyAtoms()
    # For tie-breaking: sort by chiral centers, then by heavy atom count (descending)
    return ccc, heavy


def selection_key_lowest(i):
    ccc = count_chiral_centers(substrate_mols[i])
    heavy = substrate_mols[i].GetNumHeavyAtoms()
    # For lowest: sort by chiral centers ascending, then by heavy atom count descending
    return ccc, -heavy
