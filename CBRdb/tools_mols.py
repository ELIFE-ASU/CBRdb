from rdkit import Chem as Chem
from rdkit import RDLogger
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

def sanitize_mol(mol):
    """
    Standardizes and normalizes the given molecule.

    Parameters:
    mol (rdkit.Chem.Mol): The molecule to be sanitized.

    Returns:
    rdkit.Chem.Mol: The sanitized molecule.
    """
    # Standardize the molecule
    mol.UpdatePropertyCache(strict=False)
    Chem.SetConjugation(mol)
    Chem.SetHybridization(mol)
    # Normalize the molecule
    Chem.SanitizeMol(mol, sanitizeOps=(Chem.SANITIZE_ALL ^ Chem.SANITIZE_CLEANUP ^ Chem.SANITIZE_PROPERTIES))
    # Update the properties
    mol.UpdatePropertyCache(strict=False)
    return mol


def standardize_mol(mol):
    """
    Standardizes the given molecule by sanitizing, normalizing, and kekulizing it.

    Parameters:
    mol (rdkit.Chem.Mol): The molecule to be standardized.

    Returns:
    rdkit.Chem.Mol: The standardized molecule.
    """
    # Sanitize the molecule
    mol = sanitize_mol(mol)
    # Normalize the molecule
    rdMolStandardize.NormalizeInPlace(mol)
    # Kekulize the molecule
    Chem.Kekulize(mol)
    # Update the properties
    mol.UpdatePropertyCache(strict=False)
    return mol


def fix_r_group(mol, target="[H:0]", re="*"):
    """
    Replaces all occurrences of a specified target atom or group in the molecule with another specified atom or group.

    Parameters:
    mol (rdkit.Chem.Mol): The molecule to be modified.
    target (str): The atom or group to be replaced (default is "[H:0]").
    re (str): The replacement atom or group (default is "*").

    Returns:
    rdkit.Chem.Mol: The standardized molecule after replacement.
    """
    # Convert the molecule to a SMILES string
    smi = Chem.MolToSmiles(mol)
    # Replace all occurrences of the target with the replacement
    smi = smi.replace(target, re)
    # Convert the SMILES string back to a molecule without sanitizing
    mol = Chem.MolFromSmiles(smi, sanitize=False)
    # Standardize the molecule
    return standardize_mol(mol)


def get_chirality(mol):
    """
    Determines the number of chiral centers in the given molecule.

    Parameters:
    mol (rdkit.Chem.Mol): The molecule to be analyzed.

    Returns:
    int: The number of chiral centers in the molecule.
    float: NaN if an error occurs during the analysis.
    """
    try:
        return len(
            Chem.FindMolChiralCenters(mol,
                                      useLegacyImplementation=False,
                                      includeUnassigned=True,
                                      includeCIP=False))
    except:
        return float('NaN')


def mol_replacer(smi, target="[H]"):
    """
    Replaces all occurrences of the '*' atom in the given SMILES string with the specified target atom or group.

    Parameters:
    smi (str): The SMILES string representing the molecule.
    target (str): The atom or group to replace '*' with (default is "[H]").

    Returns:
    str: The modified SMILES string with all '*' atoms replaced by the target.
    """
    # Convert the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smi)
    # Create a molecule for the '*' atom
    star = Chem.MolFromSmiles('*')
    # Create a molecule for the target atom or group
    tar = Chem.MolFromSmiles(target)
    # Replace all occurrences of '*' with the target in the molecule
    mol_re = Chem.ReplaceSubstructs(mol, star, tar, replaceAll=True)
    # Convert the modified molecule back to a SMILES string and return it
    return Chem.MolToSmiles(mol_re[0])


def get_mol_descriptors(mol):
    """
    Calculates various molecular descriptors for the given molecule.

    Parameters:
    mol (rdkit.Chem.Mol): The molecule to be analyzed.

    Returns:
    tuple: A tuple containing the molecular formula (str), exact molecular weight (float),
           number of heavy atoms (int), and number of chiral centers (int or float).
    """
    mol = sanitize_mol(mol)
    return (rdMolDescriptors.CalcMolFormula(mol),
            rdMolDescriptors.CalcExactMolWt(mol),
            rdMolDescriptors.CalcNumHeavyAtoms(mol),
            get_chirality(mol))
