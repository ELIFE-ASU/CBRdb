from rdkit import Chem as Chem
from rdkit import RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)


def standardize_mol(mol):
    # Standardize the molecule
    mol.UpdatePropertyCache(strict=False)
    Chem.SetConjugation(mol)
    Chem.SetHybridization(mol)
    # Normalize the molecule
    Chem.SanitizeMol(mol, sanitizeOps=(Chem.SANITIZE_ALL ^ Chem.SANITIZE_CLEANUP ^ Chem.SANITIZE_PROPERTIES))
    rdMolStandardize.NormalizeInPlace(mol)
    return mol


def get_chirality(mol):
    try:
        return len(
            Chem.FindMolChiralCenters(mol,
                                      useLegacyImplementation=False,
                                      includeUnassigned=True,
                                      includeCIP=False))
    except:
        return float('NaN')
