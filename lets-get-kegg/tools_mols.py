from rdkit import Chem as Chem
from rdkit import RDLogger
from rdkit.Chem import rdMolDescriptors
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
    # kekulize the molecule
    Chem.Kekulize(mol)
    # Update the properties
    mol.UpdatePropertyCache(strict=False)
    return mol


def fix_r_group(smi, target="[H:0]", re="*"):
    # replace all occurrences of [H:0] with *
    smi = smi.replace(target, re)
    # Get the molecule
    mol = Chem.MolFromSmiles(smi, sanitize=False)
    # Standardize the molecule
    mol = standardize_mol(mol)
    # Return the smiles
    return Chem.MolToSmiles(mol)


def get_chirality(mol):
    try:
        return len(
            Chem.FindMolChiralCenters(mol,
                                      useLegacyImplementation=False,
                                      includeUnassigned=True,
                                      includeCIP=False))
    except:
        return float('NaN')


def mol_replacer(smi):
    mol = Chem.MolFromSmiles(smi)
    star = Chem.MolFromSmiles('*')
    h = Chem.MolFromSmiles('[H]')
    mol_re = Chem.ReplaceSubstructs(mol, star, h, replaceAll=True)
    return Chem.MolToSmiles(mol_re[0])


def get_mol_descriptors(mol):
    return (rdMolDescriptors.CalcMolFormula(mol),
            rdMolDescriptors.CalcExactMolWt(mol),
            rdMolDescriptors.CalcNumRings(mol),
            get_chirality(mol))
