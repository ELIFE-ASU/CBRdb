import pytest
from rdkit import Chem

# Import functions to test from the main script.
# If the main script is not importable as a module, you may need to copy the functions here or refactor the code.

# Import functions to test from the main script.
from AtomTracking4CBRdb import (
    canonicalize_smiles,
    smiles_to_mols,
    label_atom,
    trace_atom,
    remove_atom_map_numbers,
)

def test_canonicalize_smiles(smiles, expected):
    """
    Test canonicalization of SMILES strings.
    Returns canonical SMILES or None for invalid input.
    """
    assert canonicalize_smiles(smiles) == expected

def test_smiles_to_mols_valid():
    """
    Test conversion of valid dot-separated SMILES to list of RDKit Mol objects.
    """
    smiles = "CCO.CN"
    mols = smiles_to_mols(smiles)
    assert len(mols) == 2
    assert all(isinstance(m, Chem.Mol) for m in mols)

def test_smiles_to_mols_invalid():
    """
    Test conversion of SMILES with one valid and one invalid fragment.
    Only valid fragments should be returned as Mol objects.
    """
    smiles = "CCO.invalid"
    mols = smiles_to_mols(smiles)
    assert len(mols) == 1
    assert Chem.MolToSmiles(mols[0]) == "CCO"

def test_label_atom_and_remove_atom_map_numbers():
    """
    Test labeling an atom with a map number and removing all atom map numbers.
    """
    mol = Chem.MolFromSmiles("CCO")
    labeled = label_atom(mol, 1, 42)
    # Atom 1 should have map number 42
    assert labeled.GetAtomWithIdx(1).GetAtomMapNum() == 42
    # Remove atom map numbers and check SMILES
    smiles = remove_atom_map_numbers(labeled)
    assert smiles == "CCO"

def test_trace_atom_found():
    """
    Test tracing an atom with a specific map number in a list of molecules.
    Should return the molecule containing the atom with the map number.
    """
    mol1 = Chem.MolFromSmiles("CCO")
    mol2 = Chem.MolFromSmiles("CCN")
    labeled = label_atom(mol2, 1, 99)
    mols = [mol1, labeled]
    found = trace_atom(99, mols)
    assert found is not None
    assert Chem.MolToSmiles(found) == Chem.MolToSmiles(labeled)

def test_trace_atom_not_found():
    """
    Test tracing an atom map number that does not exist in any molecule.
    Should return None.
    """
    mol1 = Chem.MolFromSmiles("CCO")
    mol2 = Chem.MolFromSmiles("CCN")
    mols = [mol1, mol2]
    found = trace_atom(123, mols)
    assert found is None

def test_remove_atom_map_numbers_none():
    """
    Test that remove_atom_map_numbers returns None when input is None.
    """
    assert remove_atom_map_numbers(None) is None

@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C(C(=O)O)N", "NCC(=O)O"),  # glycine
        ("[CH3:1][CH2:2][OH:3]", "CCO"),  # atom-mapped ethanol
        ("", None),
        (None, None),
        ("invalid", None),
    ],
)
def test_canonicalize_smiles(smiles, expected):
    """Test canonicalization of SMILES strings."""
    assert canonicalize_smiles(smiles) == expected