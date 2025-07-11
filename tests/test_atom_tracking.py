from rdkit import Chem

from CBRdb import (
    canonicalize_smiles,
    smiles_to_mols,
    label_atom,
    trace_atom,
    remove_atom_map_numbers,
)


def test_canonicalize_smiles_valid():
    """
    Test that the `canonicalize_smiles` function correctly canonicalizes a SMILES string
    and removes atom map numbers. Verifies that the output matches the canonical form
    of equivalent SMILES strings without atom mapping.
    """
    # Should canonicalize and remove atom map numbers
    smi = "[C:1](O)C"
    result = canonicalize_smiles(smi)
    assert result == canonicalize_smiles("COC") or result == canonicalize_smiles("CCO")


def test_canonicalize_smiles_invalid():
    """
    Test the `canonicalize_smiles` function with invalid input values.

    This test verifies that the function returns `None` when provided with:
    - An empty string as input.
    - `None` as input.
    - An invalid SMILES string.

    Expected behavior: `canonicalize_smiles` should return `None` for all invalid inputs.
    """
    assert canonicalize_smiles("") is None
    assert canonicalize_smiles(None) is None
    assert canonicalize_smiles("not_a_smiles") is None


def test_smiles_to_mols_single():
    """
    Test the `smiles_to_mols` function with a single SMILES string.

    This test verifies that:
    - Passing a single SMILES string ("CCO") to `smiles_to_mols` returns a list containing one molecule.
    - The returned object in the list is an instance of `Chem.Mol`.
    """
    mols = smiles_to_mols("CCO")
    assert len(mols) == 1
    assert isinstance(mols[0], Chem.Mol)


def test_smiles_to_mols_multiple():
    """
    Test that the `smiles_to_mols` function correctly parses a SMILES string containing multiple molecules.

    This test verifies that:
    - The function splits the input SMILES string "CCO.CN" into two separate molecules.
    - The returned list contains exactly two elements.
    - Each element in the list is an instance of `Chem.Mol`.
    """
    mols = smiles_to_mols("CCO.CN")
    assert len(mols) == 2
    assert all(isinstance(m, Chem.Mol) for m in mols)


def test_smiles_to_mols_invalid_fragment():
    """
    Test that `smiles_to_mols` correctly handles an invalid SMILES fragment.

    This test verifies that when an invalid fragment is included in the input string,
    the function returns only the valid molecules, ignoring the invalid part.

    Assertions:
        - The returned list contains only one molecule.
        - The SMILES string of the returned molecule matches the valid input fragment ("CCO").
    """
    mols = smiles_to_mols("CCO.not_a_smiles")
    assert len(mols) == 1
    assert Chem.MolToSmiles(mols[0]) == Chem.MolToSmiles(Chem.MolFromSmiles("CCO"))


def test_label_atom_and_trace_atom():
    """
    Test the functionality of labeling an atom in a molecule and tracing it by its atom map number.

    This test performs the following steps:
    1. Creates an ethanol molecule from its SMILES representation.
    2. Labels the atom at index 1 with an atom map number 42 using the `label_atom` function.
    3. Asserts that the labeled atom with atom map number 42 exists in the molecule.
    4. Uses the `trace_atom` function to find the atom with atom map number 42 in a list containing the labeled molecule.
    5. Asserts that the result is not None and that the atom with atom map number 42 is present in the traced molecule.

    Functions tested:
    - `label_atom(mol, atom_idx, atom_map_num)`: Labels a specific atom in a molecule with a given atom map number.
    - `trace_atom(atom_map_num, mols)`: Finds and returns the molecule containing the atom with the specified atom map number.

    Assertions:
    - The atom with the specified atom map number is present after labeling.
    - The atom can be traced successfully in the list of molecules.
    """
    mol = Chem.MolFromSmiles("CCO")
    labeled = label_atom(mol, 1, 42)
    found = False
    for atom in labeled.GetAtoms():
        if atom.GetAtomMapNum() == 42:
            found = True
    assert found
    # trace_atom should find the labeled atom
    mols = [labeled]
    result = trace_atom(42, mols)
    assert result is not None
    assert any(atom.GetAtomMapNum() == 42 for atom in result.GetAtoms())


def test_trace_atom_not_found():
    """
    Test that trace_atom returns None when the specified atom index does not exist in the molecule.

    This test creates a simple ethanol molecule (CCO), places it in a list, and attempts to trace an atom
    with index 99, which is out of bounds. The expected behavior is that trace_atom returns None.
    """
    mol = Chem.MolFromSmiles("CCO")
    mols = [mol]
    assert trace_atom(99, mols) is None


def test_remove_atom_map_numbers():
    """
    Test the remove_atom_map_numbers function to ensure that it removes atom map numbers from a molecule.
    Verifies that the resulting SMILES string does not contain atom map numbers and matches the canonical SMILES for ethanol ("CCO").
    """
    mol = Chem.MolFromSmiles("[CH3:5][CH2:7][OH:9]")
    result = remove_atom_map_numbers(mol)
    # Should be canonical and have no atom map numbers
    assert ":" not in result
    assert result == canonicalize_smiles("CCO")


def test_remove_atom_map_numbers_none():
    """
    Test that remove_atom_map_numbers returns None when given None as input.
    """
    assert remove_atom_map_numbers(None) is None
