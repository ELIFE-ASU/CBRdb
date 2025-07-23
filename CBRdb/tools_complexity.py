import traceback
from typing import Dict, Any, Optional

import networkx as nx
import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.GraphDescriptors import BertzCT
from rdkit.Chem.rdchem import Mol


def count_unique_bonds(mol: Mol) -> int:
    """
    Counts the number of unique bonds in a molecule.

    This function iterates over all the bonds in the given RDKit molecule object,
    identifies unique bonds based on the atom types (symbols) and bond type,
    and returns the count of these unique bonds.

    A bond is considered unique if the pair of atom types (sorted alphabetically)
    and the bond type are distinct.

    Parameters:
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The RDKit molecule object whose bonds are to be analyzed.

    Returns:
    -------
    int
        The number of unique bonds in the molecule.
    """
    unique_bonds = set()  # A set to store unique bonds
    for bond in mol.GetBonds():
        # Get the atom types (symbols) of the bonded atoms and sort them
        atom_types = tuple(sorted([bond.GetBeginAtom().GetSymbol(), bond.GetEndAtom().GetSymbol()]))
        # Get the bond type (e.g., single, double, etc.)
        bond_type = bond.GetBondType()
        # Add the unique bond (atom types and bond type) to the set
        unique_bonds.add((atom_types, bond_type))
    # Return the count of unique bonds
    return len(unique_bonds)


def molecular_weight(mol: Mol) -> float:
    """
    Calculates the molecular weight of a molecule.

    This function uses RDKit's molecular descriptors to compute the
    molecular weight of the given molecule.

    Parameters:
    -----------
    mol : rdkit.Chem.rdchem.Mol
        The RDKit molecule object for which the molecular weight is to be calculated.

    Returns:
    --------
    float
        The molecular weight of the molecule.
    """
    return Descriptors.MolWt(mol)


def bertz_complexity(mol: Mol) -> float:
    """
    Calculates the Bertz complexity of a molecule.

    The Bertz complexity is a molecular descriptor that quantifies the
    structural complexity of a molecule. It is based on graph theory
    and considers factors such as the number of atoms, bonds, and
    branching in the molecular structure.

    Parameters:
    -----------
    mol : rdkit.Chem.rdchem.Mol
        The RDKit molecule object for which the Bertz complexity is to be calculated.

    Returns:
    --------
    float
        The Bertz complexity of the molecule.
    """
    return BertzCT(mol)


def wiener_index(mol: Mol) -> int:
    """
    Calculates the Wiener index of a molecule.

    The Wiener index is a topological descriptor that represents the sum of
    all shortest path distances between pairs of vertices in a molecular graph.
    It is used in cheminformatics to study molecular structure-activity
    relationships and predict physicochemical properties.

    Parameters:
    -----------
    mol : rdkit.Chem.rdchem.Mol
        The RDKit molecule object for which the Wiener index is to be calculated.

    Returns:
    --------
    int
        The Wiener index of the molecule.
    """
    distance_matrix = Chem.rdmolops.GetDistanceMatrix(mol)
    graph = nx.Graph()

    # Add nodes to the graph
    for i in range(len(distance_matrix)):
        graph.add_node(i)

    # Add edges with weights based on the distance matrix
    for i in range(len(distance_matrix)):
        for j in range(i + 1, len(distance_matrix)):
            graph.add_edge(i, j, weight=distance_matrix[i, j])

    return nx.wiener_index(graph)


def balaban_index(mol: Mol) -> float:
    """
    Calculates the Balaban index of a molecule.

    The Balaban index is a topological descriptor that provides a measure of
    molecular connectivity. It is used in cheminformatics to study molecular
    structure-activity relationships and predict physicochemical properties.

    Parameters:
    -----------
    mol : rdkit.Chem.rdchem.Mol
        The RDKit molecule object for which the Balaban index is to be calculated.

    Returns:
    --------
    float
        The Balaban index of the molecule.
    """
    return Descriptors.BalabanJ(mol)


def randic_index(mol: Mol) -> float:
    """
    Calculates the Randic index of a molecule.

    Randic index is a topological descriptor calculated by summing the
    inverse square roots of the product of the degrees of connected atom pairs.
    It is used for the study of molecular structure-activity
    relationships and the prediction of physicochemical properties.

    :param mol: An RDKit molecule object.
    :return: The Randic index.
    """
    adj_matrix = Chem.rdmolops.GetAdjacencyMatrix(mol)
    degrees = [sum(row) for row in adj_matrix]
    randic_sum = 0
    for i, row in enumerate(adj_matrix):
        for j, val in enumerate(row):
            if val == 1:
                randic_sum += 1 / (degrees[i] * degrees[j]) ** 0.5
    return randic_sum / 2


def kirchhoff_index(mol: Mol) -> float:
    """
    Calculates the Kirchhoff index of a molecule.

    Kirchhoff index is a topological index calculated as the sum of the effective
    resistances between all pairs of vertices in the molecular graph. It is used for
    predicting physicochemical properties and molecular activities.

    :param mol: An RDKit molecule object.
    :return: The Kirchhoff index.
    """
    adjacency_matrix = Chem.rdmolops.GetAdjacencyMatrix(mol).astype(np.float64)
    degree_matrix = np.diag(np.sum(adjacency_matrix, axis=1))
    laplacian_matrix = degree_matrix - adjacency_matrix
    pseudo_inverse_laplacian = np.linalg.pinv(laplacian_matrix)
    diagonal_elements = np.diagonal(pseudo_inverse_laplacian)

    kirchhoff_sum = 0
    for i in range(len(diagonal_elements)):
        for j in range(i + 1, len(diagonal_elements)):
            kirchhoff_sum += (diagonal_elements[i] + diagonal_elements[j] - 2 * pseudo_inverse_laplacian[i, j])

    return kirchhoff_sum


def spacial_score(mol: Mol, normalise: bool = False) -> float:
    """
    Calculates the spacial score of a molecule. https://github.com/frog2000/Spacial-Score

    Spacial score is a descriptor that quantifies the spatial arrangement
    of atoms in a molecule. It can be used to predict various molecular properties.

    :param mol: An RDKit molecule object.
    :param normalise: A boolean indicating whether to normalise the score.
    :return: The spacial score of the molecule.
    """
    return rdkit.Chem.SpacialScore.SPS(mol, normalise)


def get_mol_descriptors(mol: Mol, missing_val: Optional[Any] = None) -> Dict[str, Any]:
    """
    Calculates molecular descriptors for a given molecule. Please note that there are a lot of descriptors.

    https://greglandrum.github.io/rdkit-blog/posts/2022-12-23-descriptor-tutorial.html

    This function iterates over all available molecular descriptors in RDKit,
    calculates each descriptor for the provided molecule, and stores the results
    in a dictionary. If a descriptor calculation fails, a specified missing value
    is assigned.

    :param mol: An RDKit molecule object.
    :param missing_val: The value to assign if a descriptor calculation fails. Default is None.
    :return: A dictionary with descriptor names as keys and their calculated values as values.
    """
    res = {}
    for nm, fn in Descriptors._descList:
        try:
            res[nm] = fn(mol)
        except:
            traceback.print_exc()
            res[nm] = missing_val
    return res


def get_chirality(mol: Mol) -> int:
    """
    Determine the chirality of a molecule.

    This function calculates the number of chiral centres in a given RDKit molecule object.

    Parameters:
        mol (rdkit.Chem.rdchem.Mol): An RDKit molecule object.

    Returns:
        int: The number of chiral centres in the molecule.
    """
    nc = len(Chem.FindMolChiralCenters(mol,
                                       useLegacyImplementation=False,
                                       includeUnassigned=True,
                                       includeCIP=False))
    return nc
