import math
import traceback
from collections import defaultdict
from typing import Dict, Any, Optional

import networkx as nx
import numpy as np
import rdkit
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.GraphDescriptors import BertzCT
from rdkit.Chem.SpacialScore import SPS
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


def bertz(mol: Mol) -> float:
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


def get_mol_descriptors(mol: Mol, missingval: Optional[Any] = None) -> Dict[str, Any]:
    """
    Calculates molecular descriptors for a given molecule. Please note that there are a lot of descriptors.

    https://greglandrum.github.io/rdkit-blog/posts/2022-12-23-descriptor-tutorial.html

    This function iterates over all available molecular descriptors in RDKit,
    calculates each descriptor for the provided molecule, and stores the results
    in a dictionary. If a descriptor calculation fails, a specified missing value
    is assigned.

    :param mol: An RDKit molecule object.
    :param missingval: The value to assign if a descriptor calculation fails. Default is None.
    :return: A dictionary with descriptor names as keys and their calculated values as values.
    """
    res = {}
    for nm, fn in Descriptors._descList:
        try:
            res[nm] = fn(mol)
        except:
            traceback.print_exc()
            res[nm] = missingval
    return res


def tanimoto_similarity(mol1: Mol, mol2: Mol) -> float:
    """
    Calculates the Tanimoto similarity between two molecules.

    Tanimoto similarity is a measure of the similarity between two sets of
    molecular fingerprints. It is commonly used in cheminformatics to compare
    the structural similarity of molecules.

    :param mol1: An RDKit molecule object representing the first molecule.
    :param mol2: An RDKit molecule object representing the second molecule.
    :return: The Tanimoto similarity between the two molecules.
    """
    fpgen = Chem.GetRDKitFPGenerator()
    fp1 = fpgen.GetFingerprint(mol1)
    fp2 = fpgen.GetFingerprint(mol2)
    return DataStructs.TanimotoSimilarity(fp1, fp2)


def dice_morgan_similarity(mol1: Mol, mol2: Mol, radius: int = 3) -> float:
    """
    Calculates the Dice similarity between two molecules using Morgan fingerprints.

    Dice similarity is a measure of the similarity between two sets of
    molecular fingerprints. It is commonly used in cheminformatics to compare
    the structural similarity of molecules.

    :param mol1: An RDKit molecule object representing the first molecule.
    :param mol2: An RDKit molecule object representing the second molecule.
    :param radius: The radius parameter for the Morgan fingerprint. Default is 3.
    :return: The Dice similarity between the two molecules.
    """
    fpgen = Chem.GetMorganGenerator(radius=radius)
    fp1 = fpgen.GetSparseCountFingerprint(mol1)
    fp2 = fpgen.GetSparseCountFingerprint(mol2)
    return DataStructs.DiceSimilarity(fp1, fp2)


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


def fcfp4(mol: Mol) -> int:
    """
    Generates the FCFP_4 fingerprint (functional-based ECFP4) for a molecule.

    https://doi.org/10.1021/ci0503558

    This function computes the FCFP_4 fingerprint of a molecule using RDKit's
    Morgan fingerprinting method. The fingerprint is generated with a radius of 2
    and 2048 bits, focusing on functional groups by setting `useFeatures=True`.

    Parameters:
    -----------
    mol : rdkit.Chem.rdchem.Mol
        The RDKit molecule object for which the FCFP_4 fingerprint is to be generated.

    Returns:
    --------
    int
        The number of bits set to 1 in the generated FCFP_4 fingerprint.
    """
    # Generate the FCFP_4 fingerprint with functional group focus
    fp = Chem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048, useFeatures=True)
    # Return the number of bits set to 1
    return fp.GetNumOnBits()


def _determine_atom_substituents(atom_id, mol, distance_matrix):
    """
    Determines the substituents of an atom in a molecule.

    This function identifies the substituents (neighboring atoms and their shells)
    for a given atom in a molecule based on the distance matrix. It also tracks
    shared neighbors and the maximum shell distance for each substituent.

    Parameters:
    ----------
    atom_id : int
        The index of the atom in the molecule for which substituents are to be determined.
    mol : rdkit.Chem.rdchem.Mol
        The RDKit molecule object containing the atom.
    distance_matrix : np.ndarray
        The distance matrix of the molecule, where each element represents the
        shortest path distance between two atoms.

    Returns:
    -------
    tuple
        A tuple containing:
        - subs (defaultdict): A dictionary mapping substituent indices to lists of atom indices.
        - shared_neighbors (defaultdict): A dictionary tracking how many substituents each atom is involved in.
        - max_shell (defaultdict): A dictionary mapping substituent indices to their maximum shell distance.
    """
    atom_paths = distance_matrix[atom_id]
    # Determine the direct neighbors of the atom
    neighbors = [n for n, i in enumerate(atom_paths) if i == 1]
    # Store the ids of the neighbors (substituents)
    subs = defaultdict(list)
    # Track in how many substituents an atom is involved (can happen in rings)
    shared_neighbors = defaultdict(int)
    # Determine the max path length for each substituent
    max_shell = defaultdict(int)
    for n in neighbors:
        subs[n].append(n)
        shared_neighbors[n] += 1
        max_shell[n] = 0
    # Second shell of neighbors
    min_dist = 2
    # Max distance from atom
    max_dist = int(np.max(atom_paths))
    for d in range(min_dist, max_dist + 1):
        new_shell = [n for n, i in enumerate(atom_paths) if i == d]
        for a_idx in new_shell:
            atom = mol.GetAtomWithIdx(a_idx)
            # Find neighbors of the current atom that are part of the substituent already
            for n in atom.GetNeighbors():
                n_idx = n.GetIdx()
                for k, v in subs.items():
                    # Check if the neighbor is in the substituent, not in the same shell,
                    # and the current atom hasn't been added yet
                    if n_idx in v and n_idx not in new_shell and a_idx not in v:
                        subs[k].append(a_idx)
                        shared_neighbors[a_idx] += 1
                        max_shell[k] = d

    return subs, shared_neighbors, max_shell


def _get_chemical_non_equivs(atom: Chem.rdchem.Atom, mol: Mol) -> float:
    """
    Calculates the chemical non-equivalence of an atom in a molecule.

    This function determines the number of unique substituent groups attached to an atom
    in a molecule. It uses the distance matrix of the molecule and the atom's substituents
    to identify unique groups based on their atomic symbols.

    Parameters:
    ----------
    atom : Chem.rdchem.Atom
        The RDKit atom object for which chemical non-equivalence is to be calculated.
    mol : Mol
        The RDKit molecule object containing the atom.

    Returns:
    -------
    float
        The number of unique substituent groups attached to the atom.
    """
    # Initialize a list to store substituents for up to 4 groups
    substituents = [[] for _ in range(4)]

    # Get the distance matrix of the molecule
    distance_matrix = Chem.GetDistanceMatrix(mol)

    # Determine the substituents of the atom
    atom_substituents = _determine_atom_substituents(atom.GetIdx(), mol, distance_matrix)[0]
    try:
        # Populate the substituents list with atomic symbols of neighboring atoms
        for item, key in enumerate(atom_substituents):
            for subatom in atom_substituents[key]:
                substituents[item].append(mol.GetAtomWithIdx(subatom).GetSymbol())

        # Calculate the number of unique substituent groups
        return float(len(set(tuple(sub) for sub in substituents if sub)))
    except Exception as e:
        print(f"Error calculating chemical non-equivalence for atom {atom.GetIdx()}: {e}")
        print(traceback.format_exc())
        return 0.0


def _get_bottcher_local_diversity(atom: Chem.rdchem.Atom) -> float:
    """
    Calculates the Bottcher local diversity of an atom.

    This function determines the diversity of an atom based on the unique
    neighboring atom types. It adds an additional value of 1.0 if the atom's
    own type is not present among its neighbors.

    Parameters:
    ----------
    atom : Chem.rdchem.Atom
        The RDKit atom object for which the Bottcher local diversity is to be calculated.

    Returns:
    -------
    float
        The Bottcher local diversity score of the atom.
    """
    # Get the set of unique symbols of neighboring atoms
    neighbors = {neighbor.GetSymbol() for neighbor in atom.GetNeighbors()}
    # Calculate diversity, adding 1.0 if the atom's symbol is not in its neighbors
    return len(neighbors) + (1.0 if atom.GetSymbol() not in neighbors else 0.0)


def _get_num_isomeric_possibilities(atom: Chem.rdchem.Atom) -> float:
    """
    Determines the number of isomeric possibilities for an atom.

    This function checks if the atom has a '_CIPCode' property, which indicates
    the presence of stereochemical information. If the property exists, the atom
    has two isomeric possibilities (e.g., R/S or E/Z). Otherwise, it has only one.

    Parameters:
    ----------
    atom : Chem.rdchem.Atom
        The RDKit atom object for which the number of isomeric possibilities is to be determined.

    Returns:
    -------
    float
        The number of isomeric possibilities: 2.0 if the '_CIPCode' property exists, otherwise 1.0.
    """
    return 2.0 if atom.HasProp('_CIPCode') else 1.0


def _get_num_valence_electrons(atom: Chem.rdchem.Atom, pt: Chem.rdchem.PeriodicTable) -> float:
    """
    Calculates the number of valence electrons for a given atom.

    This function uses the periodic table to determine the number of outer-shell
    (valence) electrons for the specified atom based on its atomic number.

    Parameters:
    ----------
    atom : Chem.rdchem.Atom
        The RDKit atom object for which the number of valence electrons is to be calculated.
    pt : Chem.rdchem.PeriodicTable
        The RDKit periodic table object used to retrieve atomic properties.

    Returns:
    -------
    float
        The number of valence electrons for the given atom.
    """
    return float(pt.GetNOuterElecs(pt.GetAtomicNumber(atom.GetSymbol())))


def _get_bottcher_bond_index(atom: Chem.rdchem.Atom) -> float:
    """
    Calculates the Bottcher bond index for a given atom.

    This function computes a ranking value based on the bond types connected to the atom.
    Each bond type is assigned a specific weight, and additional adjustments are made
    for aromatic bonds involving carbon or nitrogen atoms.

    Parameters:
    ----------
    atom : Chem.rdchem.Atom
        The RDKit atom object for which the Bottcher bond index is to be calculated.

    Returns:
    -------
    float
        The Bottcher bond index for the given atom.

    Raises:
    ------
    ValueError
        If an unsupported bond type is encountered.
    """
    b_sub_i_ranking = 0.0
    bond_weights = {
        'SINGLE': 1.0,
        'DOUBLE': 2.0,
        'TRIPLE': 3.0,
        'QUADRUPLE': 4.0,
        'QUINTUPLE': 5.0,
        'HEXTUPLE': 6.0
    }
    bonds = [str(bond.GetBondType()) for bond in atom.GetBonds()]
    for bond in bonds:
        b_sub_i_ranking += bond_weights.get(bond, 0.0)
        if bond not in bond_weights and bond != 'AROMATIC':
            raise ValueError(f"Unsupported bond type {bond}")

    if 'AROMATIC' in bonds:
        if atom.GetSymbol() == 'C':
            b_sub_i_ranking += 3.0
        elif atom.GetSymbol() == 'N':
            b_sub_i_ranking += 2.0
    return b_sub_i_ranking


def bottcher(mol: Mol) -> float:
    """
    Calculates the Bottcher complexity of a molecule.

    https://github.com/boskovicgroup/bottchercomplexity
    https://doi.org/10.1021/acs.jcim.5b00723

    The Bottcher complexity is a molecular descriptor that quantifies the structural
    complexity of a molecule. It considers factors such as chemical non-equivalence,
    local diversity, isomeric possibilities, valence electrons, and bond indices
    for each atom in the molecule.

    Parameters:
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The RDKit molecule object for which the Bottcher complexity is to be calculated.

    Returns:
    -------
    float
        The Bottcher complexity of the molecule.
    """
    complexity = 0.0
    # Assign stereochemistry to the molecule
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)
    pt = Chem.GetPeriodicTable()

    # Filter atoms to correct for symmetry
    atoms_corrected_for_symmetry = []
    atom_stereo_classes = set()
    for atom in mol.GetAtoms():
        cip_rank = atom.GetProp('_CIPRank')
        if cip_rank not in atom_stereo_classes:
            atoms_corrected_for_symmetry.append(atom)
            atom_stereo_classes.add(cip_rank)

    # Calculate complexity
    for atom in atoms_corrected_for_symmetry:
        d = _get_chemical_non_equivs(atom, mol)  # Chemical non-equivalence
        e = _get_bottcher_local_diversity(atom)  # Local diversity
        s = _get_num_isomeric_possibilities(atom)  # Isomeric possibilities
        v = _get_num_valence_electrons(atom, pt)  # Number of valence electrons
        b = _get_bottcher_bond_index(atom)  # Bond index
        # Update complexity using the calculated factors
        complexity += d * e * s * math.log(v * b, 2)

    return complexity


def proudfoot(mol: Mol) -> float:
    """
    Calculates the Proudfoot complexity of a molecule.

    https://doi.org/10.1016/j.bmcl.2017.03.008

    The Proudfoot complexity is a molecular descriptor that quantifies the structural
    complexity of a molecule. It is based on the distribution of molecular paths,
    atomic complexity, molecular complexity, log-sum complexity, and structural entropy.

    Parameters:
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The RDKit molecule object for which the Proudfoot complexity is to be calculated.

    Returns:
    -------
    float
        The Proudfoot complexity of the molecule.
    """
    # Generate the Morgan fingerprint for the molecule with a radius of 2
    fingerprint = rdMolDescriptors.GetMorganFingerprint(mol, 2)
    paths = fingerprint.GetNonzeroElements()

    # Calculate path frequencies per atom environment
    atom_paths = defaultdict(list)
    for path, count in paths.items():
        # Determine the atoms involved in the path
        atoms_in_path = path % mol.GetNumAtoms()
        atom_paths[atoms_in_path].append(count)

    # Step 1: Calculate atomic complexity (C_A)
    c_a_values = {}
    for atom, path_counts in atom_paths.items():
        total_paths = sum(path_counts)
        # Calculate the fraction of each path
        path_fractions = [count / total_paths for count in path_counts]
        # Compute atomic complexity using Shannon entropy
        ca = -sum(p * math.log2(p) for p in path_fractions) + math.log2(total_paths)
        c_a_values[atom] = ca

    # Step 2: Calculate molecular complexity (C_M)
    c_m = sum(c_a_values.values())

    # Step 3: Calculate log-sum complexity (C_M*)
    c_m_star = math.log2(sum(2 ** ca for ca in c_a_values.values()))

    # Step 4: Calculate structural entropy complexity (C_SE)
    atom_types = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    total_atoms = len(atom_types)
    # Calculate the frequency of each atom type
    type_frequencies = {atype: atom_types.count(atype) / total_atoms for atype in set(atom_types)}
    # Compute structural entropy using Shannon entropy
    c_se = -sum(freq * math.log2(freq) for freq in type_frequencies.values())

    return c_m


def mc1(mol: Mol) -> float:
    """
    Calculates the molecular connectivity index (MC1) of a molecule.

    https://pubs.acs.org/doi/full/10.1021/acs.jcim.5c00334

    The MC1 index is a measure of the proportion of non-divalent nodes
    in a molecule. It is calculated as 1 minus the ratio of divalent nodes
    to the total number of atoms in the molecule.

    Parameters:
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The RDKit molecule object for which the MC1 index is to be calculated.

    Returns:
    -------
    float
        The MC1 index of the molecule.
    """
    total_atoms = len(mol.GetAtoms())
    divalent_nodes = sum(1 for atom in mol.GetAtoms() if atom.GetDegree() == 2)
    return 1.0 - (divalent_nodes / total_atoms)


def mc2(mol: Mol) -> int:
    """
    Calculates a molecular connectivity index (MC2) for a molecule.

    https://pubs.acs.org/doi/full/10.1021/acs.jcim.5c00334

    This function identifies carbon-oxygen (C=O) double bonds in the molecule
    and checks if the carbon atom in the bond is connected to a nitrogen (N) or
    oxygen (O) atom (excluding the double-bonded oxygen). It then counts the
    number of non-divalent atoms that are not part of these specific C=O-X
    double bonds.

    Parameters:
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The RDKit molecule object for which the MC2 index is to be calculated.

    Returns:
    -------
    int
        The count of non-divalent atoms not involved in C=O-X double bonds.
    """
    double_bond_set = set()

    for bond in mol.GetBonds():
        # Check for a C=O double bond
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            atoms = [bond.GetBeginAtom(), bond.GetEndAtom()]
            if {atom.GetAtomicNum() for atom in atoms} == {6, 8}:  # C and O
                carbon = next(atom for atom in atoms if atom.GetAtomicNum() == 6)
                oxygen = next(atom for atom in atoms if atom.GetAtomicNum() == 8)

                # Check carbon's neighbors for N or O (excluding the double-bonded O)
                if any(neighbor.GetIdx() != oxygen.GetIdx() and neighbor.GetAtomicNum() in [7, 8]
                       for neighbor in carbon.GetNeighbors()):
                    double_bond_set.update([carbon.GetIdx(), oxygen.GetIdx()])

    # Count non-divalent atoms not in C=O-X double bonds
    return sum(1 for atom in mol.GetAtoms() if atom.GetDegree() != 2 and atom.GetIdx() not in double_bond_set)


def get_all_mol_descriptors(mol):
    """
    Calculates a comprehensive set of molecular descriptors for a given molecule.

    This function computes various molecular descriptors using RDKit and custom functions.
    It combines descriptors such as molecular formula, molecular weight, number of heavy atoms,
    and several topological indices into a single dictionary.

    Parameters:
    -----------
    mol : rdkit.Chem.rdchem.Mol
        The RDKit molecule object for which descriptors are to be calculated.

    Returns:
    --------
    dict
        A dictionary containing molecular descriptors. The keys represent descriptor names, and
        the values are the corresponding calculated values. The descriptors include:
        - 'formula': Molecular formula.
        - 'molecular_weight': Exact molecular weight.
        - 'n_heavy_atoms': Number of heavy atoms.
        - 'unique_bonds': Count of unique bonds.
        - 'bertz': Bertz complexity.
        - 'wiener_index': Wiener index.
        - 'balaban_index': Balaban index.
        - 'randic_index': Randic index.
        - 'kirchhoff_index': Kirchhoff index.
        - 'spacial_score': Spacial score.
        - 'chirality': Number of chiral centers.
        - 'fcfp4': FCFP_4 fingerprint.
        - 'bottcher': Bottcher complexity.
        - 'proudfoot': Proudfoot complexity.
        - 'mc1': Molecular connectivity index (MC1).
        - 'mc2': Molecular connectivity index (MC2).
        - Additional descriptors calculated using `get_mol_descriptors`.

    Notes:
    ------
    - The function combines descriptors from RDKit and custom calculations into a single dictionary.
    - If any descriptor calculation fails, it may not be included in the output dictionary.
    """
    out_dict = {'formula': rdMolDescriptors.CalcMolFormula(mol),
                "molecular_weight": rdMolDescriptors.CalcExactMolWt(mol),
                "n_heavy_atoms": rdMolDescriptors.CalcNumHeavyAtoms(mol),
                'unique_bonds': count_unique_bonds(mol),}
                # 'bertz': bertz(mol),
                # 'wiener_index': wiener_index(mol),
                # 'balaban_index': balaban_index(mol),
                # 'randic_index': randic_index(mol),
                # 'kirchhoff_index': kirchhoff_index(mol),
                # 'spacial_score': spacial_score(mol),
                # 'chirality': get_chirality(mol),
                # 'fcfp4': fcfp4(mol),
                # 'bottcher': bottcher(mol),
                # 'proudfoot': proudfoot(mol),
                # 'mc1': mc1(mol),
                # 'mc2': mc2(mol)}

    # batch_dict = get_mol_descriptors(mol)
    # Combine the descriptors into a single dictionary
    # out_dict.update(batch_dict)
    return out_dict
