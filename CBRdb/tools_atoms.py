import math
import os
import re
import shutil
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd
from ase import Atoms
from ase.calculators.orca import ORCA
from ase.calculators.orca import OrcaProfile
from ase.io import read
from ase.optimize import BFGS
from ase.thermochemistry import IdealGasThermo
from ase.units import Hartree
from ase.vibrations import Vibrations
from mace.calculators import mace_omol
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import PointGroupAnalyzer
from rdkit import Chem as Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import Mol

from .tools_mols import standardize_mol


def smi_to_atoms(smiles: str) -> tuple[Atoms, int, int]:
    """
    Convert a SMILES string to an ASE Atoms object, charge, and spin multiplicity.

    This function takes a SMILES string, sanitizes it, adds explicit hydrogens,
    standardizes the molecule, and converts it to an ASE Atoms object. It also
    calculates the charge and spin multiplicity of the molecule.

    Parameters:
    -----------
    smiles : str
        The SMILES string representing the molecule.

    Returns:
    --------
    tuple[Atoms, int, int]
        A tuple containing:
        - Atoms: The ASE Atoms object representing the molecule.
        - int: The formal charge of the molecule.
        - int: The spin multiplicity of the molecule.

    Raises:
    -------
    ValueError
        If the SMILES string cannot be parsed into a valid molecule.
    """
    mol = Chem.MolFromSmiles(smiles, sanitize=True)  # Parse the SMILES string into an RDKit molecule.
    mol = Chem.AddHs(mol)  # Add explicit hydrogens to the molecule.
    mol = standardize_mol(mol)  # Standardize the molecule structure.
    if mol is None:
        raise ValueError(f"Failed to parse SMILES string: {smiles}")  # Raise an error if the molecule is invalid.

    return mol_to_atoms(mol), get_charge(mol), get_spin_multiplicity(mol)  # Convert to Atoms and calculate properties.


def mol_to_atoms(mol: Mol, optimise: bool = True) -> Atoms:
    """
    Convert an RDKit molecule to an ASE Atoms object via an SDF file.

    Parameters:
    -----------
    mol : rdkit.Chem.rdchem.Mol
        RDKit molecule object to be converted.
    optimise : bool, optional
        Whether to optimise the molecule geometry before conversion (default is True).

    Returns:
    --------
    ase.Atoms
        Atoms object representing the molecule.
    """
    mol = standardize_mol(mol)
    mol = Chem.AddHs(mol)
    # If optimisation is enabled, embed and optimise the molecule using RDKit
    if optimise:
        AllChem.EmbedMolecule(mol, maxAttempts=5000, useRandomCoords=True, randomSeed=0xf00d)
        AllChem.MMFFOptimizeMolecule(mol)

    # Create a temporary file to store the molecule in SDF format
    with tempfile.NamedTemporaryFile(suffix='.sdf', delete=False) as temp_file:
        sdf_path = temp_file.name  # Get the path to the temporary file

    # Write the RDKit molecule to the SDF file
    writer = Chem.SDWriter(sdf_path)
    writer.write(mol)
    writer.close()

    # Read the SDF file and convert it to an ASE Atoms object
    atoms = read(sdf_path)

    # Remove the temporary SDF file to clean up
    os.remove(sdf_path)

    # Return the ASE Atoms object
    return atoms


def get_charge(mol: Mol) -> int:
    """
    Calculate the formal charge of a molecule.

    Parameters:
    -----------
    mol : rdkit.Chem.rdchem.Mol
        An RDKit molecule object

    Returns:
    --------
    int
        The formal charge of the molecule.
    """
    return Chem.GetFormalCharge(mol)


def _calc_unpaired(capacity: int, electrons: int) -> int:
    """
    Calculate the number of unpaired electrons in an atomic subshell.

    This function determines the number of unpaired electrons in a subshell
    based on its electron capacity and the number of electrons present.

    Parameters:
    -----------
    capacity : int
        The maximum number of electrons the subshell can hold.
    electrons : int
        The number of electrons present in the subshell.

    Returns:
    --------
    int
        The number of unpaired electrons in the subshell.
        If the number of electrons is less than or equal to half the capacity,
        all electrons are unpaired. Otherwise, the number of unpaired electrons
        is calculated based on the pairing rules.
    """
    orbitals = capacity // 2  # Calculate the number of orbitals in the subshell.
    return electrons if electrons <= orbitals else 2 * orbitals - electrons


def _aufbau_multiplicity(z: int) -> int:
    """
    Calculate the spin multiplicity of an atom based on the Aufbau principle.

    This function determines the spin multiplicity (2S+1) of an atom by filling
    its subshells with electrons according to the Aufbau principle. It calculates
    the number of unpaired electrons in the subshells and adds 1 to compute the
    spin multiplicity.

    Parameters:
    -----------
    z : int
        The atomic number of the atom.

    Returns:
    --------
    int
        The spin multiplicity of the atom (2S+1).

    Notes:
    ------
    - Subshells are filled in order of increasing energy levels.
    - The function uses the `_calc_unpaired` helper function to calculate the
      number of unpaired electrons in each subshell.
    """
    subshells = [
        ('1s', 2), ('2s', 2), ('2p', 6), ('3s', 2), ('3p', 6),
        ('4s', 2), ('3d', 10), ('4p', 6), ('5s', 2), ('4d', 10),
        ('5p', 6), ('6s', 2), ('4f', 14), ('5d', 10), ('6p', 6),
        ('7s', 2), ('5f', 14), ('6d', 10), ('7p', 6),
    ]
    remaining, unpaired = z, 0
    for _, cap in subshells:
        if remaining == 0:
            break
        n = min(cap, remaining)
        remaining -= n
        unpaired += _calc_unpaired(cap, n)
    return unpaired + 1  # 2S+1


def get_spin_multiplicity(mol: Chem.Mol) -> int:
    """
    Determine the spin multiplicity of a molecule.

    This function calculates the spin multiplicity (2S+1) of a molecule based on its structure.
    It first checks for explicit spin multiplicity properties, handles isolated atoms using
    predefined exceptions and the Aufbau principle, and finally calculates the multiplicity
    based on the number of unpaired electrons (radicals) in the molecule.

    Parameters:
    -----------
    mol : rdkit.Chem.Mol
        An RDKit molecule object.

    Returns:
    --------
    int
        The spin multiplicity of the molecule.

    Notes:
    ------
    - Spin multiplicity is calculated as (2S+1), where S is the total spin.
    - For isolated atoms, predefined exceptions and the Aufbau principle are used.
    - For molecules, the radical count determines the spin multiplicity.
    """
    mol = Chem.AddHs(mol)  # Add explicit hydrogens to the molecule.

    # 1 – Check for explicit spin multiplicity properties.
    for key in ("spinMultiplicity", "SpinMultiplicity"):
        if mol.HasProp(key):
            return int(mol.GetProp(key))  # Return the explicitly defined spin multiplicity.

    # 2 – Handle isolated atoms.
    if mol.GetNumAtoms() == 1:
        _EXCEPTIONS = {24: 7, 29: 2, 42: 7, 47: 2}  # Predefined exceptions for Cr, Cu, Mo, Ag.
        _GROUND_STATE_MULTIPLICITY = {
            z: _EXCEPTIONS.get(z, _aufbau_multiplicity(z))  # Use exceptions or the Aufbau principle.
            for z in range(1, 118 + 1)  # Atomic numbers from 1 to 118.
        }
        z = mol.GetAtomWithIdx(0).GetAtomicNum()  # Get the atomic number of the isolated atom.
        return _GROUND_STATE_MULTIPLICITY.get(z, 1)  # Return the ground-state multiplicity.

    # 3 – Calculate spin multiplicity for molecules based on radical count.
    n_rad = sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())  # Count unpaired electrons (radicals).
    return (n_rad + 1) if n_rad else 1  # Return (radical count + 1) or 1 if no radicals are present.


def standardise_smiles(smi):
    """
    Standardise a SMILES (Simplified Molecular Input Line Entry System) string.

    This function takes a SMILES string, adds explicit hydrogens to the molecule,
    and returns a standardised version of the SMILES string with canonical, isomeric,
    and Kekulé representations.

    Parameters:
    -----------
    smi : str
        The input SMILES string to be standardised.

    Returns:
    --------
    str
        The standardised SMILES string.

    Raises:
    -------
    ValueError
        If the input SMILES string is invalid and cannot be parsed.
    """
    # Convert the SMILES string to an RDKit molecule object and add explicit hydrogens
    mol = Chem.AddHs(Chem.MolFromSmiles(smi, sanitize=True))

    # Raise an error if the molecule could not be created
    if not mol:
        raise ValueError(f"Invalid SMILES: {smi}")

    # Return the standardised SMILES string with specified options
    return Chem.MolToSmiles(mol, isomericSmiles=True, kekuleSmiles=True, canonical=True)


def orca_calc_preset(orca_path=None,
                     directory=None,
                     calc_type='DFT',
                     xc='wB97X',
                     charge=0,
                     multiplicity=1,
                     basis_set='def2-SVP',
                     n_procs=10,
                     f_solv=False,
                     f_disp=False,
                     atom_list=None,
                     calc_extra=None,
                     blocks_extra=None,
                     scf_option=None):
    """
    Create and configure an ORCA calculator preset for quantum chemistry calculations.

    Parameters:
    -----------
    orca_path : str, optional
        Path to the ORCA executable. If None, it will attempt to read from the environment variable 'ORCA_PATH'.
    directory : str, optional
        Directory where the calculation will be performed. Defaults to a temporary directory.
    calc_type : str, optional
        Type of calculation to perform (e.g., 'DFT', 'MP2', 'CCSD', 'QM/XTB2'). Default is 'DFT'.
    xc : str, optional
        Exchange-correlation functional to use. Default is 'wB97X'.
    charge : int, optional
        Total charge of the system. Default is 0.
    multiplicity : int, optional
        Spin multiplicity of the system. Default is 1.
    basis_set : str, optional
        Basis set to use for the calculation. Default is 'def2-SVP'.
    n_procs : int, optional
        Number of processors to use. Default is 10.
    f_solv : bool or str, optional
        Solvent model to use. If True, defaults to 'WATER'. Default is False (no solvent).
    f_disp : bool or str, optional
        Dispersion correction to use. If True, defaults to 'D4'. Default is False (no dispersion correction).
    atom_list : list, optional
        List of atoms for QM/MM calculations. Only used if `calc_type` is 'QM/XTB2'. Default is None.
    calc_extra : str, optional
        Additional calculation options to include in the ORCA input. Default is None.
    blocks_extra : str, optional
        Additional ORCA input blocks to include. Default is None.
    scf_option : str, optional
        Additional SCF options to include in the ORCA input. Default is None.

    Returns:
    --------
    ORCA
        Configured ORCA calculator object.
    """
    if orca_path is None:
        # Try and read the path from the environment
        orca_path = os.environ.get('ORCA_PATH')
    if directory is None:
        # Create a temporary directory for the calculation
        directory = os.path.join(tempfile.mkdtemp(), 'orca')

    # Create an ORCA profile with the specified command
    profile = OrcaProfile(command=orca_path)

    # Configure the number of processors
    if n_procs > 1:
        inpt_procs = '%pal nprocs {} end'.format(n_procs)
    else:
        inpt_procs = ''

    # Configure the solvent model
    if f_solv is not None and f_solv is not False:
        if f_solv:
            f_solv = 'WATER'
        inpt_solv = '''
                                              %CPCM SMD TRUE
                                                  SMDSOLVENT "{}"
                                              END'''.format(f_solv)
    else:
        inpt_solv = ''

    # Configure the dispersion correction
    if f_disp is None or f_disp is False:
        inpt_disp = ''
    else:
        if f_disp:
            f_disp = 'D4'
        inpt_disp = f_disp

    # Configure QM/MM atom list for QM/XTB2 calculations
    if atom_list is not None and calc_type == 'QM/XTB2':
        inpt_xtb = '''
                                              %QMMM QMATOMS {{}} END END
                                              '''.format(str(atom_list).strip('[').strip(']'))
    else:
        inpt_xtb = ''

    # Add any additional input blocks
    if blocks_extra is None:
        blocks_extra = ''

    # Combine all input blocks
    inpt_blocks = inpt_procs + inpt_solv + blocks_extra

    # Configure the main calculation input based on the calculation type
    if calc_type == 'DFT':
        inpt_simple = '{} {} {}'.format(xc, inpt_disp, basis_set)
    elif calc_type == 'MP2':
        inpt_simple = 'DLPNO-{} {} {}/C'.format(calc_type, basis_set, basis_set)
    elif calc_type == 'CCSD':
        inpt_simple = 'DLPNO-{}(T) {} {}/C'.format(calc_type, basis_set, basis_set)
    elif calc_type == 'QM/XTB2':
        inpt_simple = '{} {} {} {}'.format(calc_type, xc, inpt_disp, basis_set)
        inpt_blocks = inpt_procs + inpt_solv + inpt_xtb
    else:
        inpt_simple = '{} {}'.format(calc_type, basis_set)

    if multiplicity > 1:
        if calc_type == 'DFT' or calc_type == 'QM/XTB2':
            inpt_simple = 'UKS  ' + inpt_simple
        elif calc_type == 'MP2' or calc_type == 'CCSD':
            inpt_simple = 'UKS ' + inpt_simple

    # Add the SCF option if provided
    if scf_option is not None:
        inpt_simple += ' ' + scf_option

    # Add any extra calculation options
    if calc_extra is not None:
        inpt_simple += ' ' + calc_extra

    # Create and return the ORCA calculator object
    calc = ORCA(
        profile=profile,
        charge=charge,
        mult=multiplicity,
        directory=directory,
        orcasimpleinput=inpt_simple,
        orcablocks=inpt_blocks
    )
    return calc


def optimise_atoms(atoms,
                   charge=0,
                   multiplicity=1,
                   orca_path=None,
                   xc='r2SCAN-3c',
                   basis_set='def2-QZVP',
                   tight_opt=False,
                   tight_scf=False,
                   f_solv=False,
                   f_disp=False,
                   n_procs=10):
    """
    Optimise the geometry of a molecule using the ORCA quantum chemistry package.

    This function sets up an ORCA calculation to optimise the geometry of a molecule
    represented by an ASE `Atoms` object. It supports various calculation options,
    including tight optimisation, solvent effects, and dispersion corrections.

    Parameters:
    -----------
    atoms : ase.Atoms
        An ASE `Atoms` object representing the molecule to be optimised.
    charge : int, optional
        Total charge of the molecule. Default is 0.
    multiplicity : int, optional
        Spin multiplicity of the molecule. Default is 1.
    orca_path : str, optional
        Path to the ORCA executable. If None, it will attempt to read from the environment variable 'ORCA_PATH'.
    xc : str, optional
        Exchange-correlation functional to use. Default is 'r2SCAN-3c'.
    basis_set : str, optional
        Basis set to use for the calculation. Default is 'def2-QZVP'.
    tight_opt : bool, optional
        Whether to use tight geometry optimisation. Default is False.
    tight_scf : bool, optional
        Whether to use tight SCF convergence criteria. Default is False.
    f_solv : bool, optional
        Whether to include solvent effects in the calculation. Default is False.
    f_disp : bool, optional
        Whether to include dispersion corrections in the calculation. Default is False.
    n_procs : int, optional
        Number of processors to use for the calculation. Default is 10.

    Returns:
    --------
    ase.Atoms
        An ASE `Atoms` object representing the optimised geometry of the molecule.

    Raises:
    -------
    ValueError
        If the ORCA path cannot be determined or the calculation fails.
    """
    # Determine the ORCA path
    if orca_path is None:
        # Try to read the path from the environment variable
        orca_path = os.environ.get('ORCA_PATH')
    else:
        # Convert the provided path to an absolute path
        orca_path = os.path.abspath(orca_path)

    if tight_opt:
        # Set up geometry optimization and frequency calculation parameters
        opt_option = 'TIGHTOPT'
    else:
        # Set up frequency calculation parameters only
        opt_option = 'OPT'

    if tight_scf:
        # Set up tight SCF convergence parameters
        calc_extra = f'{opt_option} TIGHTSCF'
    else:
        # Use default SCF convergence parameters
        calc_extra = f'{opt_option}'

    # Create a temporary working directory
    with tempfile.TemporaryDirectory() as temp_dir:
        orca_file = os.path.join(temp_dir, "orca.xyz")

        # Set up the ORCA calculator with the specified parameters
        calc = orca_calc_preset(orca_path=orca_path,
                                directory=temp_dir,
                                charge=charge,
                                multiplicity=multiplicity,
                                xc=xc,
                                basis_set=basis_set,
                                n_procs=n_procs,
                                f_solv=f_solv,
                                f_disp=f_disp,
                                calc_extra=calc_extra)
        # Assign the calculator to the molecule
        atoms.calc = calc

        # Trigger the calculation to optimise the geometry
        _ = atoms.get_potential_energy()

        # Load the optimised geometry from the ORCA output file
        return read(orca_file, format="xyz")


def load_ir_data(filename):
    """
    Load IR spectrum data from ORCA output file into a pandas DataFrame.

    This function reads an ORCA output file containing IR frequency data,
    extracts the relevant information, and returns it as a pandas DataFrame.

    Parameters:
    -----------
    filename : str
        Path to the IR frequency output file.

    Returns:
    --------
    pd.DataFrame
        A DataFrame containing the following columns:
        - 'Mode': Mode number (int).
        - 'Frequency (cm^-1)': Frequency in inverse centimeters (float).
        - 'Epsilon': Epsilon value (float).
        - 'Intensity (km/mol)': Intensity in km/mol (float).

    Raises:
    -------
    ValueError
        If the IR spectrum data cannot be found in the file.
    """
    # Read the file
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Find the start of the IR spectrum data
    start_idx = None
    for i, line in enumerate(lines):
        if 'Mode   freq       eps      Int' in line:
            start_idx = i + 2  # Skip header and separator line
            break

    if start_idx is None:
        raise ValueError("Could not find IR spectrum data in the file")

    # Extract data
    data = []
    for i in range(start_idx, len(lines)):
        line = lines[i].strip()
        if not line or line.startswith('*') or line.startswith('The first'):
            break

        # Parse the line using regex to handle varying whitespace
        match = re.match(r'\s*(\d+):\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)', line)
        if match:
            mode = int(match.group(1))  # Mode number
            freq = float(match.group(2))  # Frequency in cm^-1
            eps = float(match.group(3))  # Epsilon value
            intensity = float(match.group(4))  # Intensity in km/mol
            data.append([mode, freq, eps, intensity])

    # Create DataFrame
    df = pd.DataFrame(data, columns=['Mode', 'Frequency (cm^-1)', 'Epsilon', 'Intensity (km/mol)'])

    return df


def load_raman_data(filename):
    """
    Load Raman spectrum data from ORCA output file into a pandas DataFrame.

    This function reads an ORCA output file containing Raman frequency data,
    extracts the relevant information, and returns it as a pandas DataFrame.

    Parameters:
    -----------
    filename : str
        Path to the Raman frequency output file.

    Returns:
    --------
    pd.DataFrame
        A DataFrame containing the following columns:
        - 'Mode': Mode number (int).
        - 'Frequency (cm^-1)': Frequency in inverse centimeters (float).
        - 'Activity': Activity (float).
        - 'Depolarization': Depolarization value (float).

    Raises:
    -------
    ValueError
        If the Raman spectrum data cannot be found in the file.
    """
    # Read the file
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Find the start of the Raman spectrum data
    start_idx = None
    for i, line in enumerate(lines):
        if 'Mode    freq (cm**-1)   Activity   Depolarization' in line:
            start_idx = i + 2  # Skip header and separator line
            break

    if start_idx is None:
        raise ValueError("Could not find Raman spectrum data in the file")

    # Extract data
    data = []
    for i in range(start_idx, len(lines)):
        line = lines[i].strip()
        if not line or line.startswith('The first'):
            break

        # Parse the line using regex to handle varying whitespace
        match = re.match(r'\s*(\d+):\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)', line)
        if match:
            mode = int(match.group(1))  # Mode number
            freq = float(match.group(2))  # Frequency in cm^-1
            activity = float(match.group(3))  # Activity in A^4 amu^-1
            depolarization = float(match.group(4))  # Depolarization value
            data.append([mode, freq, activity, depolarization])

    # Create DataFrame
    df = pd.DataFrame(data, columns=['Mode', 'Frequency (cm^-1)', 'Activity (A^4 amu^-1)', 'Depolarization'])

    return df


def load_vib_data(filename):
    """
    Load vibrational spectrum data from a file into a pandas DataFrame.

    Parameters:
    -----------
    filename : str
        Path to the file containing vibrational spectrum data.

    Returns:
    --------
    pd.DataFrame
        A DataFrame containing the following columns:
        - 'Mode': Mode number (int).
        - 'Frequency (cm^-1)': Frequency in inverse centimeters (float).

    Raises:
    -------
    ValueError
        If the vibrational spectrum data cannot be found in the file.
    """
    # Read the file
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Find the start of the IR spectrum data
    start_idx = None
    for i, line in enumerate(lines):
        if 'Mode   freq       eps      Int' in line:
            start_idx = i + 2  # Skip header and separator line
            break

    if start_idx is None:
        raise ValueError("Could not find IR spectrum data in the file")

    # Extract data
    data = []
    for i in range(start_idx, len(lines)):
        line = lines[i].strip()
        if not line or line.startswith('*') or line.startswith('The first'):
            break

        # Parse the line using regex to handle varying whitespace
        match = re.match(r'\s*(\d+):\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)', line)
        if match:
            mode = int(match.group(1))  # Mode number
            freq = float(match.group(2))  # Frequency in cm^-1
            data.append([mode, freq])

    # Create DataFrame
    df = pd.DataFrame(data, columns=['Mode', 'Frequency (cm^-1)'])

    return df


def calculate_vib_spectrum(atoms,
                           charge=0,
                           multiplicity=1,
                           orca_path=None,
                           xc='r2SCAN-3c',
                           basis_set='def2-QZVP',
                           tight_opt=False,
                           tight_scf=False,
                           f_solv=False,
                           f_disp=False,
                           n_procs=10):
    """
    Calculate vibrational spectrum data using the ORCA quantum chemistry package.

    This function sets up and performs a vibrational spectrum calculation for a molecule
    represented by an ASE `Atoms` object. It computes IR, Raman, and vibrational spectrum data.

    Parameters:
    -----------
    atoms : ase.Atoms
        An ASE `Atoms` object representing the molecule.
    charge : int, optional
        Total charge of the molecule. Default is 0.
    multiplicity : int, optional
        Spin multiplicity of the molecule. Default is 1.
    orca_path : str, optional
        Path to the ORCA executable. If None, it will attempt to read from the environment variable 'ORCA_PATH'.
    xc : str, optional
        Exchange-correlation functional to use. Default is 'r2SCAN-3c'.
    basis_set : str, optional
        Basis set to use for the calculation. Default is 'def2-QZVP'.
    tight_opt : bool, optional
        Whether to use tight geometry optimisation. Default is False.
    tight_scf : bool, optional
        Whether to use tight SCF convergence criteria. Default is False.
    f_solv : bool, optional
        Whether to include solvent effects in the calculation. Default is False.
    f_disp : bool, optional
        Whether to include dispersion corrections in the calculation. Default is False.
    n_procs : int, optional
        Number of processors to use for the calculation. Default is 10.

    Returns:
    --------
    tuple
        A tuple containing three pandas DataFrames:
        - data_ir : pd.DataFrame
            IR spectrum data.
        - data_raman : pd.DataFrame
            Raman spectrum data.
        - data_vib : pd.DataFrame
            Vibrational spectrum data.

    Raises:
    -------
    ValueError
        If the ORCA path cannot be determined or the calculation fails.
    """
    # Determine the ORCA path
    if orca_path is None:
        # Try to read the path from the environment variable
        orca_path = os.environ.get('ORCA_PATH')
    else:
        # Convert the provided path to an absolute path
        orca_path = os.path.abspath(orca_path)

    # Get the total number of electrons in the system
    total_electrons = get_total_electrons(atoms)
    # Prevent too many processors being used
    if n_procs > total_electrons:
        n_procs = round_to_nearest_two(total_electrons - 2)

    # Set optimization flags
    opt_flag = 'TIGHTOPT' if tight_opt else 'OPT'
    if len(atoms) == 1:  # Skip optimization for single atoms
        opt_flag = ''

    # Set SCF flags
    scf_flag = 'TIGHTSCF' if tight_scf else ''
    calc_extra = f'{opt_flag} {scf_flag} FREQ'.strip()

    blocks_extra = '''
                          %ELPROP
                              POLAR 1
                          END'''

    # Create a temporary working directory
    with tempfile.TemporaryDirectory() as temp_dir:
        orca_file = os.path.join(temp_dir, 'orca.out')

        # Set up the ORCA calculator with the specified parameters
        calc = orca_calc_preset(orca_path=orca_path,
                                directory=temp_dir,
                                charge=charge,
                                multiplicity=multiplicity,
                                xc=xc,
                                basis_set=basis_set,
                                n_procs=n_procs,
                                f_solv=f_solv,
                                f_disp=f_disp,
                                calc_extra=calc_extra,
                                blocks_extra=blocks_extra)

        # Attach the calculator to the ASE Atoms object
        atoms.calc = calc
        try:
            # Perform the calculation (this will write the output to the ORCA file)
            _ = atoms.get_potential_energy()

            # Load IR spectrum data
            data_ir = load_ir_data(orca_file)

            # Load Raman spectrum data
            data_raman = load_raman_data(orca_file)

            # Load vibrational spectrum data
            data_vib = load_vib_data(orca_file)
        except Exception as e:
            data_ir = None
            data_raman = None
            data_vib = None
            print(f"Failed to perform vibrational spectrum calculation: {e}")

        return data_ir, data_raman, data_vib


def get_total_electrons(atoms: Atoms) -> int:
    """
    Calculate the total number of electrons in a molecule.

    This function computes the total number of electrons in a molecule
    represented by an ASE `Atoms` object. It sums the atomic numbers (Z)
    of all atoms in the molecule and adjusts for the explicit charge
    provided in the `Atoms.info` dictionary.

    Parameters:
    -----------
    atoms : ase.Atoms
        An ASE `Atoms` object representing the molecule.

    Returns:
    --------
    int
        The total number of electrons in the molecule, corrected for its charge.
    """
    # Sum atomic numbers (Z) for every atom in the molecule
    n_electrons = int(np.sum(atoms.get_atomic_numbers()))

    # Correct for explicit total charge, if provided in the `Atoms.info` dictionary
    charge = atoms.info.get('charge', 0.0)
    n_electrons -= int(round(charge))

    return n_electrons


def round_to_nearest_two(number):
    """
    Round a number to the nearest multiple of 2.
    If the result would be 0, return 1 instead.

    Parameters:
    -----------
    number : float or int
        The number to be rounded

    Returns:
    --------
    int
        The nearest multiple of 2, or 1 if result would be 0
    """
    # Round to nearest multiple of 2
    result = round(number / 2) * 2

    # If result is 0, set it to 1
    if result == 0:
        result = 1

    return result


def calculate_ccsd_energy(atoms,
                          charge=0,
                          multiplicity=1,
                          orca_path=None,
                          basis_set='def2-TZVPP',
                          n_procs=10):
    """
    Perform a CCSD (Coupled Cluster Single and Double) energy calculation using the ORCA quantum chemistry package.

    This function sets up and executes a CCSD energy calculation for a molecule represented by an ASE `Atoms` object.
    It ensures that the number of processors used does not exceed the total number of electrons in the system.

    Parameters:
    -----------
    atoms : ase.Atoms
        An ASE `Atoms` object representing the molecule.
    charge : int, optional
        Total charge of the molecule. Default is 0.
    multiplicity : int, optional
        Spin multiplicity of the molecule. Default is 1.
    orca_path : str, optional
        Path to the ORCA executable. If None, it will attempt to read from the environment variable 'ORCA_PATH'.
    basis_set : str, optional
        Basis set to use for the calculation. Default is 'def2-TZVPP'.
    n_procs : int, optional
        Number of processors to use for the calculation. Default is 10.

    Returns:
    --------
    float
        The CCSD energy of the molecule in eV.

    Raises:
    -------
    ValueError
        If the number of processors exceeds the adjusted limit based on the total number of electrons.
    """
    # If no ORCA path is provided, try to read it from the environment variable
    orca_path = os.path.abspath(orca_path or os.getenv('ORCA_PATH', 'orca'))

    # Get the total number of electrons in the system
    total_electrons = get_total_electrons(atoms)
    # Prevent too many processors being used
    if n_procs > total_electrons:
        n_procs = round_to_nearest_two(total_electrons - 2)

    # Create a temporary directory for the ORCA calculation
    with tempfile.TemporaryDirectory() as temp_dir:
        # Set up the ORCA calculator with the specified parameters
        calc = orca_calc_preset(orca_path=orca_path,
                                directory=temp_dir,
                                calc_type='CCSD',
                                charge=charge,
                                multiplicity=multiplicity,
                                basis_set=basis_set,
                                n_procs=n_procs)
        # Attach the ORCA calculator to the ASE Atoms object
        atoms.calc = calc

        # Perform the energy calculation
        return atoms.get_potential_energy()


def grab_value(orca_file, term, splitter):
    """
    Extract a specific numerical value from an ORCA output file.

    This function reads an ORCA output file in reverse order, searches for a specific term,
    and extracts the numerical value associated with it. The value is converted from Hartree
    units to eV using the ASE `Hartree` constant.

    Parameters:
    -----------
    orca_file : str
        Path to the ORCA output file.
    term : str
        The term to search for in the file.
    splitter : str
        The delimiter used to split the line containing the term.

    Returns:
    --------
    float or None
        The extracted value in eV, or None if the term is not found.
    """
    with open(orca_file, 'r') as f:
        for line in reversed(f.readlines()):
            if term in line:
                return float(line.split(splitter)[-1].split('Eh')[0]) * Hartree
        return None


def calculate_free_energy(atoms,
                          charge=0,
                          multiplicity=1,
                          temperature=None,
                          pressure=None,
                          orca_path=None,
                          xc='r2SCAN-3c',
                          basis_set='def2-QZVP',
                          tight_opt=False,
                          tight_scf=False,
                          f_solv=False,
                          f_disp=False,
                          n_procs=10,
                          use_ccsd=False,
                          ccsd_energy=None):
    """
    Calculate the Gibbs free energy, enthalpy, and entropy of a molecule.

    This function performs a quantum chemistry calculation using the ORCA package to compute
    the Gibbs free energy, enthalpy, and entropy of a molecule represented by an ASE `Atoms` object.
    It supports temperature and pressure adjustments, CCSD energy calculations, and various ORCA options.

    Parameters:
    -----------
    atoms : ase.Atoms
        An ASE `Atoms` object representing the molecule.
    charge : int, optional
        Total charge of the molecule. Default is 0.
    multiplicity : int, optional
        Spin multiplicity of the molecule. Default is 1.
    temperature : float, optional
        Temperature in Kelvin for the calculation. Default is None.
    pressure : float, optional
        Pressure in atm for the calculation. Default is None.
    orca_path : str, optional
        Path to the ORCA executable. If None, it will attempt to read from the environment variable 'ORCA_PATH'.
    xc : str, optional
        Exchange-correlation functional to use. Default is 'r2SCAN-3c'.
    basis_set : str, optional
        Basis set to use for the calculation. Default is 'def2-QZVP'.
    tight_opt : bool, optional
        Whether to use tight geometry optimization. Default is False.
    tight_scf : bool, optional
        Whether to use tight SCF convergence criteria. Default is False.
    f_solv : bool, optional
        Whether to include solvent effects in the calculation. Default is False.
    f_disp : bool, optional
        Whether to include dispersion corrections in the calculation. Default is False.
    n_procs : int, optional
        Number of processors to use for the calculation. Default is 10.
    use_ccsd : bool, optional
        Whether to use CCSD energy calculations. Default is False.
    ccsd_energy : float, optional
        Precomputed CCSD energy in eV. If None, CCSD energy will be calculated if `use_ccsd` is True.

    Returns:
    --------
    tuple
        A tuple containing:
        - energy : float
            The Gibbs free energy in eV.
        - enthalpy : float
            The enthalpy in eV.
        - entropy : float
            The entropy correction in eV.

    Raises:
    -------
    ValueError
        If the CCSD energy calculation fails or the ORCA setup is invalid.
    """
    # Determine the ORCA path
    orca_path = os.path.abspath(orca_path or os.getenv('ORCA_PATH', 'orca'))

    # Set optimization flags
    opt_flag = 'TIGHTOPT' if tight_opt else 'OPT'
    if len(atoms) == 1:  # Skip optimization for single atoms
        opt_flag = ''

    # Set SCF flags
    scf_flag = 'TIGHTSCF' if tight_scf else ''
    calc_extra = f'{opt_flag} {scf_flag} FREQ'.strip()

    # Set up the %thermo block for this temperature and pressure
    if temperature is not None and pressure is None:
        blocks_extra = f'''
                                  %freq
                                      Temp {temperature}
                                  end
                                  '''
    elif pressure is not None and temperature is None:
        blocks_extra = f'''
                                          %freq
                                              Pressure {pressure}
                                          end
                                          '''
    elif pressure is None and temperature is not None:
        blocks_extra = f'''
                                          %freq
                                              Temp {temperature}
                                              Pressure {pressure}
                                          end
                                          '''
    else:
        blocks_extra = None

    # Perform CCSD energy calculation if required and not provided
    if use_ccsd and ccsd_energy is None:
        ccsd_energy = calculate_ccsd_energy(atoms,
                                            orca_path=orca_path,
                                            charge=charge,
                                            multiplicity=multiplicity,
                                            n_procs=n_procs)
        if ccsd_energy is None:
            raise ValueError("CCSD energy calculation failed. Please check the ORCA setup.")

    # Create a temporary directory for the calculation
    with tempfile.TemporaryDirectory() as temp_dir:
        orca_file = os.path.join(temp_dir, 'orca.out')

        # Set up the ORCA calculator
        calc = orca_calc_preset(orca_path=orca_path,
                                directory=temp_dir,
                                charge=charge,
                                multiplicity=multiplicity,
                                xc=xc,
                                basis_set=basis_set,
                                n_procs=n_procs,
                                f_solv=f_solv,
                                f_disp=f_disp,
                                calc_extra=calc_extra,
                                blocks_extra=blocks_extra)
        atoms.calc = calc

        # Trigger the calculation
        _ = atoms.get_potential_energy()

        # Extract entropy correction
        entropy = grab_value(orca_file, 'Total entropy correction', '...')

        # Calculate Gibbs free energy based on CCSD or DFT results
        if use_ccsd:
            g_e_ele = grab_value(orca_file, 'G-E(el)', '...')
            g_e_solv = grab_value(orca_file, 'Free-energy (cav+disp)', ':') if f_solv else 0.0
            energy = ccsd_energy + g_e_ele + g_e_solv
        else:
            energy = grab_value(orca_file, 'Final Gibbs free energy', '...')

        # Return energy, enthalpy, and entropy
        return energy, energy - entropy, entropy


def list_to_str(lst):
    """
    Convert a list of items to a comma-separated string.

    This function takes a list of items, converts each item to a string,
    and joins them with commas to create a single string.

    Parameters:
    -----------
    lst : list
        A list of items to be converted to a string.

    Returns:
    --------
    str
        A comma-separated string representation of the list.
    """
    lst = [str(item) for item in lst]  # Convert each item in the list to a string.
    return ', '.join(lst)  # Join the string representations with commas.


def calculate_hessian(atoms,
                      charge=0,
                      multiplicity=1,
                      orca_path=None,
                      xc='r2SCAN-3c',
                      basis_set='def2-QZVP',
                      tight_opt=False,
                      tight_scf=False,
                      f_solv=False,
                      f_disp=False,
                      n_procs=10):
    """
    Perform a Hessian matrix calculation using the ORCA quantum chemistry package.

    This function sets up and executes a Hessian matrix calculation for a molecule
    represented by an ASE `Atoms` object. It optimizes the geometry and computes
    the Hessian matrix, which is used for vibrational analysis.

    Parameters:
    -----------
    atoms : ase.Atoms
        An ASE `Atoms` object representing the molecule.
    charge : int, optional
        Total charge of the molecule. Default is 0.
    multiplicity : int, optional
        Spin multiplicity of the molecule. Default is 1.
    orca_path : str, optional
        Path to the ORCA executable. If None, it will attempt to read from the environment variable 'ORCA_PATH'.
    xc : str, optional
        Exchange-correlation functional to use. Default is 'r2SCAN-3c'.
    basis_set : str, optional
        Basis set to use for the calculation. Default is 'def2-QZVP'.
    tight_opt : bool, optional
        Whether to use tight geometry optimization. Default is False.
    tight_scf : bool, optional
        Whether to use tight SCF convergence criteria. Default is False.
    f_solv : bool, optional
        Whether to include solvent effects in the calculation. Default is False.
    f_disp : bool, optional
        Whether to include dispersion corrections in the calculation. Default is False.
    n_procs : int, optional
        Number of processors to use for the calculation. Default is 10.

    Returns:
    --------
    tuple
        A tuple containing:
        - atoms : ase.Atoms
            The optimized geometry of the molecule.
        - hessian_file : str
            Path to the file containing the Hessian matrix.

    Raises:
    -------
    ValueError
        If the ORCA path cannot be determined or the calculation fails.
    """
    # Determine the ORCA path
    if orca_path is None:
        # Try to read the path from the environment variable
        orca_path = os.environ.get('ORCA_PATH')
    else:
        # Convert the provided path to an absolute path
        orca_path = os.path.abspath(orca_path)

    if tight_opt:
        # Set up geometry optimization and frequency calculation parameters
        opt_option = 'TIGHTOPT'
    else:
        # Set up frequency calculation parameters only
        opt_option = 'OPT'

    if tight_scf:
        # Set up tight SCF convergence parameters
        calc_extra = f'{opt_option} TIGHTSCF FREQ'
    else:
        # Use default SCF convergence parameters
        calc_extra = f'{opt_option} FREQ'

    # Create a temporary directory for the ORCA calculation
    with tempfile.TemporaryDirectory() as temp_dir:

        # Set up the ORCA calculator with the specified parameters
        calc = orca_calc_preset(orca_path=orca_path,
                                directory=temp_dir,
                                charge=charge,
                                multiplicity=multiplicity,
                                xc=xc,
                                basis_set=basis_set,
                                n_procs=n_procs,
                                f_solv=f_solv,
                                f_disp=f_disp,
                                calc_extra=calc_extra)

        # Attach the ORCA calculator to the ASE Atoms object
        atoms.calc = calc

        # Perform the energy calculation
        _ = atoms.get_potential_energy()

        # Load the optimized geometry from the ORCA output file
        atoms_file = os.path.join(temp_dir, "orca.xyz")
        hessian_file = os.path.join(temp_dir, "orca.hess")
        return read(atoms_file, format="xyz"), hessian_file


def _calculate_free_energy_batch(atoms,
                                 hessian,
                                 temp,
                                 pressure,
                                 charge=0,
                                 multiplicity=1,
                                 orca_path=None,
                                 xc='r2SCAN-3c',
                                 basis_set='def2-QZVP',
                                 f_solv=False,
                                 f_disp=False,
                                 n_procs=10,
                                 ccsd_energy=None):
    orca_path = os.path.abspath(orca_path or os.getenv('ORCA_PATH', 'orca'))

    # Set up the %thermo block for this temperature and pressure
    blocks_extra = f'''
    %GEOM
         INHESSNAME "{hessian}"
    END
    
    %freq
        Temp {temp}
        Pressure {pressure}
    end
    '''

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = 'tmp'
        os.makedirs(temp_dir, exist_ok=True)
        orca_file = os.path.join(temp_dir, "orca.out")

        calc = orca_calc_preset(orca_path=orca_path,
                                directory=temp_dir,
                                charge=charge,
                                multiplicity=multiplicity,
                                xc=xc,
                                basis_set=basis_set,
                                n_procs=n_procs,
                                f_solv=f_solv,
                                f_disp=f_disp,
                                calc_extra='PRINTTHERMOCHEM',
                                blocks_extra=blocks_extra)

        atoms.calc = calc
        # Trigger the calculation to optimise the system.
        try:
            _ = atoms.get_potential_energy()
        except:
            pass

        entropy = grab_value(orca_file, 'Total entropy correction', '...')

        if ccsd_energy is not None:
            g_e_ele = grab_value(orca_file, 'G-E(el)', '...')
            g_e_solv = grab_value(orca_file, 'Free-energy (cav+disp)', ':') if f_solv else 0
            energy = ccsd_energy + g_e_ele + g_e_solv
        else:
            energy = grab_value(orca_file, 'Final Gibbs free energy', '...')

        return energy, energy - entropy, entropy


def calculate_free_energy_batch(atoms,
                                hessian,
                                temp,
                                pressure,
                                charge=0,
                                multiplicity=1,
                                orca_path=None,
                                xc='r2SCAN-3c',
                                basis_set='def2-QZVP',
                                tight_opt=False,
                                tight_scf=False,
                                f_solv=False,
                                f_disp=False,
                                n_procs=10,
                                ccsd_energy=None):
    orca_path = os.path.abspath(orca_path or os.getenv('ORCA_PATH', 'orca'))

    # Set up the %thermo block for this temperature and pressure
    blocks_extra = f'''
    %GEOM
         INHESSNAME "{hessian}"
    END

    %freq
        Temp {temp}
        Pressure {pressure}
    end
    '''

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = 'tmp'
        os.makedirs(temp_dir, exist_ok=True)
        orca_file = os.path.join(temp_dir, "orca.out")

        calc = orca_calc_preset(orca_path=orca_path,
                                directory=temp_dir,
                                charge=charge,
                                multiplicity=multiplicity,
                                xc=xc,
                                basis_set=basis_set,
                                n_procs=n_procs,
                                f_solv=f_solv,
                                f_disp=f_disp,
                                calc_extra='PRINTTHERMOCHEM',
                                blocks_extra=blocks_extra)

        atoms.calc = calc
        # Trigger the calculation to optimise the system.
        try:
            _ = atoms.get_potential_energy()
        except:
            pass

        entropy = grab_value(orca_file, 'Total entropy correction', '...')

        if ccsd_energy is not None:
            g_e_ele = grab_value(orca_file, 'G-E(el)', '...')
            g_e_solv = grab_value(orca_file, 'Free-energy (cav+disp)', ':') if f_solv else 0
            energy = ccsd_energy + g_e_ele + g_e_solv
        else:
            energy = grab_value(orca_file, 'Final Gibbs free energy', '...')

        return energy, energy - entropy, entropy


def get_formation_references(mol):
    """
    Generate a list of reference molecules for calculating formation energies.

    This function determines the elemental composition of a molecule, checks for unsupported elements,
    and generates a list of standard reference molecules based on the elemental composition.

    Parameters:
    -----------
    mol : rdkit.Chem.Mol
        An RDKit molecule object representing the molecule.

    Returns:
    --------
    list[tuple[str, float]]
        A list of tuples where each tuple contains:
        - str: The SMILES string of the reference molecule.
        - float: The count of the reference molecule needed for the calculation.

    Raises:
    -------
    ValueError
        If the molecule contains unsupported elements.

    Notes:
    ------
    - Supported elements are defined in the `supported_set`.
    - Reference molecules are predefined for common elements like H, Cl, N, O, I, F, and Br.
      For other elements, the atomic symbol is used as the reference.
    """
    # Get elemental composition
    atom_counts = defaultdict(int)
    symbols = set()
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        atom_counts[symbol] += 1.0
        symbols.add(symbol)

    # Define the set of supported elements
    supported_set = {'H', 'Ba', 'Cu', 'S', 'K', 'Cl', 'Rb', 'B', 'N', 'Se', 'Te', 'O', 'Fe', 'Co', 'Mg', 'Ge', 'I',
                     'Tl', 'Pt', 'Xe', 'Zn', 'Gd', 'Cd', 'C', 'Al', 'F', 'Li', 'Ca', 'Ni', 'Th', 'Sr', 'Sn', 'Au', 'Ag',
                     'V', 'Pu', 'Sb', 'Mn', 'Cr', 'Mo', 'Ra', 'Pb', 'Bi', 'As', 'Hg', 'Si', 'Br', 'Rn', 'P', 'W', 'Be',
                     'Na'}
    # Raise an error if unsupported elements are found
    if not symbols.issubset(supported_set):
        raise ValueError(f"Unsupported elements in the molecule: {symbols - supported_set}")

    # Standard reference molecules
    references = []
    # Loop over the atom_counts dictionary
    for symbol, count in atom_counts.items():
        if count > 0:
            if symbol == 'H':
                references.append(("[H][H]", atom_counts['H'] / 2.0))
            elif symbol == 'Cl':
                references.append(("Cl-Cl", atom_counts['Cl'] / 2.0))
            elif symbol == 'N':
                references.append(("N#N", atom_counts['N'] / 2.0))
            elif symbol == 'O':
                references.append(("O=O", atom_counts['O'] / 2.0))
            elif symbol == 'I':
                references.append(("I-I", atom_counts['I'] / 2.0))
            elif symbol == 'F':
                references.append(("F-F", atom_counts['F'] / 2.0))
            elif symbol == 'Br':
                references.append(("Br-Br", atom_counts['Br'] / 2.0))
            else:
                references.append((f"[{symbol}]", atom_counts[f"{symbol}"]))

    return references


def calculate_free_energy_formation(mol,
                                    orca_path=None,
                                    xc='r2SCAN-3c',
                                    basis_set='def2-QZVP',
                                    tight_opt=False,
                                    tight_scf=False,
                                    f_solv=False,
                                    f_disp=False,
                                    n_procs=10,
                                    use_ccsd=False,
                                    debug=False):
    """
    Calculate the Gibbs free energy of formation for a molecule.

    This function computes the Gibbs free energy of formation for a molecule
    represented by an RDKit `Mol` object. It calculates the free energy, enthalpy,
    and entropy of the molecule and its reference molecules, and determines the
    formation energy by subtracting the reference contributions.

    Parameters:
    -----------
    mol : rdkit.Chem.Mol
        An RDKit molecule object representing the molecule.
    orca_path : str, optional
        Path to the ORCA executable. If None, it will attempt to read from the environment variable 'ORCA_PATH'.
    xc : str, optional
        Exchange-correlation functional to use. Default is 'r2SCAN-3c'.
    basis_set : str, optional
        Basis set to use for the calculation. Default is 'def2-QZVP'.
    tight_opt : bool, optional
        Whether to use tight geometry optimisation. Default is False.
    tight_scf : bool, optional
        Whether to use tight SCF convergence criteria. Default is False.
    f_solv : bool, optional
        Whether to include solvent effects in the calculation. Default is False.
    f_disp : bool, optional
        Whether to include dispersion corrections in the calculation. Default is False.
    n_procs : int, optional
        Number of processors to use for the calculation. Default is 10.
    use_ccsd : bool, optional
        Whether to use CCSD energy calculations. Default is False.
    debug : bool, optional
        Whether to print debug information during the calculation. Default is False.

    Returns:
    --------
    tuple
        A tuple containing:
        - d_free : float
            The Gibbs free energy of formation in eV.
        - d_enthalpy : float
            The enthalpy of formation in eV.
        - d_entropy : float
            The entropy correction of formation in eV.

    Raises:
    -------
    ValueError
        If the calculation fails or the molecule contains unsupported elements.
    """
    mol = Chem.AddHs(mol)  # Add explicit hydrogens to the molecule.
    atoms, charge, multiplicity = mol_to_atoms(mol), get_charge(mol), get_spin_multiplicity(
        mol)  # Convert to ASE Atoms and calculate properties.
    free, enthalpy, entropy = calculate_free_energy(atoms,
                                                    charge=charge,
                                                    multiplicity=multiplicity,
                                                    orca_path=orca_path,
                                                    xc=xc,
                                                    basis_set=basis_set,
                                                    tight_opt=tight_opt,
                                                    tight_scf=tight_scf,
                                                    f_solv=f_solv,
                                                    f_disp=f_disp,
                                                    n_procs=n_procs,
                                                    use_ccsd=use_ccsd)  # Calculate free energy, enthalpy, and entropy.
    if debug:
        print(f"Mol Free: {free}, Enthalpy: {enthalpy}, Entropy: {entropy}", flush=True)

    free_atoms = 0.0
    enthalpy_atoms = 0.0
    entropy_atoms = 0.0
    # Get the formation references
    references = get_formation_references(mol)  # Generate reference molecules for formation energy calculation.
    # Loop over the references and calculate the free energy
    for ref_smi, ref_count in references:
        ref_atoms, ref_charge, ref_multiplicity = smi_to_atoms(ref_smi)  # Convert reference SMILES to ASE Atoms.
        ref_free, ref_enthalpy, ref_entropy = calculate_free_energy(ref_atoms,
                                                                    charge=charge,
                                                                    multiplicity=multiplicity,
                                                                    orca_path=orca_path,
                                                                    xc=xc,
                                                                    basis_set=basis_set,
                                                                    tight_opt=tight_opt,
                                                                    tight_scf=tight_scf,
                                                                    f_solv=f_solv,
                                                                    f_disp=f_disp,
                                                                    n_procs=n_procs,
                                                                    use_ccsd=use_ccsd)  # Calculate free energy for reference molecules.
        if debug:
            print(f"Ref: {ref_smi}, count: {ref_count}, charge: {ref_charge}, multi: {ref_multiplicity}", flush=True)
            print(f"Free: {ref_free}, Enthalpy: {ref_enthalpy}, Entropy: {ref_entropy}", flush=True)
        free_atoms += ref_free * ref_count  # Accumulate free energy contributions from reference molecules.
        enthalpy_atoms += ref_enthalpy * ref_count  # Accumulate enthalpy contributions from reference molecules.
        entropy_atoms += ref_entropy * ref_count  # Accumulate entropy contributions from reference molecules.
    d_free = free - free_atoms  # Calculate the Gibbs free energy of formation.
    d_enthalpy = enthalpy - enthalpy_atoms  # Calculate the enthalpy of formation.
    d_entropy = entropy - entropy_atoms  # Calculate the entropy correction of formation.
    if debug:
        print(f"Atoms Free: {free_atoms}, Enthalpy: {enthalpy_atoms}, Entropy: {entropy_atoms}", flush=True)
        print(f"Deltas Free: {d_free}, Enthalpy: {d_enthalpy}, Entropy: {d_entropy}", flush=True)
    return d_free, d_enthalpy, d_entropy  # Return the formation energy components.


def extract_conformer_info(filepath: Union[str, Path]) -> pd.DataFrame:
    """
    Extract conformer information from an ORCA output file.

    This function reads an ORCA output file and parses the ensemble table to extract
    conformer data, including conformer index, energy, and percentage of the total.

    Parameters:
    -----------
    filepath : Union[str, Path]
        Path to the ORCA output file containing the ensemble table.

    Returns:
    --------
    pd.DataFrame
        A pandas DataFrame containing the following columns:
        - 'Conformer': Conformer index (int).
        - 'Energy_kcal_mol': Energy in kcal/mol (float).
        - 'Percent_total': Percentage of the total (float).

    Raises:
    -------
    ValueError
        If the ensemble table cannot be located in the file.
    """
    # Compile a regex pattern to match a data line in the ensemble table
    line_pat = re.compile(
        r"""^\s*
            (?P<conformer>\d+)\s+          # integer index
            (?P<energy>-?\d+\.\d+)\s+      # energy in kcal/mol
            \d+\s+                         # degeneracy (ignored)
            (?P<ptotal>\d+\.\d+)\s+        # % total
            \d+\.\d+\s*?$                  # % cumulative (ignored)
        """,
        re.VERBOSE,
    )

    # Compile a regex pattern to locate the table header
    header_pat = re.compile(r"Conformer\s+Energy.*% total", re.I)

    # Initialize variables for parsing
    rows = []
    in_table = False

    # Open the file and read its contents
    with open(filepath, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            # Check for the table header to start reading data
            if not in_table and header_pat.search(line):
                in_table = True  # Start reading on the next lines
                continue

            if in_table:
                # Stop reading when the table ends
                if line.strip() == "" or line.strip().startswith("Conformers"):
                    break
                # Match a data line and extract values
                m = line_pat.match(line)
                if m:
                    rows.append(
                        (
                            int(m["conformer"]),
                            float(m["energy"]),
                            float(m["ptotal"]),
                        )
                    )

    # Raise an error if no data was found
    if not rows:
        raise ValueError(
            "Could not locate ensemble table. Check that the file is complete."
        )

    # Return the extracted data as a pandas DataFrame
    return pd.DataFrame(
        rows, columns=["Conformer", "Energy_kcal_mol", "Percent_total"]
    )


def calculate_goat(atoms,
                   charge=0,
                   multiplicity=1,
                   orca_path=None,
                   n_procs=10):
    """
    Perform a GOAT (Global Optimization of Atomic Topologies) calculation using ORCA.

    This function sets up and executes a GOAT calculation to optimize molecular conformers
    and extract conformer information from the ORCA output file.

    Parameters:
    -----------
    atoms : ase.Atoms
        ASE Atoms object representing the molecule to be optimized.
    charge : int, optional
        Total charge of the molecule. Default is 0.
    multiplicity : int, optional
        Spin multiplicity of the molecule. Default is 1.
    orca_path : str, optional
        Path to the ORCA executable. If None, it will attempt to read from the environment variable 'ORCA_PATH'.
    n_procs : int, optional
        Number of processors to use for the calculation. Default is 10.

    Returns:
    --------
    tuple
        - atoms : list of ase.Atoms
            List of ASE Atoms objects representing the optimized conformers.
        - df : pandas.DataFrame
            DataFrame containing conformer information, including:
            - 'Conformer': Conformer index (int).
            - 'Energy_kcal_mol': Energy in kcal/mol (float).
            - 'Percent_total': Percentage of the total (float).
    """
    # Determine the ORCA path
    if orca_path is None:
        # Try to read the path from the environment variable
        orca_path = os.environ.get('ORCA_PATH')
    else:
        # Convert the provided path to an absolute path
        orca_path = os.path.abspath(orca_path)

    # Create an ORCA profile with the specified command
    profile = OrcaProfile(command=orca_path)

    # Configure the number of processors
    if n_procs > 1:
        inpt_procs = '%pal nprocs {} end'.format(n_procs)
    else:
        inpt_procs = ''

    # Create a temporary working directory
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create and configure the ORCA calculator object
        calc = ORCA(
            profile=profile,
            charge=charge,
            mult=multiplicity,
            directory=temp_dir,
            orcasimpleinput='GOAT XTB',
            orcablocks=inpt_procs
        )
        # Assign the calculator to the ASE Atoms object
        atoms.calc = calc

        # Trigger the calculation to optimize the geometry
        _ = atoms.get_potential_energy()

        # Define paths for the output files
        xyz_file = os.path.join(temp_dir, "orca.finalensemble.xyz")  # Path to the final ensemble file
        orca_file = os.path.join(temp_dir, "orca.out")  # Path to the ORCA output file

        # Extract conformer information from the ORCA output file
        df = extract_conformer_info(orca_file)

        # Read the optimized conformers from the ensemble file
        atoms = read(xyz_file, format="xyz", index=':')

        # Return the optimized conformers and conformer information
        return atoms, df


def get_symmetry_number(atoms, tolerance=0.3):
    if getattr(atoms, "pbc", None) is not None and any(atoms.pbc):
        # Molecules should be non-periodic for a correct point-group analysis
        atoms = atoms.copy()
        atoms.pbc = (False, False, False)

    # Convert ASE -> pymatgen Molecule
    pmg_mol = AseAtomsAdaptor.get_molecule(atoms)

    # Analyze point group
    pga = PointGroupAnalyzer(pmg_mol, tolerance=tolerance)
    # e.g., "C2v", "D6h", "Td", "Oh", "Cinfv", "Dinfh"
    symbol = getattr(pga, "sch_symbol", None) or pga.get_pointgroup()

    # Normalize common infinite-group spellings
    s = str(symbol).replace("∞", "inf").replace("*", "inf").replace("Inf", "inf")

    # Linear groups
    if s.lower() in {"cinfv", "cinfv"}:
        return 1
    if s.lower() in {"dinfh", "dinfh"}:
        return 2

    # Trivial groups
    if s in {"C1", "Ci", "Cs"}:
        return 1

    # Polyhedral (Platonic) groups
    if s in {"T", "Td", "Th"}:
        return 12
    if s in {"O", "Oh"}:
        return 24
    if s in {"I", "Ih"}:
        return 60

    # Cyclic / dihedral / improper groups with finite n
    # Match patterns like C3, C3v, C3h
    m_c = re.fullmatch(r"C(\d+)(?:[vh])?", s)
    if m_c:
        return int(m_c.group(1))

    # Match patterns like D5, D5d, D5h
    m_d = re.fullmatch(r"D(\d+)(?:[dh])?", s)
    if m_d:
        return 2 * int(m_d.group(1))

    # Improper S_{2n} groups -> sigma = n (only even-ordered S groups exist)
    m_s = re.fullmatch(r"S(\d+)", s)
    if m_s:
        order = int(m_s.group(1))
        if order % 2 == 0:
            return order // 2  # S_{2n} -> n
        # Odd S_n shouldn't occur for rigid molecules; fall through just in case.

    # As a safe fallback, count proper operations from the full symmetry-ops set.
    # (This works for all *finite* groups).
    try:
        ops = pga.get_symmetry_operations()
        proper = [op for op in ops if np.linalg.det(op.rotation_matrix) > 0]
        if len(proper) > 0:
            return int(len(proper))
    except Exception:
        pass

    # Last resort
    return 1


def classify_geometry(atoms, tol_ratio=1e-3, tol_abs=1e-6):
    n = len(atoms)
    if n == 1:
        return 'monatomic'
    if n == 2:
        return 'linear'

    # Principal moments of inertia (amu·Å^2)
    inertia = np.sort(np.asarray(atoms.get_moments_of_inertia()))
    in_min, in_mid, in_max = inertia

    # Characteristic inertia scale ~ M * <r^2> (helps set an absolute floor)
    masses = atoms.get_masses()
    pos = atoms.get_positions()  # assumes not wrapping across PBC
    com = np.average(pos, axis=0, weights=masses)
    r = pos - com
    r2_weighted_mean = np.sum(masses * np.sum(r ** 2, axis=1)) / np.sum(masses)
    in_char = np.sum(masses) * r2_weighted_mean

    is_linear = (in_min / (in_max if in_max > 0 else 1.0) < tol_ratio) and (in_min < tol_abs * in_char)
    return 'linear' if is_linear else 'nonlinear'


def multiplicity_to_total_spin(multiplicity, check_integer=True, tol=1e-8):
    m = float(multiplicity)
    if not math.isfinite(m) or m < 1:
        raise ValueError("Multiplicity must be a finite number >= 1.")

    if check_integer:
        if abs(m - round(m)) > tol:
            raise ValueError("Multiplicity must be an integer (1, 2, 3, ...).")
        m = round(m)

    return (m - 1.0) / 2.0


def free_energy_mace(atoms,
                     charge=0,
                     multiplicity=1,
                     optimise=True,
                     f_max=0.01,
                     temperature=298.15,
                     pressure=101325.0,
                     calc_model='extra_large',
                     calc_device="cuda"):
    calc = mace_omol(model=calc_model, device=calc_device)
    atoms.calc = calc
    atoms.info["charge"] = charge
    atoms.info["spin"] = multiplicity
    if optimise:
        BFGS(atoms,
             logfile=None,
             trajectory=None).run(fmax=f_max)

    energy = atoms.get_potential_energy()

    with tempfile.TemporaryDirectory() as run_dir:
        vib = Vibrations(atoms, name=run_dir)
        vib.run()
        vib_energies = vib.get_energies()
        if os.path.exists(run_dir):
            shutil.rmtree(run_dir)

        thermo = IdealGasThermo(
            vib_energies=vib_energies,
            potentialenergy=energy,
            atoms=atoms,
            geometry=classify_geometry(atoms),
            symmetrynumber=get_symmetry_number(atoms),
            spin=multiplicity_to_total_spin(multiplicity),
        )
        free_energy = thermo.get_gibbs_energy(temperature=temperature,
                                              pressure=pressure)
        free_enthalpy = thermo.get_enthalpy(temperature=temperature)
        free_entropy = thermo.get_entropy(temperature=temperature,
                                          pressure=pressure)
        return free_energy, free_enthalpy, free_entropy


def calculate_free_energy_formation_mace(mol,
                                         optimise=True,
                                         f_max=0.01,
                                         temperature=298.15,
                                         pressure=101325.0,
                                         calc_model='extra_large',
                                         calc_device="cuda"):
    mol = Chem.AddHs(mol)
    atoms = mol_to_atoms(mol)
    charge = get_charge(mol)
    multiplicity = get_spin_multiplicity(mol)
    free, enthalpy, entropy = free_energy_mace(atoms,
                                               charge=charge,
                                               multiplicity=multiplicity,
                                               optimise=optimise,
                                               f_max=f_max,
                                               temperature=temperature,
                                               pressure=pressure,
                                               calc_model=calc_model,
                                               calc_device=calc_device)
    free_atoms = 0.0
    enthalpy_atoms = 0.0
    entropy_atoms = 0.0
    # Get the formation references
    references = get_formation_references(mol)
    for ref_smi, ref_count in references:
        ref_atoms, ref_charge, ref_multiplicity = smi_to_atoms(ref_smi)
        ref_free, ref_enthalpy, ref_entropy = free_energy_mace(ref_atoms,
                                                               charge=ref_charge,
                                                               multiplicity=ref_multiplicity,
                                                               optimise=optimise,
                                                               f_max=f_max,
                                                               temperature=temperature,
                                                               pressure=pressure,
                                                               calc_model=calc_model,
                                                               calc_device=calc_device)
        free_atoms += ref_free * ref_count
        enthalpy_atoms += ref_enthalpy * ref_count
        entropy_atoms += ref_entropy * ref_count

    d_free = free - free_atoms
    d_enthalpy = enthalpy - enthalpy_atoms
    d_entropy = entropy - entropy_atoms
    return d_free, d_enthalpy, d_entropy
