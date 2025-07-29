import os
import re
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Union

import pandas as pd
from ase.atoms import Atoms
from ase.calculators.orca import ORCA
from ase.calculators.orca import OrcaProfile
from ase.io import read
from ase.units import Hartree
from rdkit import Chem as Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdchem import Mol

from .tools_mols import standardize_mol


def smi_to_atoms(smiles: str) -> tuple[Atoms, int, int]:
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    mol = Chem.AddHs(mol)
    mol = standardize_mol(mol)
    if mol is None:
        raise ValueError(f"Failed to parse SMILES string: {smiles}")

    return mol_to_atoms(mol), get_charge(mol), get_spin_multiplicity(mol)


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


def get_spin_multiplicity(mol: Mol) -> int:
    """
    Calculate the spin multiplicity of a molecule based on the number of radical electrons.

    Spin multiplicity = 2 S + 1, where S is the total spin quantum number.
    For radical electrons, S = n/2 where n is the number of unpaired electrons.

    Parameters:
    -----------
    mol : rdkit.Chem.rdchem.Mol
        An RDKit molecule object

    Returns:
    --------
    int
        The spin multiplicity of the molecule (1 for singlet, 2 for doublet, etc.)
    """
    # Get the number of radical electrons
    num_radical_electrons = Descriptors.NumRadicalElectrons(mol)

    # Calculate spin multiplicity: 2 S + 1 = n + 1, where n is the number of unpaired electrons
    multiplicity = num_radical_electrons + 1

    return multiplicity


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

    Args:
        filename (str): Path to the IR frequency output file

    Returns:
        pd.DataFrame: DataFrame containing Mode, Frequency (cm^-1), Epsilon, and Intensity
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
            mode = int(match.group(1))
            freq = float(match.group(2))
            eps = float(match.group(3))
            intensity = float(match.group(4))
            data.append([mode, freq, eps, intensity])

    # Create DataFrame
    df = pd.DataFrame(data, columns=['Mode', 'Frequency (cm^-1)', 'Epsilon', 'Intensity (km/mol)'])

    return df


def load_raman_data(filename):
    """
    Load Raman spectrum data from ORCA output file into a pandas DataFrame.

    Args:
        filename (str): Path to the Raman frequency output file

    Returns:
        pd.DataFrame: DataFrame containing Mode, Frequency, Activity, and Depolarisation
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
            mode = int(match.group(1))
            freq = float(match.group(2))
            activity = float(match.group(3))
            depolarization = float(match.group(4))
            data.append([mode, freq, activity, depolarization])

    # Create DataFrame
    df = pd.DataFrame(data, columns=['Mode', 'Frequency (cm^-1)', 'Intensity (km/mol)', 'Depolarization'])

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
        - 'Epsilon': Default value set to 1.0 (float).
        - 'Intensity (km/mol)': Default value set to 1.0 (float).

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
            eps = 1.0  # Default epsilon value
            intensity = 1.0  # Default intensity value
            data.append([mode, freq, eps, intensity])

    # Create DataFrame
    df = pd.DataFrame(data, columns=['Mode', 'Frequency (cm^-1)', 'Epsilon', 'Intensity (km/mol)'])

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

        # Perform the calculation (this will write the output to the ORCA file)
        _ = atoms.get_potential_energy()

        # Load IR spectrum data
        data_ir = load_ir_data(orca_file)

        # Load Raman spectrum data
        data_raman = load_raman_data(orca_file)

        # Load vibrational spectrum data
        data_vib = load_vib_data(orca_file)

        return data_ir, data_raman, data_vib


def calculate_ccsd_energy(atoms,
                          charge=0,
                          multiplicity=1,
                          orca_path=None,
                          basis_set='def2-TZVPP',
                          n_procs=10):
    # If no ORCA path is provided, try to read it from the environment variable
    if orca_path is None:
        orca_path = os.environ.get('ORCA_PATH')
    else:
        orca_path = os.path.abspath(orca_path)

    # Create a temporary directory for the ORCA calculation
    with tempfile.TemporaryDirectory() as temp_dir:
        # temp_dir= 'tmp'
        # os.makedirs(temp_dir, exist_ok=True)
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
    with open(orca_file, 'r') as f:
        for line in reversed(f.readlines()):
            if term in line:
                return float(line.split(splitter)[-1].split('Eh')[0]) * Hartree
        return None


def calculate_free_energy(atoms,
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
                          use_ccsd=False,
                          ccsd_energy=None):
    orca_path = os.path.abspath(orca_path or os.getenv('ORCA_PATH', 'orca'))
    opt_flag = 'TIGHTOPT' if tight_opt else 'OPT'
    scf_flag = 'TIGHTSCF' if tight_scf else ''
    calc_extra = f'{opt_flag} {scf_flag} FREQ'.strip()
    if use_ccsd and ccsd_energy is None:
        ccsd_energy = calculate_ccsd_energy(atoms,
                                            orca_path=orca_path,
                                            charge=charge,
                                            multiplicity=multiplicity,
                                            n_procs=n_procs)
        if ccsd_energy is None:
            raise ValueError("CCSD energy calculation failed. Please check the ORCA setup.")
    with tempfile.TemporaryDirectory() as temp_dir:
        orca_file = os.path.join(temp_dir, 'orca.out')
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
        atoms.calc = calc
        _ = atoms.get_potential_energy()

        entropy = grab_value(orca_file, 'Total entropy correction', '...')

        if use_ccsd:
            g_e_ele = grab_value(orca_file, 'G-E(el)', '...')
            g_e_solv = grab_value(orca_file, 'Free-energy (cav+disp)', ':') if f_solv else 0.0
            energy = ccsd_energy + g_e_ele + g_e_solv
        else:
            energy = grab_value(orca_file, 'Final Gibbs free energy', '...')

        return energy, energy - entropy, entropy


def list_to_str(lst):
    lst = [str(item) for item in lst]
    return ', '.join(lst)


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

        # Load the optimised geometry from the ORCA output file
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
    # Get elemental composition
    atom_counts = defaultdict(int)
    symbols = set()
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        atom_counts[symbol] += 1.0
        symbols.add(symbol)

    supported_set = {'H', 'Ba', 'Cu', 'S', 'K', 'Cl', 'Rb', 'B', 'N', 'Se', 'Te', 'O', 'Fe', 'Co', 'Mg', 'Ge', 'I',
                     'Tl', 'Pt', 'Xe', 'Zn', 'Gd', 'Cd', 'C', 'Al', 'F', 'Li', 'Ca', 'Ni', 'Th', 'Sr', 'Sn', 'Au', 'Ag',
                     'V', 'Pu', 'Sb', 'Mn', 'Cr', 'Mo', 'Ra', 'Pb', 'Bi', 'As', 'Hg', 'Si', 'Br', 'Rn', 'P', 'W', 'Be',
                     'Na'}
    if not symbols.issubset(supported_set):
        raise ValueError(f"Unsupported elements in the molecule: {symbols - supported_set}")

    # Standard reference molecules
    references = []
    # loop over the atom_counts dictionary
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
                                    use_ccsd=False):
    mol = Chem.AddHs(mol)

    atoms, charge, multiplicity = mol_to_atoms(mol), get_charge(mol), get_spin_multiplicity(mol)
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
                                                    use_ccsd=use_ccsd)
    print(f"Mol Free: {free}, Enthalpy: {enthalpy}, Entropy: {entropy}", flush=True)

    free_atoms = 0.0
    enthalpy_atoms = 0.0
    entropy_atoms = 0.0

    # Get the formation references
    references = get_formation_references(mol)
    # Loop over the references and calculate the free energy
    for ref_smi, ref_count in references:
        ref_atoms, ref_charge, ref_multiplicity = smi_to_atoms(ref_smi)
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
                                                                    use_ccsd=use_ccsd)
        print(f"Reference: {ref_smi}, Count: {ref_count}", flush=True)
        print(f"Free: {ref_free}, Enthalpy: {ref_enthalpy}, Entropy: {ref_entropy}", flush=True)
        free_atoms += ref_free * ref_count
        enthalpy_atoms += ref_enthalpy * ref_count
        entropy_atoms += ref_entropy * ref_count
    print(f"Atoms Free: {free_atoms}, Enthalpy: {enthalpy_atoms}, Entropy: {entropy_atoms}", flush=True)
    d_free = free - free_atoms
    d_enthalpy = enthalpy - enthalpy_atoms
    d_entropy = entropy - entropy_atoms
    print(f"Deltas Free: {d_free}, Enthalpy: {d_enthalpy}, Entropy: {d_entropy}", flush=True)
    print(f"Half values Free: {d_free / 2.0}, Enthalpy: {d_enthalpy / 2.0}, Entropy: {d_entropy / 2.0}", flush=True)
    return d_free, d_enthalpy, d_entropy


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
