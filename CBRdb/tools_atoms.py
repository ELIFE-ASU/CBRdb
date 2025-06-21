import os
import tempfile

from ase.atoms import Atoms
from ase.calculators.orca import ORCA
from ase.calculators.orca import OrcaProfile
from ase.io import read
from ase.units import Hartree
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdchem import Mol

from .tools_mols import standardize_mol


def smi_to_atoms(smiles: str) -> Atoms:
    """
    Convert a SMILES string to an ASE Atoms object via an SDF file.

    Parameters:
    -----------
    smiles : str
        SMILES string representing the molecule.

    Returns:
    --------
    Atoms
        Atoms object representing the molecule.

    Raises:
    -------
    ValueError
        If SMILES parsing fails.
    """
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    mol = standardize_mol(mol)
    if mol is None:
        raise ValueError(f"Failed to parse SMILES string: {smiles}")

    return mol_to_atoms(mol)


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
                     xc='B3LYP',
                     charge=0,
                     multiplicity=1,
                     basis_set='6-311G',  # 'cc-pVDZ', '6-31+G(d,p)',
                     nprocs=1,
                     f_solv=True,
                     f_disp=True,
                     atom_list=None,
                     calc_extra=None,
                     blocks_extra=None,
                     scf_option=None):
    if orca_path is None:
        # try and read the path from the environment
        orca_path = os.environ.get('ORCA_PATH')
    if directory is None:
        directory = os.path.join(tempfile.mkdtemp(), 'orca')

    profile = OrcaProfile(command=orca_path)

    if nprocs > 1:
        inpt_procs = '%pal nprocs {} end'.format(nprocs)
    else:
        inpt_procs = ''

    if f_solv is not None and f_solv is not False:
        if f_solv:
            f_solv = 'WATER'
        inpt_solv = '''
        %CPCM SMD TRUE
            SMDSOLVENT "{}"
        END'''.format(f_solv)
    else:
        inpt_solv = ''

    if f_disp is None or f_disp is False:
        inpt_disp = ''
    else:
        if f_disp:
            f_disp = 'D4'
        inpt_disp = f_disp

    if atom_list is not None and calc_type == 'QM/XTB2':
        inpt_xtb = '''
        %QMMM QMATOMS {{}} END END
        '''.format(str(atom_list).strip('[').strip(']'))
    else:
        inpt_xtb = ''

    if blocks_extra is None:
        blocks_extra = ''

    inpt_blocks = inpt_procs + inpt_solv + blocks_extra

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

    # Add the scf option
    if scf_option is not None:
        inpt_simple += ' ' + scf_option

    # Add the extra options
    if calc_extra is not None:
        inpt_simple += ' ' + calc_extra

    calc = ORCA(
        profile=profile,
        charge=charge,
        mult=multiplicity,
        directory=directory,
        orcasimpleinput=inpt_simple,
        orcablocks=inpt_blocks
    )
    return calc


def calculate_free_energy(mol,
                          orca_path=None,
                          xc='r2SCAN-3c',  # wB97X def2-TZVP def2/J RIJCOSX
                          basis_set='def2-QZVP',  # def2-QZVP H–Rn aug-cc-pVTZ H–Ar, Sc–Kr, Ag, Au
                          calc_extra='TIGHTOPT FREQ',
                          f_solv=False,
                          f_disp=False,
                          nprocs=10):
    # Good settings r2SCAN-3c def2-QZVP
    # Expensive settings wB97X aug-cc-pVTZ

    if orca_path is None:
        # try and read the path from the environment
        orca_path = os.environ.get('ORCA_PATH')
    else:
        orca_path = os.path.abspath(orca_path)

    # Sanitise the molecule
    mol = Chem.AddHs(mol)
    Chem.SanitizeMol(mol)

    # Get the formal charge of the molecule
    charge = get_charge(mol)

    # Get the spin multiplicity (note that this assumes the molecule has been correctly assigned a spin state)
    multiplicity = get_spin_multiplicity(mol)

    # Convert the RDKit molecule to an ASE Atoms object
    atoms = mol_to_atoms(mol)

    # Create a working directory
    temp_dir = tempfile.mkdtemp()
    orca_file = os.path.join(temp_dir, 'orca.out')

    # Set up the ORCA calculator
    calc = orca_calc_preset(orca_path=orca_path,
                            directory=temp_dir,
                            charge=charge,
                            multiplicity=multiplicity,
                            xc=xc,
                            basis_set=basis_set,
                            nprocs=nprocs,
                            f_solv=f_solv,
                            f_disp=f_disp,
                            calc_extra=calc_extra,
                            )

    # Attach the calculator to the atom object
    atoms.calc = calc

    # Perform the calculation
    _ = atoms.get_potential_energy()
    # Read the output file to extract the free energy
    with open(orca_file, 'r') as f:
        for line in reversed(f.readlines()):
            if 'Final Gibbs free energy' in line:
                energy = float(line.split('...')[-1].split('Eh')[0])
                # Convert from Hartree to eV
                energy *= Hartree
                break

    return energy
