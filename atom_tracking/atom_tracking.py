"""============================================================================
                        CHEMICAL REACTION ATOM TRACKING SYSTEM
===============================================================================

DESCRIPTION:
This script performs automated atom tracking in chemical reactions using 
advanced atom mapping algorithms. It identifies specific atoms in substrate 
molecules and traces their fate through chemical transformations to determine 
which product molecules contain the tracked atoms.

PURPOSE:
- Track the fate of specific atoms through chemical reactions
- Handle complex reaction equations with coefficients and cofactors
- Provide robust atom mapping with multiple fallback strategies
- Support parallel processing for large reaction datasets
- Generate visualizations of mapped reactions

WORKFLOW OVERVIEW:
1. REACTION PARSING: Parse reaction equations, remove coefficients, handle cofactors
2. MOLECULE SELECTION: Choose substrate molecules based on various criteria
3. ATOM SELECTION: Select specific atoms for tracking (highest degree by default)
4. ATOM MAPPING: Use RXNMapper (primary) with localmapper fallback for token limits
5. ATOM TRACKING: Trace atoms through the reaction using intelligent detection
6. RESULT OUTPUT: Save tracking results and handle errors gracefully

KEY FEATURES:
- Dual mapping system: RXNMapper → localmapper fallback for token limit errors
- Intelligent atom map number detection with multiple strategies
- Comprehensive cofactor filtering (static lists + dynamic detection)
- Robust error handling and validation throughout
- Parallel processing support with file locking
- Optional reaction visualization
- Detailed debugging output for troubleshooting

TECHNICAL COMPONENTS:
- RDKit: Molecular processing and SMILES handling
- RXNMapper: Primary atom mapping (512 token limit)
- localmapper: Fallback mapping for large reactions
- pandas: Data processing and CSV operations
- FileLock: Thread-safe file operations

COMMAND LINE USAGE:
python AtomTracking_NEW.py --row_index <index> [OPTIONS]

Required Arguments:
  --row_index INT        : Row index of reaction to process

Optional Arguments:
  --static_cofactors     : Remove cofactors using predefined list
  --dynamic_cofactors    : Remove compounds present on both reaction sides
  --highest_CCC          : Select substrate with most chiral centers
  --lowest_CCC           : Select substrate with fewest chiral centers
  --visualize            : Generate reaction visualization

INPUT FILES:
- Last_Missing.csv      : Reaction equations to process
- CBRdb_C.csv          : Compound reference with SMILES structures
- CBRdb_CanonicalIndex.csv : Canonical SMILES index (auto-generated)

OUTPUT FILES:
- AtomTracking_P3.csv   : Successful tracking results
- Error Reactions.csv   : Failed reactions with error details

STEP-BY-STEP PROCESS:
1. Load and validate input files
2. Parse reaction equation and extract compound IDs
3. Remove numerical/variable coefficients from equations
4. Apply cofactor filtering if requested
5. Convert compound IDs to SMILES structures
6. Select target molecule and atom for tracking
7. Perform atom mapping with fallback strategies
8. Track atom through reaction products
9. Identify final product containing tracked atom
10. Save results or error information

ERROR HANDLING:
- Empty/invalid SMILES validation
- Token limit errors with automatic fallback
- Missing compound ID resolution
- Mapping failure recovery strategies
- Comprehensive error logging

PARALLEL PROCESSING:
The script is designed for parallel execution across multiple reaction rows.
File locking ensures safe concurrent writes to output files.

============================================================================"""

import argparse
import os
import sys

import matplotlib.pyplot as plt
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rxnmapper import RXNMapper

# --- Import localmapper for fallback mapping ---
try:
    from localmapper import localmapper

    LOCALMAPPER_AVAILABLE = True
    print("[INFO] localmapper package is available for fallback mapping")
except ImportError:
    LOCALMAPPER_AVAILABLE = False
    print("[WARNING] localmapper package not available - no fallback mapping will be used")
    print("[INFO] To install localmapper, run: pip install localmapper")
# --- File locking for safe parallel writes ---
try:
    from filelock import FileLock
except ImportError:
    import subprocess

    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'filelock'])
    from filelock import FileLock


# -------------------- Helper Functions --------------------
def canonicalize_smiles(smiles):
    """
    Convert a SMILES string to its canonical form, removing all atom map numbers.
    Returns None if input is invalid or cannot be parsed.
    """
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def smiles_to_mols(smiles):
    """
    Convert a dot-separated SMILES string into a list of RDKit Mol objects.
    Only valid fragments are returned as Mol objects.
    """
    if not smiles or not smiles.strip():
        print(f"[WARNING] Empty or whitespace-only SMILES string provided to smiles_to_mols")
        return []

    # Split by '.' and filter out empty strings
    fragments = [smi.strip() for smi in smiles.strip().split('.') if smi.strip()]

    if not fragments:
        print(f"[WARNING] No valid SMILES fragments found in: '{smiles}'")
        return []

    mols = []
    for smi in fragments:
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                mols.append(mol)
            else:
                print(f"[WARNING] Invalid SMILES fragment: '{smi}'")
        except Exception as e:
            print(f"[ERROR] Failed to parse SMILES fragment '{smi}': {e}")

    return mols


def label_atom(mol, atom_idx, map_num):
    """
    Label a specific atom in a molecule with a given atom map number.
    Returns a new RDKit Mol object with the label applied.
    """
    mol = Chem.RWMol(mol)
    if 0 <= atom_idx < mol.GetNumAtoms():
        mol.GetAtomWithIdx(atom_idx).SetAtomMapNum(map_num)
    return mol


def trace_atom(atom_map_num, mols):
    """
    Search a list of molecules for an atom with a specific atom map number.
    Returns the first molecule containing the atom, or None if not found.
    """
    for i, mol in enumerate(mols):
        for j, atom in enumerate(mol.GetAtoms()):
            if atom.GetAtomMapNum() == atom_map_num:
                print(f"[DEBUG] Found atom map {atom_map_num} in molecule {i}, atom {j} ({atom.GetSymbol()})")
                return mol
    print(f"[DEBUG] Atom map {atom_map_num} not found in any of the {len(mols)} molecules")
    return None


def remove_atom_map_numbers(mol):
    """
    Remove all atom map numbers from a molecule and return its canonical SMILES.
    Returns None if input is None.
    """
    if mol is None:
        return None
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def visualize_mapped_rxn(mapped_rxn_smiles, reaction_id, save_images=False):
    """
    Visualize the mapped reaction SMILES. If save_images is True, saves the image
    in a folder named 'Reactions Visualizations' in the same directory as the output CSV,
    using the reaction ID as the filename.
    """
    try:
        rxn = Chem.rdChemReactions.ReactionFromSmarts(mapped_rxn_smiles, useSmiles=True)
        img = Draw.ReactionToImage(rxn, subImgSize=(600, 600))
        vis_dir = os.path.join(os.path.dirname(output_file), "Reactions Visualizations")
        os.makedirs(vis_dir, exist_ok=True)
        img_path = os.path.join(vis_dir, f"{reaction_id}.png")
        plt.imshow(img)
        plt.axis('off')
        plt.title(f"Reaction ID: {reaction_id}")
        plt.tight_layout()
        plt.savefig(img_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"[INFO] Saved reaction visualization to {img_path}")
    except Exception as e:
        print(f"[ERROR] Could not visualize reaction {reaction_id}: {e}")


# lookup_compound_id must be defined after df_canonical_index is available
lookup_compound_id = None

# -------------------- File Locator --------------------
script_dir = os.path.dirname(os.path.abspath(__file__))


def find_file(filename):
    """
    Search for a file in the current script directory and up to two parent directories.
    Returns the absolute path if found, otherwise raises FileNotFoundError.
    """
    for base in [script_dir, os.path.join(script_dir, ".."), os.path.join(script_dir, "..", "..")]:
        candidate = os.path.abspath(os.path.join(base, filename))
        if os.path.exists(candidate):
            return candidate
    raise FileNotFoundError(f"{filename} not found in current or parent directories.")


# ------------------- canonical indeces in one csv file -------------------------
def get_or_build_canonical_index(df_molecules, output_file):
    """
    Load a canonical index from CSV if it exists, otherwise build it from the molecule DataFrame.
    Each fragment of each molecule is canonicalized and indexed.
    Returns a DataFrame with columns: compound_id, original_smiles, fragment_smiles, canonical_smiles.
    """
    if os.path.exists(output_file):
        print(f"[INFO] Canonical index found at {output_file}, loading it.")
        return pd.read_csv(output_file)
    else:
        print(f"[INFO] Canonical index not found. Building it...")
        index = []
        for _, row in df_molecules.iterrows():
            cid = row['compound_id']
            smiles = row['smiles']
            if pd.isna(smiles):
                continue
            fragments = smiles.strip().split('.')
            for frag in fragments:
                mol = Chem.MolFromSmiles(frag)
                if mol:
                    canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
                    index.append({
                        'compound_id': cid,
                        'original_smiles': smiles,
                        'fragment_smiles': frag,
                        'canonical_smiles': canonical
                    })
                else:
                    print(f"[WARNING] Invalid fragment: {frag} for {cid}")
        df_index = pd.DataFrame(index)
        df_index.to_csv(output_file, index=False)
        print(f"[INFO] Canonical index saved to {output_file}")
        return df_index


# -------------------- Inputs --------------------

parser = argparse.ArgumentParser()
parser.add_argument('--row_index', type=int, required=True)
parser.add_argument('--static_cofactors', action='store_true', help='Ignores Cofactors based on a static list')
parser.add_argument('--dynamic_cofactors', action='store_true', help='Ignores Compounds present on both sides.')
parser.add_argument('--highest_CCC', action='store_true', help='Select substrate with highest number of chiral centers')
parser.add_argument('--lowest_CCC', action='store_true', help='Select substrate with lowest number of chiral centers')
parser.add_argument('--visualize', action='store_true', help='If set, visualize the mapped reaction for this row index')

args = parser.parse_args()

reactions_file = 'CBRdb_R.csv'
molecule_reference_file = 'CBRdb_C.csv'
canonical_index_file = os.path.join(script_dir, 'CBRdb_CanonicalIndex.csv')

output_file = os.path.join(script_dir, 'AtomTracking.csv')
error_output_file = os.path.join(script_dir, 'Error Reactions.csv')
lock_file = output_file + '.lock'
error_lock_file = error_output_file + '.lock'

try:
    reactions_file_path = find_file(reactions_file)
    molecule_ref_path = find_file(molecule_reference_file)
except FileNotFoundError as e:
    print(e)
    sys.exit(1)

try:
    df_input_reactions = pd.read_csv(reactions_file_path)
    df_molecules = pd.read_csv(molecule_ref_path)
except Exception as e:
    print(f"Error loading input files: {e}")
    sys.exit(1)

if args.row_index >= len(df_input_reactions):
    print(f"Row index {args.row_index} out of bounds.")
    sys.exit(1)

df_canonical_index = get_or_build_canonical_index(df_molecules, canonical_index_file)


# Now that df_canonical_index is available, define lookup_compound_id
def lookup_compound_id(canonical_smiles):
    """
    Look up the compound ID for a given canonical SMILES string in the canonical index DataFrame.
    Returns the compound ID if found, otherwise 'Unknown'.
    """
    match = df_canonical_index[df_canonical_index['canonical_smiles'] == canonical_smiles]
    if not match.empty:
        return match['compound_id'].values[0]
    return 'Unknown'


# ----------------------Selection-------------------------
selected_row = df_input_reactions.iloc[args.row_index]
selected_id = selected_row['id']
raw_reaction = selected_row['reaction']

rxn_mapper = RXNMapper()

# Initialize localmapper if available
local_mapper = None
if LOCALMAPPER_AVAILABLE:
    try:
        local_mapper = localmapper()
        print("[INFO] localmapper initialized successfully")
    except Exception as e:
        print(f"[WARNING] Failed to initialize localmapper: {e}")
        LOCALMAPPER_AVAILABLE = False


def perform_atom_mapping(reaction_smiles, reaction_id):
    """
    Perform atom mapping with RXNMapper first, fallback to localmapper if token limit exceeded.
    Returns tuple: (mapped_rxn, mapper_used, error_message)
    """
    # Try RXNMapper first
    try:
        print(f"[INFO] Attempting mapping with RXNMapper for reaction {reaction_id}")
        mapped_rxn = rxn_mapper.get_attention_guided_atom_maps([reaction_smiles])[0]['mapped_rxn']
        print(f"[INFO] RXNMapper successful for reaction {reaction_id}")
        return mapped_rxn, "RXNMapper", None
    except Exception as e:
        error_msg = str(e)
        print(f"[WARNING] RXNMapper failed for reaction {reaction_id}: {error_msg}")

        # Check if error is due to token limit
        token_limit_keywords = ["tokens", "token", "512", "should be at most", "token limit", "maximum length"]
        is_token_limit_error = any(keyword in error_msg.lower() for keyword in token_limit_keywords)

        if is_token_limit_error:
            print(f"[INFO] Token limit exceeded for reaction {reaction_id}, attempting localmapper fallback")

            if LOCALMAPPER_AVAILABLE and local_mapper is not None:
                try:
                    print(f"[INFO] Attempting mapping with localmapper for reaction {reaction_id}")
                    mapped_rxn = local_mapper.get_atom_map(reaction_smiles)
                    print(f"[INFO] localmapper successful for reaction {reaction_id}")
                    return mapped_rxn, "localmapper_fallback", None
                except Exception as local_e:
                    local_error_msg = str(local_e)
                    print(f"[ERROR] localmapper also failed for reaction {reaction_id}: {local_error_msg}")
                    return None, "both_failed", f"RXNMapper: {error_msg}; localmapper: {local_error_msg}"
            else:
                print(f"[ERROR] localmapper not available for fallback")
                return None, "rxnmapper_failed_no_fallback", f"RXNMapper: {error_msg}; localmapper not available"
        else:
            # For non-token-limit errors, don't try localmapper
            print(f"[ERROR] RXNMapper failed with non-token-limit error: {error_msg}")
            return None, "rxnmapper_failed_other", error_msg


# -------------------Cofactor list--------------------
cofactor_ids = [
    'C00002',  # Adenosine triphosphate (ATP)
    'C00003',  # Nicotinamide Adenine Dinucleotide (NAD+)(Vitamin B3)
    'C00004',  # Reduced Nicotinamide-adenine dinucleotide (NADH)(Vitamin B3)
    'C00005',  # Reduced nicotinamide adenine dinucleotide phosphate (NADPH)(Vitamin B3)
    'C00006',  # Nicotinamide adenine dinucleotide phosphate (NADP+)(Vitamin B3)
    'C00008',  # Adenosine diphosphate(ADP)
    'C00010',  # Coenzyme A (Vitamin B5)
    'C00016',  # Flavin adenine dinucleotide (FAD)(Vitamin B2)
    'C00018',  # Pyridoxal phosphate (Vitamin B6)
    'C00019',  # S-Adenosyl methionine (SAM)
    'C00020',  # Adenosine monophosphate (AMP)
    'C00032',  # Heme B
    'C00034',  # Manganese (Mn)
    'C00038',  # Zinc (Zn2+)
    'C00051',  # Glutathione (GSH)
    'C00053',  # 3'-Phosphoadenylyl sulfate (PAPS)
    'C00061',  # Flavin mononucleotide (FMN)(Vitamin B2)
    'C00063',  # Cytidine triphosphate (CTP)
    'C00068',  # Thiamin diphosphate (Vitamin B1)
    'C00072',  # Ascorbic acid (Vitamin C)
    'C00076',  # Calcium Ion (Ca2+)
    'C00101',  # Tetrahydrofolate (Vitamin B9)
    'C00113',  # Pyrroloquinoline quinone
    'C00120',  # Biotin (Vitamin B7)
    'C00143',  # Methylenetetrahydrofolate (Vitamin B9)
    'C00175',  # Cobalt ion (Co2+)
    'C00194',  # Cobamide coenzyme (Vitamin B12)
    'C00272',  # Tetrahydrobiopterin
    'C00305',  # Magnesium ion (Mg2+)
    'C00415',  # Dihydrofolic acid (Vitamin B9)
    'C00828',  # Menaquinone (Vitamin K2)
    'C00862',  # Methanofuran
    'C01217',  # Tetrahydromethanopterin
    'C01352',  # Flavin adenine dinucleotide (FADH2)
    'C02059',  # Phylloquinone (Vitamin K1)
    'C03576',  # Coenzyme M
    'C04628',  # Coenzyme B
    'C05924',  # Molybdopterin
    'C06453',  # Methylcobalamin (Vitamin B12)
    'C11378',  # Ubiquinone-10 (Coenzyme Q10)
    'C14818',  # Iron (Fe2+)
    'C14819',  # Iron (Fe3+)
    'C15670',  # Heme A
    'C15672',  # Heme O
    'C15817',  # Heme C
    'C16241',  # Lipoic acid
    'C18237',  # Molybdenum cofactor
    'C19153',  # Coenzyme F420-0
    'C19609',  # Nickel ion (Ni2+)
    'C19610',  # Manganese (Mn2+)
    'C22424',  # Copper ion (Cu2+)
]

# -------------------- Parse Reaction --------------------
if '<=>' not in raw_reaction:
    print(f"Reaction format error: missing '<=>'")
    sys.exit(1)

substrates_raw, products_raw = raw_reaction.split('<=>')
substrate_ids = [x.strip().split()[-1] for x in substrates_raw.strip().split('+')]
product_ids = [x.strip().split()[-1] for x in products_raw.strip().split('+')]

compound_smiles_map = {
    row['compound_id']: row['smiles'] for _, row in df_molecules.iterrows()
    if pd.notna(row['smiles'])
}


# --- Toggle for ignoring compounds present on both sides ---
def omit_shared_compounds(sub_ids, prod_ids):
    shared = set(sub_ids) & set(prod_ids)
    if shared:
        print(f"[INFO] Omitting compounds present on both sides: {shared}")
    return [cid for cid in sub_ids if cid not in shared], [cid for cid in prod_ids if cid not in shared]


if args.dynamic_cofactors:
    substrate_ids, product_ids = omit_shared_compounds(substrate_ids, product_ids)

# Filter out empty SMILES and compounds not found in mapping
substrate_smiles_list = []
for cid in substrate_ids:
    print(f"[DEBUG] Processing substrate compound: {cid}")
    if cid in compound_smiles_map:
        smiles = compound_smiles_map[cid]
        print(f"[DEBUG] Found SMILES for {cid}: '{smiles}'")
        if smiles and smiles.strip():  # Check if SMILES is not empty or just whitespace
            substrate_smiles_list.append(smiles.strip())
        else:
            print(f"[WARNING] Empty SMILES found for substrate compound {cid}")
    else:
        print(f"[WARNING] Substrate compound {cid} not found in reference file")

product_smiles_list = []
for cid in product_ids:
    print(f"[DEBUG] Processing product compound: {cid}")
    if cid in compound_smiles_map:
        smiles = compound_smiles_map[cid]
        print(f"[DEBUG] Found SMILES for {cid}: '{smiles}'")
        if smiles and smiles.strip():  # Check if SMILES is not empty or just whitespace
            product_smiles_list.append(smiles.strip())
        else:
            print(f"[WARNING] Empty SMILES found for product compound {cid}")
    else:
        print(f"[WARNING] Product compound {cid} not found in reference file")

substrate_smiles = '.'.join(substrate_smiles_list)
product_smiles = '.'.join(product_smiles_list)

# -------------------- Filter Cofactors --------------------
static_cofactors = args.static_cofactors

# --- Remove stoichiometry numbers from reaction string ---
import re


def remove_coefficients(reaction_string):
    """Remove numerical coefficients and variable coefficients from the reaction string"""
    # Remove numbers or variable coefficients followed by space before compound names
    # This handles: "2 compound", "n compound", "m compound", etc.
    reaction_string = re.sub(r'\b(\d+|n|m)\s+', '', reaction_string)

    # Remove variable coefficients that appear as separate terms
    # This handles: "n", "m", "m-1", "n-x", "x", etc. when they appear between + signs
    # But be careful not to remove parts of compound names
    reaction_string = re.sub(r'\+\s*([nm](-\d+)?|[nm]-[xn]|x)\s*\+', ' + ', reaction_string)
    reaction_string = re.sub(r'^([nm](-\d+)?|[nm]-[xn]|x)\s*\+', '', reaction_string)
    reaction_string = re.sub(r'\+\s*([nm](-\d+)?|[nm]-[xn]|x)\s*$', '', reaction_string)

    # Handle cases where coefficient terms are at the beginning/end without +
    reaction_string = re.sub(r'^([nm](-\d+)?|[nm]-[xn]|x)\s+', '', reaction_string)
    reaction_string = re.sub(r'\s+([nm](-\d+)?|[nm]-[xn]|x)$', '', reaction_string)

    # Special handling for complex coefficient expressions like "n-x+x"
    # Replace "n-x+x" with empty string (as it equals n)
    reaction_string = re.sub(r'\bn-x\+x\b', 'n', reaction_string)

    # Clean up any double spaces and multiple + signs
    reaction_string = re.sub(r'\s+', ' ', reaction_string).strip()
    reaction_string = re.sub(r'\+\s*\+', '+', reaction_string)
    reaction_string = re.sub(r'^\+\s*', '', reaction_string)
    reaction_string = re.sub(r'\s*\+$', '', reaction_string)

    return reaction_string


def replace_compounds_with_smiles(reaction_string, compound_mapping):
    """Replace compound IDs with their corresponding SMILES"""
    if not reaction_string or not reaction_string.strip():
        print(f"Warning: Empty or whitespace-only reaction string")
        return ""

    # Split by reaction arrow to handle reactants and products separately
    if '<=>' in reaction_string:
        reactants, products = reaction_string.split('<=>')
    elif '=>' in reaction_string:
        reactants, products = reaction_string.split('=>')
    elif '->' in reaction_string:
        reactants, products = reaction_string.split('->')
    else:
        # If no arrow found, treat entire string as reactants
        reactants = reaction_string
        products = ''

    def replace_in_side(side):
        if not side or not side.strip():
            return ""

        # Split by + to get individual compounds
        compounds = [comp.strip() for comp in side.split('+') if comp.strip()]
        smiles_list = []

        for compound in compounds:
            compound = compound.strip()  # Extra strip for safety

            # Skip empty strings and coefficient-only terms
            if not compound or compound in ['n', 'm', 'x', 'm-1', 'n-x'] or re.match(r'^[nm](-\d+)?$|^[nm]-[xn]$|^x$',
                                                                                     compound):
                print(f"Debug: Skipping coefficient/empty term: '{compound}'")
                continue

            # Remove any remaining coefficient patterns from individual compound names
            original_compound = compound
            compound = re.sub(r'^(\d+|[nm](-\d+)?|[nm]-[xn]|x)\s+', '', compound)

            if compound and compound in compound_mapping:
                smiles = compound_mapping[compound]
                if smiles and smiles.strip():  # Check if SMILES is not empty
                    smiles_list.append(smiles.strip())
                else:
                    print(f"Warning: Empty SMILES found for compound '{compound}'")
            elif compound:
                # If compound not found in mapping, keep original and print warning
                print(f"Warning: Compound '{compound}' not found in reference file")
                smiles_list.append(compound)
            else:
                print(f"Warning: After coefficient removal, compound became empty. Original: '{original_compound}'")

        # Join with . instead of +
        result = '.'.join(smiles_list)
        if not result:
            print(f"Warning: No valid compounds found in side: '{side}'")
        return result

    reactants_smiles = replace_in_side(reactants)
    if products:
        products_smiles = replace_in_side(products)
        return f"{reactants_smiles}>>{products_smiles}"
    else:
        return reactants_smiles


# Save reaction after removing stoichiometry
reaction_no_stoich = remove_coefficients(raw_reaction)


# --- Save reaction after cofactor filtering ---

def get_cofactor_filtered_reaction(substrate_ids, product_ids, cofactor_ids, remove_cofactors):
    # Remove cofactors from both sides if remove_cofactors is True
    if not remove_cofactors:
        return ''
    filtered_subs = [cid for cid in substrate_ids if cid not in cofactor_ids]
    filtered_prods = [cid for cid in product_ids if cid not in cofactor_ids]
    # If no change, return empty string
    if filtered_subs == substrate_ids and filtered_prods == product_ids:
        return 'No cofactor found'
    sub_names = [cid for cid in filtered_subs]
    prod_names = [cid for cid in filtered_prods]
    return ' + '.join(sub_names) + ' <=> ' + ' + '.join(prod_names)


# Show cofactor-filtered reaction if either static or dynamic cofactor flag is set
remove_cofactors_flag = args.static_cofactors or args.dynamic_cofactors
reaction_after_cofactor = get_cofactor_filtered_reaction(substrate_ids, product_ids, cofactor_ids,
                                                         remove_cofactors_flag)

if static_cofactors:
    cofactor_smiles = set()
    for cid in cofactor_ids:
        if cid in compound_smiles_map:
            canonical = canonicalize_smiles(compound_smiles_map[cid])
            if canonical:
                cofactor_smiles.add(canonical)


    def filter_cofactors(smiles_string):
        if not smiles_string:
            return ""
        return '.'.join([
            smi for smi in smiles_string.strip().split('.')
            if canonicalize_smiles(smi) not in cofactor_smiles
        ])


    substrate_smiles = filter_cofactors(substrate_smiles)
    product_smiles = filter_cofactors(product_smiles)

# Add debug information about the SMILES being processed
print(f"[DEBUG] Substrate SMILES: '{substrate_smiles}'")
print(f"[DEBUG] Product SMILES: '{product_smiles}'")

# Check if we have valid SMILES before processing
if not substrate_smiles or not substrate_smiles.strip():
    print(f"[ERROR] No valid substrate SMILES found for reaction {selected_id}. Halting.")
    sys.exit(0)

if not product_smiles or not product_smiles.strip():
    print(f"[ERROR] No valid product SMILES found for reaction {selected_id}. Halting.")
    sys.exit(0)

substrate_mols = smiles_to_mols(substrate_smiles)

if not substrate_mols:
    print(f"No non-cofactor substrates found for reaction {selected_id}. Halting.")
    sys.exit(0)


# -------------------- Auto-select Molecule --------------------

# --- Chiral center counting function ---
def count_chiral_centers(mol):
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    return len(chiral_centers)


# -------------------- Atom Tracking Function --------------------
def perform_atom_tracking(substrate_mols, product_smiles, mol_idx, atom_idx, selected_id, mapper_name=""):
    """
    Perform atom tracking for a given substrate molecule and atom.
    Returns tuple: (product_cid, atom_mapped_rxn, tracking_successful, mapper_used, mapping_error)
    """
    # Store original atom properties for comparison
    original_mol = substrate_mols[mol_idx]
    target_atom = original_mol.GetAtomWithIdx(atom_idx)
    target_symbol = target_atom.GetSymbol()
    target_degree = target_atom.GetDegree()

    # Re-label the atom in case it was modified
    substrate_mols[mol_idx] = label_atom(substrate_mols[mol_idx], atom_idx, 1)

    # Generate reaction SMILES
    reactant_smiles = Chem.MolToSmiles(substrate_mols[mol_idx], canonical=False)
    reaction_smiles = f"{reactant_smiles}>>{product_smiles}"

    print(f"[DEBUG] {mapper_name} Reaction SMILES for mapping: {reaction_smiles}")

    # Perform atom mapping
    mapped_rxn, mapper_used, mapping_error = perform_atom_mapping(reaction_smiles, selected_id)

    if mapped_rxn:
        atom_mapped_rxn = mapped_rxn
        try:
            reactant_part, product_part = mapped_rxn.split('>>')
            print(f"[DEBUG] {mapper_name} Reactant part: {reactant_part}")
            print(f"[DEBUG] {mapper_name} Product part: {product_part}")

            # Parse reactant molecules to find the map number assigned to our target atom
            reactant_mols = smiles_to_mols(reactant_part)
            target_map_num = None

            # Debug: Show all atom map numbers in reactant molecules first
            print(f"[DEBUG] {mapper_name} Number of reactant molecules parsed: {len(reactant_mols)}")
            reactant_map_nums = []
            for i, mol in enumerate(reactant_mols):
                mol_map_nums = [atom.GetAtomMapNum() for atom in mol.GetAtoms() if atom.GetAtomMapNum() > 0]
                print(f"[DEBUG] {mapper_name} Reactant molecule {i} has atom map numbers: {mol_map_nums}")
                reactant_map_nums.extend(mol_map_nums)
            print(f"[DEBUG] {mapper_name} All atom map numbers in reactants: {set(reactant_map_nums)}")

            # Find which map number was assigned to our target atom in the mapped reaction
            print(
                f"[DEBUG] {mapper_name} Looking for target atom: {target_symbol} with degree {target_degree} at position {atom_idx}")

            # Strategy 1: Look for our labeled atom (map number 1) in reactants first
            for r_idx, r_mol in enumerate(reactant_mols):
                for atom in r_mol.GetAtoms():
                    if atom.GetAtomMapNum() == 1:
                        target_map_num = 1
                        print(f"[DEBUG] {mapper_name} Found our labeled atom (map 1) preserved in reactants")
                        break
                if target_map_num is not None:
                    break

            # Strategy 2: If map 1 not found, find atoms matching our target properties
            if target_map_num is None:
                print(f"[DEBUG] {mapper_name} Map number 1 not preserved, searching by atom properties...")
                # Look for atoms with same symbol and degree in the expected molecule
                if mol_idx < len(reactant_mols):
                    r_mol = reactant_mols[mol_idx]
                    matching_atoms = []
                    for atom in r_mol.GetAtoms():
                        if (atom.GetSymbol() == target_symbol and
                                atom.GetDegree() == target_degree and
                                atom.GetAtomMapNum() > 0):
                            matching_atoms.append(atom.GetAtomMapNum())

                    if matching_atoms:
                        target_map_num = matching_atoms[0]  # Take the first match
                        print(
                            f"[DEBUG] {mapper_name} Found matching {target_symbol} atom with map number: {target_map_num}")

            # Strategy 3: Broader search if still not found
            if target_map_num is None:
                print(f"[DEBUG] {mapper_name} Broader search for any {target_symbol} atom...")
                for r_mol in reactant_mols:
                    for atom in r_mol.GetAtoms():
                        if atom.GetSymbol() == target_symbol and atom.GetAtomMapNum() > 0:
                            target_map_num = atom.GetAtomMapNum()
                            print(f"[DEBUG] {mapper_name} Found {target_symbol} atom with map number: {target_map_num}")
                            break
                    if target_map_num is not None:
                        break

            # Fallback to looking for map number 1
            if target_map_num is None:
                target_map_num = 1
                print(f"[DEBUG] {mapper_name} Could not determine target map number, using default: 1")

            product_mols = smiles_to_mols(product_part)
            print(f"[DEBUG] {mapper_name} Number of product molecules parsed: {len(product_mols)}")

            # Debug: Show all atom map numbers in product molecules
            all_map_nums = []
            for i, mol in enumerate(product_mols):
                mol_map_nums = [atom.GetAtomMapNum() for atom in mol.GetAtoms() if atom.GetAtomMapNum() > 0]
                print(f"[DEBUG] {mapper_name} Product molecule {i} has atom map numbers: {mol_map_nums}")
                all_map_nums.extend(mol_map_nums)

            print(f"[DEBUG] {mapper_name} All atom map numbers in products: {set(all_map_nums)}")

            # Check mapping consistency
            reactant_set = set(reactant_map_nums)
            product_set = set(all_map_nums)
            common_maps = reactant_set & product_set
            only_reactants = reactant_set - product_set
            only_products = product_set - reactant_set

            print(f"[DEBUG] {mapper_name} Map numbers in both reactants and products: {len(common_maps)} atoms")
            if only_reactants:
                print(f"[DEBUG] {mapper_name} Map numbers only in reactants: {only_reactants}")
            if only_products:
                print(f"[DEBUG] {mapper_name} Map numbers only in products: {only_products}")

            print(f"[DEBUG] {mapper_name} Looking for atom map number: {target_map_num}")

            matched = trace_atom(target_map_num, product_mols)

            if matched:
                clean_smiles = remove_atom_map_numbers(matched)
                product_cid = lookup_compound_id(clean_smiles)
                print(f"[DEBUG] {mapper_name} Matched Product SMILES: {clean_smiles} → {product_cid}")
                return product_cid, atom_mapped_rxn, True, mapper_used, None
            else:
                print(f"[DEBUG] {mapper_name} Atom with map number {target_map_num} not found in any product.")
                # Final fallback: look for any atom of the same element in products
                print(
                    f"[DEBUG] {mapper_name} Final fallback: searching for any mapped {target_symbol} atoms in products...")
                for p_idx, p_mol in enumerate(product_mols):
                    for atom in p_mol.GetAtoms():
                        if atom.GetSymbol() == target_symbol and atom.GetAtomMapNum() > 0:
                            print(
                                f"[DEBUG] {mapper_name} Found mapped {target_symbol} atom in product molecule {p_idx} with map number: {atom.GetAtomMapNum()}")
                            clean_smiles = remove_atom_map_numbers(p_mol)
                            product_cid = lookup_compound_id(clean_smiles)
                            print(
                                f"[DEBUG] {mapper_name} Fallback {target_symbol} match - Product SMILES: {clean_smiles} → {product_cid}")
                            return product_cid, atom_mapped_rxn, True, mapper_used, None
                return 'NotFound', atom_mapped_rxn, False, mapper_used, None

        except Exception as e:
            post_error = f"Post-mapping processing error: {str(e)}"
            print(f"[ERROR] {mapper_name} Error processing mapped reaction: {e}")
            return 'Error', 'Error', False, mapper_used, post_error
    else:
        print(f"[ERROR] {mapper_name} Mapping failed: {mapping_error}")
        return 'Error', 'Error', False, 'failed', mapping_error


# --- Molecule selection logic ---

def selection_key_highest(i):
    ccc = count_chiral_centers(substrate_mols[i])
    heavy = substrate_mols[i].GetNumHeavyAtoms()
    # For tie-breaking: sort by chiral centers, then by heavy atom count (descending)
    return (ccc, heavy)


def selection_key_lowest(i):
    ccc = count_chiral_centers(substrate_mols[i])
    heavy = substrate_mols[i].GetNumHeavyAtoms()
    # For lowest: sort by chiral centers ascending, then by heavy atom count descending
    return (ccc, -heavy)


if args.highest_CCC:
    # Highest chiral center count, break ties with highest heavy atom count
    mol_idx = max(range(len(substrate_mols)), key=selection_key_highest)
    reason = 'highest chiral center count'
elif args.lowest_CCC:
    # Lowest chiral center count, break ties with highest heavy atom count
    mol_idx = min(range(len(substrate_mols)), key=selection_key_lowest)
    reason = 'lowest chiral center count'
else:
    mol_idx = max(range(len(substrate_mols)), key=lambda i: substrate_mols[i].GetNumHeavyAtoms())
    reason = 'highest heavy atom count (default)'

mol = substrate_mols[mol_idx]
atom_idx = max(range(mol.GetNumAtoms()), key=lambda i: mol.GetAtomWithIdx(i).GetDegree())
substrate_mols[mol_idx] = label_atom(substrate_mols[mol_idx], atom_idx, 1)

print(f"[DEBUG] Selected molecule index: {mol_idx}, atom index: {atom_idx}, selection reason: {reason}")

# Debug: Show information about the labeled atom
labeled_mol = substrate_mols[mol_idx]
labeled_atom = labeled_mol.GetAtomWithIdx(atom_idx)
print(f"[DEBUG] Labeled atom details:")
print(f"[DEBUG]   - Atom symbol: {labeled_atom.GetSymbol()}")
print(f"[DEBUG]   - Atom degree: {labeled_atom.GetDegree()}")
print(f"[DEBUG]   - Atom map number: {labeled_atom.GetAtomMapNum()}")
print(f"[DEBUG]   - Substrate molecule SMILES with label: {Chem.MolToSmiles(labeled_mol, canonical=False)}")

# -------------------- Perform Mapping --------------------

# First attempt with initial mapping
print(f"[INFO] Starting atom tracking for reaction {selected_id}")
product_cid, atom_mapped_rxn, tracking_successful, mapper_used, mapping_error = perform_atom_tracking(
    substrate_mols, product_smiles, mol_idx, atom_idx, selected_id, "[INITIAL]"
)

# If tracking failed and we used localmapper fallback, try alternative atom selection
if not tracking_successful and mapper_used == "localmapper_fallback":
    print(f"[INFO] Initial tracking with localmapper failed, trying alternative atom selection...")

    # Try with a different atom selection strategy
    # First, try the atom with second highest degree
    degrees = [(i, substrate_mols[mol_idx].GetAtomWithIdx(i).GetDegree()) for i in
               range(substrate_mols[mol_idx].GetNumAtoms())]
    degrees.sort(key=lambda x: x[1], reverse=True)

    if len(degrees) > 1:
        alt_atom_idx = degrees[1][0]  # Second highest degree atom
        print(f"[INFO] Trying alternative atom: index {alt_atom_idx} (degree {degrees[1][1]})")

        alt_product_cid, alt_atom_mapped_rxn, alt_tracking_successful, alt_mapper_used, alt_mapping_error = perform_atom_tracking(
            substrate_mols, product_smiles, mol_idx, alt_atom_idx, selected_id, "[ALTERNATIVE]"
        )

        # If alternative tracking was successful, use those results
        if alt_tracking_successful:
            product_cid = alt_product_cid
            atom_mapped_rxn = alt_atom_mapped_rxn
            tracking_successful = True
            atom_idx = alt_atom_idx  # Update the atom index for logging
            mapper_used = alt_mapper_used
            mapping_error = alt_mapping_error
            print(f"[INFO] Alternative atom tracking successful!")
        else:
            print(f"[INFO] Alternative atom tracking also failed, keeping original results")

# Handle visualization if requested
if atom_mapped_rxn != 'Error' and args.visualize:
    try:
        visualize_mapped_rxn(atom_mapped_rxn, selected_id, save_images=True)
    except Exception as e:
        print(f"[WARNING] Visualization failed: {e}")

print(f"[INFO] Final result - Mapping: {mapper_used}, Tracking: {'successful' if tracking_successful else 'failed'}")

# Set final mapping_error if needed
if not tracking_successful and mapping_error is None:
    mapping_error = "Atom tracking failed - atom not found in products"


def save_error_reaction(error_output_file, error_lock_file, reaction_data):
    """
    Save error reaction data to the error CSV file with the same structure as the main output.
    """
    error_row = pd.DataFrame([reaction_data])

    # Use file lock to ensure safe concurrent writes for error file
    with FileLock(error_lock_file, timeout=60):
        if os.path.exists(error_output_file) and os.path.getsize(error_output_file) > 0:
            df_existing_errors = pd.read_csv(error_output_file, quoting=1)  # QUOTE_ALL for better CSV handling
            df_combined_errors = pd.concat([df_existing_errors, error_row], ignore_index=True)
            df_combined_errors = df_combined_errors.sort_values(by='reaction_id').reset_index(drop=True)
            df_combined_errors.to_csv(error_output_file, index=False, quoting=1)  # QUOTE_ALL
        else:
            error_row.to_csv(error_output_file, index=False, quoting=1)  # QUOTE_ALL

    print(f"[INFO] Saved error reaction to {error_output_file}")


def process_error_reactions_smiles(error_output_file, compound_mapping):
    """
    Process the Error Reactions CSV file and add reaction_smiles column.
    """
    if not os.path.exists(error_output_file):
        print(f"[INFO] No error reactions file found at {error_output_file}")
        return

    print(f"[INFO] Processing error reactions to add SMILES...")

    # Read the error reactions file with proper CSV handling for long strings
    try:
        df_errors = pd.read_csv(error_output_file, quoting=1)  # QUOTE_ALL for better CSV handling
    except Exception as e:
        print(f"[ERROR] Failed to read error reactions file: {e}")
        print(f"[INFO] Attempting to read with different CSV parameters...")
        try:
            # Try with different parameters if the file is malformed
            df_errors = pd.read_csv(error_output_file, quoting=1, error_bad_lines=False, warn_bad_lines=True)
        except Exception as e2:
            print(f"[ERROR] Still failed to read error reactions file: {e2}")
            return

    # Check if reaction_smiles column already exists and inform user it will be overwritten
    if 'reaction_smiles' in df_errors.columns:
        print(f"[INFO] reaction_smiles column already exists in error reactions file - overwriting...")

    # Process each reaction to generate SMILES
    reaction_smiles_list = []
    for index, row in df_errors.iterrows():
        reaction_string = None

        # Try to get reaction string from various possible columns
        if 'reaction_no_stoich' in row and pd.notna(row['reaction_no_stoich']):
            reaction_string = str(row['reaction_no_stoich'])
        elif 'reaction_id' in row and pd.notna(row['reaction_id']):
            # Try to look up the original reaction from the input file
            reaction_id = row['reaction_id']
            matching_reactions = df_input_reactions[df_input_reactions['id'] == reaction_id]
            if not matching_reactions.empty:
                original_reaction = matching_reactions['reaction'].iloc[0]
                reaction_string = remove_coefficients(original_reaction)
                print(f"[INFO] Found original reaction for {reaction_id}: {original_reaction}")

        if reaction_string:
            # Replace compounds with SMILES
            reaction_smiles = replace_compounds_with_smiles(reaction_string, compound_mapping)
            reaction_smiles_list.append(reaction_smiles)
        else:
            print(f"[WARNING] No reaction string found for error row {index}")
            reaction_smiles_list.append('')

    # Add the reaction_smiles column
    df_errors['reaction_smiles'] = reaction_smiles_list

    # Save the updated file with proper CSV quoting
    df_errors.to_csv(error_output_file, index=False, quoting=1)  # QUOTE_ALL
    print(f"[INFO] Added reaction_smiles column to {error_output_file}")
    print(f"[INFO] Processed {len(reaction_smiles_list)} error reactions")


# -------------------- Get Substrate ID --------------------
selected_mol = substrate_mols[mol_idx]
selected_clean = remove_atom_map_numbers(selected_mol)
substrate_cid = lookup_compound_id(selected_clean)
print(f"[DEBUG] Selected Substrate SMILES: {selected_clean} → {substrate_cid}")

# -------------------- Save Output --------------------

# Prepare the output data
# Ensure mapper_used is always set
if 'mapper_used' not in locals() or mapper_used is None:
    mapper_used = 'unknown'

output_data = {
    'reaction_id': selected_id,
    'tracking': f"{substrate_cid} => {product_cid}",
    'atom_mapped_rxn': atom_mapped_rxn,
    'reaction_no_stoich': reaction_no_stoich,
    'reaction_after_cofactor': reaction_after_cofactor if reaction_after_cofactor is not None else '',
    'selection_method': reason,
    'mapper_used': mapper_used,
    'cofactor_handling': (
        'static' if args.static_cofactors and not args.dynamic_cofactors else
        'dynamic' if args.dynamic_cofactors and not args.static_cofactors else
        'both' if args.static_cofactors and args.dynamic_cofactors else
        'none'
    )
}

# If there was a mapping error, add error information and save to error file
if mapping_error:
    output_data['error_message'] = mapping_error
    save_error_reaction(error_output_file, error_lock_file, output_data)
else:
    # Save successful reaction to main output file
    output_row = pd.DataFrame([output_data])

    # --- Use file lock to ensure safe concurrent writes ---
    with FileLock(lock_file, timeout=60):
        if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
            df_existing = pd.read_csv(output_file)
            df_combined = pd.concat([df_existing, output_row], ignore_index=True)
            df_combined = df_combined.sort_values(by='reaction_id').reset_index(drop=True)
            df_combined.to_csv(output_file, index=False)
        else:
            output_row.to_csv(output_file, index=False)

    print(f"[INFO] Saved tracking result to {output_file}")

# -------------------- Process Error Reactions for SMILES --------------------
# Process error reactions to add SMILES column
if os.path.exists(error_output_file) and os.path.getsize(error_output_file) > 0:
    try:
        process_error_reactions_smiles(error_output_file, compound_smiles_map)
    except Exception as e:
        print(f"[WARNING] Failed to process error reactions for SMILES: {e}")
else:
    print(f"[INFO] No error reactions file found to process for SMILES")
