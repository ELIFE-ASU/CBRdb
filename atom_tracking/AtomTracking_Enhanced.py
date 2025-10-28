"""============================================================================
                        CHEMICAL REACTION ATOM TRACKING SYSTEM
===============================================================================
DESCRIPTION:
Automated atom tracking in chemical reactions using advanced atom mapping algorithms.
Identifies specific atoms in substrate molecules and traces their fate through chemical transformations.
PURPOSE:
- Track the fate of specific atoms through chemical reactions
- Handle complex reaction equations with coefficients and cofactors
- Provide robust atom mapping with multiple fallback strategies
- Support parallel processing for large reaction datasets
- Generate visualizations of mapped reactions
WORKFLOW OVERVIEW:
1. Reaction parsing and preprocessing
2. Molecule and atom selection
3. Atom mapping (RXNMapper/localmapper fallback)
4. Atom tracking through products
5. Output results and error handling
TECHNICAL COMPONENTS:
- RDKit: Molecular processing and SMILES handling
- RXNMapper: Primary atom mapping (512 token limit)
- localmapper: Fallback mapping for large reactions
- pandas: Data processing and CSV operations
- FileLock: Thread-safe file operations
============================================================================"""
import argparse
import os
import re
import sys

import matplotlib.pyplot as plt
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rxnmapper import RXNMapper

# --- File locking for safe parallel writes ---
try:
    from filelock import FileLock
except ImportError:
    import subprocess

    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'filelock'])
    from filelock import FileLock
# --- Import localmapper for fallback mapping ---
try:
    from localmapper import localmapper

    LOCALMAPPER_AVAILABLE = True
except ImportError:
    LOCALMAPPER_AVAILABLE = False


# -------------------- Helper Functions --------------------
def canonicalize_smiles(smiles):
    """Convert a SMILES string to its canonical form, removing all atom map numbers.
    Returns None if input is invalid or cannot be parsed."""
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def smiles_to_mols(smiles):
    """Convert a dot-separated SMILES string into a list of RDKit Mol objects.
    Only valid fragments are returned as Mol objects."""
    if not smiles or not smiles.strip():
        return []
    fragments = [smi.strip() for smi in smiles.strip().split('.') if smi.strip()]
    mols = []
    for smi in fragments:
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                mols.append(mol)
        except Exception:
            pass
    return mols


def label_atom(mol, atom_idx, map_num):
    """Label a specific atom in a molecule with a given atom map number.
    Returns a new RDKit Mol object with the label applied."""
    mol = Chem.RWMol(mol)
    if 0 <= atom_idx < mol.GetNumAtoms():
        mol.GetAtomWithIdx(atom_idx).SetAtomMapNum(map_num)
    return mol


def trace_atom(atom_map_num, mols):
    """Search a list of molecules for an atom with a specific atom map number.
    Returns the first molecule containing the atom, or None if not found."""
    for mol in mols:
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() == atom_map_num:
                return mol
    return None


def remove_atom_map_numbers(mol):
    """Remove all atom map numbers from a molecule and return its canonical SMILES.
    Returns None if input is None."""
    if mol is None:
        return None
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def visualize_mapped_rxn(mapped_rxns, reaction_id, save_images=False):
    """Visualize the mapped reaction SMILES. If save_images is True, saves the image
    in a folder named 'Reactions Visualizations' using the reaction ID as the filename."""
    try:
        rxn = Chem.rdChemReactions.ReactionFromSmarts(mapped_rxns, useSmiles=True)
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
    except Exception:
        pass


def find_file(filename):
    """Search for a file in the current script directory and up to two parent directories.
    Returns the absolute path if found, otherwise raises FileNotFoundError."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    for base in [script_dir, os.path.join(script_dir, ".."), os.path.join(script_dir, "..", "..")]:
        candidate = os.path.abspath(os.path.join(base, filename))
        if os.path.exists(candidate):
            return candidate
    raise FileNotFoundError(f"{filename} not found in current or parent directories.")


def get_or_build_canonical_index(df_molecules, output_file):
    """Load a canonical index from CSV if it exists, otherwise build it from the molecule DataFrame.
    Each fragment of each molecule is canonicalized and indexed."""
    if os.path.exists(output_file):
        return pd.read_csv(output_file)
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
    df_index = pd.DataFrame(index)
    df_index.to_csv(output_file, index=False)
    return df_index


def remove_coefficients(reaction_string):
    """Remove numeric or symbolic coefficients (e.g. 2, n, m+1, n-1+x) from reactions."""
    # Remove coefficient tokens before compounds
    reaction_string = re.sub(r'\b[0-9a-zA-Z+\-*/]+\s+(?=C\d{5})', '', reaction_string)

    # Collapse multiple spaces and clean operators
    reaction_string = re.sub(r'\s+', ' ', reaction_string)
    reaction_string = re.sub(r'\+\s*\+', '+', reaction_string)
    reaction_string = re.sub(r'^\+\s*', '', reaction_string)
    reaction_string = re.sub(r'\s*\+$', '', reaction_string)

    return reaction_string.strip()


def replace_compounds_with_smiles(reaction_string, compound_mapping):
    """Replace compound IDs with their corresponding SMILES in a reaction string."""
    import re
    if not reaction_string or not reaction_string.strip():
        return ""
    if '<=>' in reaction_string:
        reactants, products = reaction_string.split('<=>')
    elif '=>' in reaction_string:
        reactants, products = reaction_string.split('=>')
    elif '->' in reaction_string:
        reactants, products = reaction_string.split('->')
    else:
        reactants = reaction_string
        products = ''

    def replace_in_side(side):
        if not side or not side.strip():
            return ""
        compounds = [comp.strip() for comp in side.split('+') if comp.strip()]
        smiles_list = []
        for compound in compounds:
            if not compound or compound in ['n', 'm', 'x', 'm-1', 'n-x'] or re.match(r'^[nm](-\d+)?$|^[nm]-[xn]$|^x$',
                                                                                     compound):
                continue
            compound = re.sub(r'^(\d+|[nm](-\d+)?|[nm]-[xn]|x)\s+', '', compound)
            if compound and compound in compound_mapping:
                smiles = compound_mapping[compound]
                if smiles and smiles.strip():
                    smiles_list.append(smiles.strip())
        return '.'.join(smiles_list)

    reactants_smiles = replace_in_side(reactants)
    if products:
        products_smiles = replace_in_side(products)
        return f"{reactants_smiles}>>{products_smiles}"
    else:
        return reactants_smiles


def save_error_reaction(error_output_file, error_lock_file, reaction_data):
    """Save error reaction data to the error CSV file with the same structure as the main output."""
    error_row = pd.DataFrame([reaction_data])
    with FileLock(error_lock_file, timeout=60):
        if os.path.exists(error_output_file) and os.path.getsize(error_output_file) > 0:
            df_existing_errors = pd.read_csv(error_output_file, quoting=1)
            df_combined_errors = pd.concat([df_existing_errors, error_row], ignore_index=True)
            df_combined_errors = df_combined_errors.sort_values(by='reaction_id').reset_index(drop=True)
            df_combined_errors.to_csv(error_output_file, index=False, quoting=1)
        else:
            error_row.to_csv(error_output_file, index=False, quoting=1)


def log_and_exit_early(error_message, selected_id=None, reaction_no_stoich=None, reason=None):
    """Log error and exit script early."""
    output_data = {
        'reaction_id': selected_id if selected_id is not None else '',
        'tracking': '',
        'mapped_rxns': '',
        'reaction_no_stoich': reaction_no_stoich if reaction_no_stoich is not None else '',
        'reaction_after_cofactor': '',
        'selection_method': reason if reason is not None else '',
        'mapper_used': '',
        'cofactor_handling': 'unknown',
        'mapping_error': error_message,
        'error_message': error_message
    }
    try:
        save_error_reaction(error_output_file, error_lock_file, output_data)
    except Exception:
        pass
    sys.exit(1)


def process_error_reactions_smiles(error_output_file, compound_mapping):
    """Process the Error Reactions CSV file and add reaction_smiles column."""
    if not os.path.exists(error_output_file):
        return
    try:
        df_errors = pd.read_csv(error_output_file, quoting=1)
    except Exception:
        try:
            df_errors = pd.read_csv(error_output_file, quoting=1, error_bad_lines=False, warn_bad_lines=True)
        except Exception:
            return
    reaction_smiles_list = []
    for _, row in df_errors.iterrows():
        reaction_string = None
        if 'reaction_no_stoich' in row and pd.notna(row['reaction_no_stoich']):
            reaction_string = str(row['reaction_no_stoich'])
        elif 'reaction_id' in row and pd.notna(row['reaction_id']):
            reaction_id = row['reaction_id']
            matching_reactions = df_input_reactions[df_input_reactions['id'] == reaction_id]
            if not matching_reactions.empty:
                original_reaction = matching_reactions['reaction'].iloc[0]
                reaction_string = remove_coefficients(original_reaction)
        if reaction_string:
            reaction_smiles = replace_compounds_with_smiles(reaction_string, compound_mapping)
            reaction_smiles_list.append(reaction_smiles)
        else:
            reaction_smiles_list.append('')
    df_errors['reaction_smiles'] = reaction_smiles_list
    df_errors.to_csv(error_output_file, index=False, quoting=1)


def count_chiral_centers(mol):
    """Count the number of chiral centers in a molecule."""
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    return len(chiral_centers)


def perform_atom_mapping(reaction_smiles, reaction_id):
    """Perform atom mapping with RXNMapper first, fallback to localmapper if token limit exceeded.
    Returns tuple: (mapped_rxn, mapper_used, error_message)"""
    try:
        mapped_rxn = rxn_mapper.get_attention_guided_atom_maps([reaction_smiles])[0]['mapped_rxn']
        return mapped_rxn, "RXNMapper", None
    except Exception as e:
        error_msg = str(e)
        token_limit_keywords = ["tokens", "token", "512", "should be at most", "token limit", "maximum length"]
        is_token_limit_error = any(keyword in error_msg.lower() for keyword in token_limit_keywords)
        if is_token_limit_error:
            if LOCALMAPPER_AVAILABLE and local_mapper is not None:
                try:
                    mapped_rxn = local_mapper.get_atom_map(reaction_smiles)
                    return mapped_rxn, "localmapper_fallback", None
                except Exception as local_e:
                    local_error_msg = str(local_e)
                    return None, "both_failed", f"RXNMapper: {error_msg}; localmapper: {local_error_msg}"
            else:
                return None, "rxnmapper_failed_no_fallback", f"RXNMapper: {error_msg}; localmapper not available"
        else:
            return None, "rxnmapper_failed_other", error_msg


def perform_atom_tracking(substrate_mols, product_smiles, mol_idx, atom_idx, selected_id, mapper_name=""):
    """Perform atom tracking for a given substrate molecule and atom.
    Returns tuple: (product_cid, atom_mapped_rxn, tracking_successful, mapper_used, mapping_error)"""
    original_mol = substrate_mols[mol_idx]
    target_atom = original_mol.GetAtomWithIdx(atom_idx)
    target_symbol = target_atom.GetSymbol()
    target_degree = target_atom.GetDegree()
    substrate_mols[mol_idx] = label_atom(substrate_mols[mol_idx], atom_idx, 1)
    reactant_smiles = Chem.MolToSmiles(substrate_mols[mol_idx], canonical=False)
    reaction_smiles = f"{reactant_smiles}>>{product_smiles}"
    mapped_rxn, mapper_used, mapping_error = perform_atom_mapping(reaction_smiles, selected_id)
    if mapped_rxn:
        atom_mapped_rxn = mapped_rxn
        try:
            reactant_part, product_part = mapped_rxn.split('>>')
            reactant_mols = smiles_to_mols(reactant_part)
            target_map_num = None
            reactant_map_nums = []
            for mol in reactant_mols:
                mol_map_nums = [atom.GetAtomMapNum() for atom in mol.GetAtoms() if atom.GetAtomMapNum() > 0]
                reactant_map_nums.extend(mol_map_nums)
            for r_mol in reactant_mols:
                for atom in r_mol.GetAtoms():
                    if atom.GetAtomMapNum() == 1:
                        target_map_num = 1
                        break
                if target_map_num is not None:
                    break
            if target_map_num is None:
                if mol_idx < len(reactant_mols):
                    r_mol = reactant_mols[mol_idx]
                    matching_atoms = []
                    for atom in r_mol.GetAtoms():
                        if (atom.GetSymbol() == target_symbol and
                                atom.GetDegree() == target_degree and
                                atom.GetAtomMapNum() > 0):
                            matching_atoms.append(atom.GetAtomMapNum())
                    if matching_atoms:
                        target_map_num = matching_atoms[0]
            if target_map_num is None:
                for r_mol in reactant_mols:
                    for atom in r_mol.GetAtoms():
                        if atom.GetSymbol() == target_symbol and atom.GetAtomMapNum() > 0:
                            target_map_num = atom.GetAtomMapNum()
                            break
                    if target_map_num is not None:
                        break
            if target_map_num is None:
                target_map_num = 1
            product_mols = smiles_to_mols(product_part)
            all_map_nums = []
            for mol in product_mols:
                mol_map_nums = [atom.GetAtomMapNum() for atom in mol.GetAtoms() if atom.GetAtomMapNum() > 0]
                all_map_nums.extend(mol_map_nums)
            matched = trace_atom(target_map_num, product_mols)
            if matched:
                clean_smiles = remove_atom_map_numbers(matched)
                product_cid = lookup_compound_id(clean_smiles)
                return product_cid, atom_mapped_rxn, True, mapper_used, None
            else:
                for p_mol in product_mols:
                    for atom in p_mol.GetAtoms():
                        if atom.GetSymbol() == target_symbol and atom.GetAtomMapNum() > 0:
                            clean_smiles = remove_atom_map_numbers(p_mol)
                            product_cid = lookup_compound_id(clean_smiles)
                            return product_cid, atom_mapped_rxn, True, mapper_used, None
                return 'NotFound', atom_mapped_rxn, False, mapper_used, None
        except Exception as e:
            post_error = f"Post-mapping processing error: {str(e)}"
            return 'Error', 'Error', False, mapper_used, post_error
    else:
        return 'Error', 'Error', False, 'failed', mapping_error


def selection_key_highest(i):
    """Selection key for substrate with highest chiral center count, then heavy atom count."""
    ccc = count_chiral_centers(substrate_mols[i])
    heavy = substrate_mols[i].GetNumHeavyAtoms()
    return (ccc, heavy)


def selection_key_lowest(i):
    """Selection key for substrate with lowest chiral center count, then heavy atom count."""
    ccc = count_chiral_centers(substrate_mols[i])
    heavy = substrate_mols[i].GetNumHeavyAtoms()
    return (ccc, -heavy)


def omit_shared_compounds(sub_ids, prod_ids):
    """Remove compounds present on both sides of the reaction."""
    shared = set(sub_ids) & set(prod_ids)
    return [cid for cid in sub_ids if cid not in shared], [cid for cid in prod_ids if cid not in shared]


def get_cofactor_filtered_reaction(substrate_ids, product_ids, cofactor_ids, remove_cofactors):
    """Return reaction string after removing cofactors if requested."""
    if not remove_cofactors:
        return ''
    filtered_subs = [cid for cid in substrate_ids if cid not in cofactor_ids]
    filtered_prods = [cid for cid in product_ids if cid not in cofactor_ids]
    if filtered_subs == substrate_ids and filtered_prods == product_ids:
        return 'No cofactor found'
    sub_names = [cid for cid in filtered_subs]
    prod_names = [cid for cid in filtered_prods]
    return ' + '.join(sub_names) + ' <=> ' + ' + '.join(prod_names)


# -------------------- Main Script --------------------
parser = argparse.ArgumentParser()
parser.add_argument('--row_index', type=int, required=True)
parser.add_argument('--static_cofactors', action='store_true', help='Ignores Cofactors based on a static list')
parser.add_argument('--dynamic_cofactors', action='store_true', help='Ignores Compounds present on both sides.')
parser.add_argument('--highest_CCC', action='store_true', help='Select substrate with highest number of chiral centers')
parser.add_argument('--lowest_CCC', action='store_true', help='Select substrate with lowest number of chiral centers')
parser.add_argument('--visualize', action='store_true', help='If set, visualize the mapped reaction for this row index')
args = parser.parse_args()
reactions_file = 'Remainings.csv'
molecule_reference_file = 'CBRdb_C.csv'
script_dir = os.path.dirname(os.path.abspath(__file__))
canonical_index_file = os.path.join(script_dir, 'CBRdb_CanonicalIndex.csv')
output_file = os.path.join(script_dir, 'CBRdb_AtomTracking.csv')
error_output_file = os.path.join(script_dir, 'Error Reactions.csv')
lock_file = output_file + '.lock'
error_lock_file = error_output_file + '.lock'
try:
    reactions_file_path = find_file(reactions_file)
    molecule_ref_path = find_file(molecule_reference_file)
except FileNotFoundError as e:
    log_and_exit_early(str(e))
try:
    df_input_reactions = pd.read_csv(reactions_file_path, low_memory=False)
    df_molecules = pd.read_csv(molecule_ref_path)
except Exception as e:
    log_and_exit_early(f"Error loading input files: {e}")
if args.row_index >= len(df_input_reactions):
    log_and_exit_early(f"Row index {args.row_index} out of bounds.")
df_canonical_index = get_or_build_canonical_index(df_molecules, canonical_index_file)


def lookup_compound_id(canonical_smiles):
    """Look up the compound ID for a given canonical SMILES string in the canonical index DataFrame."""
    match = df_canonical_index[df_canonical_index['canonical_smiles'] == canonical_smiles]
    if not match.empty:
        return match['compound_id'].values[0]
    return 'Unknown'


selected_row = df_input_reactions.iloc[args.row_index]
selected_id = selected_row['id']
raw_reaction = selected_row['reaction']
rxn_mapper = RXNMapper()
local_mapper = localmapper() if LOCALMAPPER_AVAILABLE else None
if '<=>' not in raw_reaction:
    log_and_exit_early("Reaction format error: missing '<=>'", selected_id=selected_id, reaction_no_stoich=raw_reaction)
substrates_raw, products_raw = raw_reaction.split('<=>')
substrate_ids = [x.strip().split()[-1] for x in substrates_raw.strip().split('+')]
product_ids = [x.strip().split()[-1] for x in products_raw.strip().split('+')]
compound_smiles_map = {
    row['compound_id']: row['smiles'] for _, row in df_molecules.iterrows()
    if pd.notna(row['smiles'])
}
orig_substrate_ids = substrate_ids.copy()
orig_product_ids = product_ids.copy()
dynamic_skipped = False
if args.dynamic_cofactors:
    substrate_ids, product_ids = omit_shared_compounds(substrate_ids, product_ids)
    if not substrate_ids or not product_ids:
        substrate_ids = orig_substrate_ids
        product_ids = orig_product_ids
        dynamic_skipped = True
substrate_smiles_list = [compound_smiles_map[cid].strip() for cid in substrate_ids if
                         cid in compound_smiles_map and compound_smiles_map[cid].strip()]
product_smiles_list = [compound_smiles_map[cid].strip() for cid in product_ids if
                       cid in compound_smiles_map and compound_smiles_map[cid].strip()]
substrate_smiles = '.'.join(substrate_smiles_list)
product_smiles = '.'.join(product_smiles_list)
mapped_rxns = None
mapper_used = None
mapping_error = None
if substrate_smiles and product_smiles:
    raw_reaction_smiles = f"{substrate_smiles}>>{product_smiles}"
    mapped_rxns, mapper_used, mapping_error = perform_atom_mapping(raw_reaction_smiles, selected_id)
static_cofactors = args.static_cofactors
cofactor_ids = [
    'C00002', 'C00003', 'C00004', 'C00005', 'C00006', 'C00008', 'C00010', 'C00016', 'C00018', 'C00019', 'C00020',
    'C00032',
    'C00034', 'C00038', 'C00051', 'C00053', 'C00061', 'C00063', 'C00068', 'C00072', 'C00076', 'C00101', 'C00113',
    'C00120',
    'C00143', 'C00175', 'C00194', 'C00272', 'C00305', 'C00415', 'C00828', 'C00862', 'C01217', 'C01352', 'C02059',
    'C03576',
    'C04628', 'C05924', 'C06453', 'C11378', 'C14818', 'C14819', 'C15670', 'C15672', 'C15817', 'C16241', 'C18237',
    'C19153',
    'C19609', 'C19610', 'C22424'
]
reaction_no_stoich = remove_coefficients(raw_reaction)
remove_cofactors_flag = args.static_cofactors or args.dynamic_cofactors
# For reaction_after_cofactor, use shared compounds if dynamic (and not skipped)
if args.dynamic_cofactors and not dynamic_skipped:
    shared = set(orig_substrate_ids) & set(orig_product_ids)
    cofactor_ids_for_filter = list(shared)
else:
    cofactor_ids_for_filter = cofactor_ids if args.static_cofactors else []
reaction_after_cofactor = get_cofactor_filtered_reaction(substrate_ids, product_ids, cofactor_ids_for_filter,
                                                         remove_cofactors_flag)
if dynamic_skipped:
    reaction_after_cofactor = 'Dynamic filtering skipped - would empty side(s)'
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
if not substrate_smiles or not substrate_smiles.strip():
    log_and_exit_early(f"No valid substrate SMILES found for reaction {selected_id}.", selected_id=selected_id,
                       reaction_no_stoich=raw_reaction)
if not product_smiles or not product_smiles.strip():
    log_and_exit_early(f"No valid product SMILES found for reaction {selected_id}.", selected_id=selected_id,
                       reaction_no_stoich=raw_reaction)
substrate_mols = smiles_to_mols(substrate_smiles)
is_variable_coeff_reaction = any(var in raw_reaction.lower() for var in ['n', 'm', 'x'])
if not substrate_mols:
    if is_variable_coeff_reaction:
        shared = set(substrate_ids) & set(product_ids)
        if shared:
            main_cid = list(shared)[0]
            substrate_cid = main_cid
            product_cid = main_cid
            mapped_rxns = ''
            mapper_used = 'special'
            mapping_error = 'Handled as self-reaction due to variable coefficients and invalid SMILES'
            reason = 'special_variable_coeff'
            cofactor_handling = (
                'static' if args.static_cofactors and not args.dynamic_cofactors else
                'dynamic' if args.dynamic_cofactors and not args.static_cofactors else
                'both' if args.static_cofactors and args.dynamic_cofactors else
                'none'
            )
            if dynamic_skipped:
                cofactor_handling += '_skipped'
            print(f"reaction_id: {selected_id}")
            print(f"tracking: {substrate_cid} => {product_cid}")
            print(f"mapped_rxns: {mapped_rxns if mapped_rxns is not None else ''}")
            print(f"reaction_no_stoich: {reaction_no_stoich}")
            print(f"reaction_after_cofactor: {reaction_after_cofactor if reaction_after_cofactor is not None else ''}")
            print(f"selection_method: {reason}")
            print(f"mapper_used: {mapper_used}")
            print(f"cofactor_handling: {cofactor_handling}")
            print(f"mapping_error: {mapping_error if mapping_error is not None else ''}")
            if os.path.exists(error_output_file) and os.path.getsize(error_output_file) > 0:
                try:
                    process_error_reactions_smiles(error_output_file, compound_smiles_map)
                except Exception:
                    pass
            sys.exit(0)
        else:
            log_and_exit_early(f"No non-cofactor substrates found for reaction {selected_id}.", selected_id=selected_id,
                               reaction_no_stoich=raw_reaction)
    else:
        log_and_exit_early(f"No non-cofactor substrates found for reaction {selected_id}.", selected_id=selected_id,
                           reaction_no_stoich=raw_reaction)
if args.highest_CCC:
    mol_idx = max(range(len(substrate_mols)), key=selection_key_highest)
    reason = 'highest chiral center count'
elif args.lowest_CCC:
    mol_idx = min(range(len(substrate_mols)), key=selection_key_lowest)
    reason = 'lowest chiral center count'
else:
    mol_idx = max(range(len(substrate_mols)), key=lambda i: substrate_mols[i].GetNumHeavyAtoms())
    reason = 'highest heavy atom count (default)'
mol = substrate_mols[mol_idx]
atom_idx = max(range(mol.GetNumAtoms()), key=lambda i: mol.GetAtomWithIdx(i).GetDegree())
substrate_mols[mol_idx] = label_atom(substrate_mols[mol_idx], atom_idx, 1)
product_cid, atom_mapped_rxn, tracking_successful, mapper_used, mapping_error = perform_atom_tracking(
    substrate_mols, product_smiles, mol_idx, atom_idx, selected_id, "[INITIAL]"
)
if not tracking_successful and mapper_used == "localmapper_fallback":
    degrees = [(i, substrate_mols[mol_idx].GetAtomWithIdx(i).GetDegree()) for i in
               range(substrate_mols[mol_idx].GetNumAtoms())]
    degrees.sort(key=lambda x: x[1], reverse=True)
    if len(degrees) > 1:
        alt_atom_idx = degrees[1][0]
        alt_product_cid, alt_atom_mapped_rxn, alt_tracking_successful, alt_mapper_used, alt_mapping_error = perform_atom_tracking(
            substrate_mols, product_smiles, mol_idx, alt_atom_idx, selected_id, "[ALTERNATIVE]"
        )
        if alt_tracking_successful:
            product_cid = alt_product_cid
            atom_mapped_rxn = alt_atom_mapped_rxn
            tracking_successful = True
            atom_idx = alt_atom_idx
            mapper_used = alt_mapper_used
            mapping_error = alt_mapping_error
if mapped_rxns and args.visualize:
    try:
        visualize_mapped_rxn(mapped_rxns, selected_id, save_images=True)
    except Exception:
        pass
if not tracking_successful and mapping_error is None:
    mapping_error = "Atom tracking failed - atom not found in products"
selected_mol = substrate_mols[mol_idx]
selected_clean = remove_atom_map_numbers(selected_mol)
substrate_cid = lookup_compound_id(selected_clean)
cofactor_handling = (
    'static' if args.static_cofactors and not args.dynamic_cofactors else
    'dynamic' if args.dynamic_cofactors and not args.static_cofactors else
    'both' if args.static_cofactors and args.dynamic_cofactors else
    'none'
)
if dynamic_skipped:
    cofactor_handling += '_skipped'
print(f"reaction_id: {selected_id}")
print(f"tracking: {substrate_cid} => {product_cid}")
print(f"mapped_rxns: {mapped_rxns if mapped_rxns is not None else ''}")
print(f"reaction_no_stoich: {reaction_no_stoich}")
print(f"reaction_after_cofactor: {reaction_after_cofactor if reaction_after_cofactor is not None else ''}")
print(f"selection_method: {reason}")
print(f"mapper_used: {mapper_used}")
print(f"cofactor_handling: {cofactor_handling}")
print(f"mapping_error: {mapping_error if mapping_error is not None else ''}")
if os.path.exists(error_output_file) and os.path.getsize(error_output_file) > 0:
    try:
        process_error_reactions_smiles(error_output_file, compound_smiles_map)
    except Exception:
        pass
