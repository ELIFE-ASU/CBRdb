import argparse
import os
import sys

import pandas as pd
from rdkit import Chem
from rxnmapper import RXNMapper

import CBRdb

if __name__ == "__main__":
    # lookup_compound_id must be defined after df_canonical_index is available
    lookup_compound_id = None

    script_dir = os.path.dirname(os.path.abspath(__file__))

    parser = argparse.ArgumentParser()
    parser.add_argument('--row_index', type=int, required=True)
    parser.add_argument('--static_cofactors', action='store_true', help='Ignores Cofactors based on a static list')
    parser.add_argument('--dynamic_cofactors', action='store_true', help='Ignores Compounds present on both sides.')
    parser.add_argument('--highest_CCC', action='store_true',
                        help='Select substrate with highest number of chiral centers')
    parser.add_argument('--lowest_CCC', action='store_true',
                        help='Select substrate with lowest number of chiral centers')
    parser.add_argument('--visualize', action='store_true',
                        help='If set, visualize the mapped reaction for this row index')
    args = parser.parse_args()

    reactions_file = 'CBRdb_R.csv'
    molecule_reference_file = 'CBRdb_C.csv'
    canonical_index_file = os.path.join(script_dir, 'CBRdb_CanonicalIndex.csv')
    output_file = os.path.join(script_dir, 'CBRdb_AtomTracking.csv')

    try:
        reactions_file_path = CBRdb.find_file(reactions_file)
        molecule_ref_path = CBRdb.find_file(molecule_reference_file)
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

    df_canonical_index = CBRdb.get_or_build_canonical_index(df_molecules, canonical_index_file)

    selected_row = df_input_reactions.iloc[args.row_index]
    selected_id = selected_row['id']
    raw_reaction = selected_row['reaction']

    rxn_mapper = RXNMapper()

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

    if args.dynamic_cofactors:
        substrate_ids, product_ids = CBRdb.omit_shared_compounds(substrate_ids, product_ids)

    substrate_smiles = '.'.join(
        [compound_smiles_map.get(cid, '') for cid in substrate_ids if cid in compound_smiles_map])
    product_smiles = '.'.join([compound_smiles_map.get(cid, '') for cid in product_ids if cid in compound_smiles_map])

    # -------------------- Filter Cofactors --------------------
    static_cofactors = args.static_cofactors

    # --- Remove stoichiometry numbers from reaction string ---

    # Save reaction after removing stoichiometry
    reaction_no_stoich = CBRdb.remove_stoichiometry(raw_reaction)

    # Show cofactor-filtered reaction if either static or dynamic cofactor flag is set
    remove_cofactors_flag = args.static_cofactors or args.dynamic_cofactors
    reaction_after_cofactor = CBRdb.get_cofactor_filtered_reaction(substrate_ids, product_ids, cofactor_ids,
                                                                   remove_cofactors_flag)

    if static_cofactors:
        cofactor_smiles = set()
        for cid in cofactor_ids:
            if cid in compound_smiles_map:
                canonical = CBRdb.canonicalize_smiles(compound_smiles_map[cid])
                if canonical:
                    cofactor_smiles.add(canonical)


        def filter_cofactors(smiles_string):
            if not smiles_string:
                return ""
            return '.'.join([
                smi for smi in smiles_string.strip().split('.')
                if CBRdb.canonicalize_smiles(smi) not in cofactor_smiles
            ])


        substrate_smiles = filter_cofactors(substrate_smiles)
        product_smiles = filter_cofactors(product_smiles)
    substrate_mols = CBRdb.smiles_to_mols(substrate_smiles)

    if not substrate_mols:
        print(f"No non-cofactor substrates found for reaction {selected_id}. Halting.")
        sys.exit(0)

    if args.highest_CCC:
        # Highest chiral center count, break ties with highest heavy atom count
        mol_idx = max(range(len(substrate_mols)), key=CBRdb.selection_key_highest)
        reason = 'highest chiral center count'
    elif args.lowest_CCC:
        # Lowest chiral center count, break ties with highest heavy atom count
        mol_idx = min(range(len(substrate_mols)), key=CBRdb.selection_key_lowest)
        reason = 'lowest chiral center count'
    else:
        mol_idx = max(range(len(substrate_mols)), key=lambda i: substrate_mols[i].GetNumHeavyAtoms())
        reason = 'highest heavy atom count (default)'

    mol = substrate_mols[mol_idx]
    atom_idx = max(range(mol.GetNumAtoms()), key=lambda i: mol.GetAtomWithIdx(i).GetDegree())
    substrate_mols[mol_idx] = CBRdb.label_atom(substrate_mols[mol_idx], atom_idx, 1)

    print(f"[DEBUG] Selected molecule index: {mol_idx}, atom index: {atom_idx}, selection reason: {reason}")

    # -------------------- Perform Mapping --------------------
    reactant_smiles = Chem.MolToSmiles(substrate_mols[mol_idx], canonical=False)
    reaction_smiles = f"{reactant_smiles}>>{product_smiles}"

    atom_mapped_rxn = None  # Ensure variable is always defined
    try:
        mapped_rxn = rxn_mapper.get_attention_guided_atom_maps([reaction_smiles])[0]['mapped_rxn']
        atom_mapped_rxn = mapped_rxn  # Always set this
        _, product_part = mapped_rxn.split('>>')
        product_mols = CBRdb.smiles_to_mols(product_part)
        matched = CBRdb.trace_atom(1, product_mols)

        if matched:
            clean_smiles = CBRdb.remove_atom_map_numbers(matched)
            product_cid = CBRdb.lookup_compound_id(clean_smiles)
            print(f"[DEBUG] Matched Product SMILES: {clean_smiles} → {product_cid}")
        else:
            product_cid = 'NotFound'
            print(f"[DEBUG] Atom not found in any product.")

        # Visualize if flag is set
        if args.visualize:
            CBRdb.visualize_mapped_rxn(atom_mapped_rxn, selected_id, save_images=True)
    except Exception as e:
        print(f"Mapping error for {selected_id}: {e}")
        product_cid = 'Error'
        atom_mapped_rxn = 'Error'

    # -------------------- Get Substrate ID --------------------
    selected_mol = substrate_mols[mol_idx]
    selected_clean = CBRdb.remove_atom_map_numbers(selected_mol)
    substrate_cid = CBRdb.lookup_compound_id(selected_clean)
    print(f"[DEBUG] Selected Substrate SMILES: {selected_clean} → {substrate_cid}")

    # -------------------- Save Output --------------------
    output_row = pd.DataFrame([{
        'reaction_id': selected_id,
        'tracking': f"{substrate_cid} => {product_cid}",
        'atom_mapped_rxn': atom_mapped_rxn,
        'reaction_no_stoich': reaction_no_stoich,
        'reaction_after_cofactor': reaction_after_cofactor if reaction_after_cofactor is not None else '',
        'selection_method': reason,
        'cofactor_handling': (
            'static' if args.static_cofactors and not args.dynamic_cofactors else
            'dynamic' if args.dynamic_cofactors and not args.static_cofactors else
            'both' if args.static_cofactors and args.dynamic_cofactors else
            'none'
        )
    }])

    if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
        df_existing = pd.read_csv(output_file)
        df_combined = pd.concat([df_existing, output_row], ignore_index=True)
        df_combined = df_combined.sort_values(by='reaction_id').reset_index(drop=True)
        df_combined.to_csv(output_file, index=False)
    else:
        output_row.to_csv(output_file, index=False)

    print(f"[INFO] Saved tracking result to {output_file}")
