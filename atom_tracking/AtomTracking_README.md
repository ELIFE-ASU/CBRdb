# Chemical Reaction Atom Tracking System - Complete Guide

## Table of Contents
1. [Overview](#overview)
2. [What This Code Does](#what-this-code-does)
3. [System Architecture](#system-architecture)
4. [Step-by-Step Workflow](#step-by-step-workflow)
5. [Key Components](#key-components)
6. [Input Requirements](#input-requirements)
7. [Command-Line Options](#command-line-options)
8. [Output Results](#output-results)
9. [Error Handling](#error-handling)
10. [Technical Details](#technical-details)
11. [Usage Examples](#usage-examples)

---

## Overview

This is an **automated atom mapping system** for chemical reactions. It answers a fundamental question in chemistry: *"Where does a specific atom in a reactant molecule end up after a chemical reaction?"*

### Real-World Application
Imagine you have a reaction like:
```
Glucose + ATP → Glucose-6-phosphate + ADP
```

This system can:
- Select a specific atom in glucose (e.g., a carbon with the most bonds)
- Track where that exact atom ends up in the products
- Handle complex reactions with multiple reactants and products
- Remove "cofactor" molecules that don't participate in the main transformation

---

## What This Code Does

### Primary Function
The script processes **one chemical reaction at a time** and performs these tasks:

1. **Reads a reaction** from a CSV file (like `CBRdb_R.csv`)
2. **Converts compound IDs to chemical structures** (SMILES format)
3. **Selects one substrate molecule** based on your criteria
4. **Picks a specific atom** in that molecule (usually the most connected one)
5. **Maps atoms** across the reaction using known mapping algorithms
6. **Tracks the selected atom** to find which product it becomes part of
7. **Reports the results** showing the transformation path

### Key Innovation
It uses **two different atom mapping tools**:
- **RXNMapper**: Primary AI-based mapper (fast, accurate, but has limits)
- **localmapper**: Backup for very large reactions that exceed RXNMapper's token limit

---

## System Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                    INPUT DATA FILES                         │
├─────────────────────────────────────────────────────────────┤
│ 1. CBRdb_R.csv       - Reaction equations                	  │
│ 2. CBRdb_C.csv       - Compound ID → SMILES mapping         │
└─────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────┐
│              PREPROCESSING & CANONICALIZATION               │
├─────────────────────────────────────────────────────────────┤
│ • Remove stoichiometric coefficients (2, n, m+1, etc.)      │
│ • Build canonical SMILES index for fast lookup              │
│ • Parse reaction sides (substrates ⇄ products)			  │
└─────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────┐
│                  COFACTOR HANDLING (Optional)               │
├─────────────────────────────────────────────────────────────┤
│ Static:  Remove known cofactors (ATP, NADH, water, etc.)    │
│ Dynamic: Remove compounds appearing on both sides           │
└─────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────┐
│                    SUBSTRATE SELECTION                      │
├─────────────────────────────────────────────────────────────┤
│ Choose ONE substrate molecule based on:                     │
│ • Highest chiral centers count (--highest_CCC)              │
│ • Lowest chiral centers count (--lowest_CCC)                │
│ • Most heavy atoms (default)                                │
└─────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────┐
│                      ATOM SELECTION                         │
├─────────────────────────────────────────────────────────────┤
│ Pick the atom with highest degree (most bonds)              │
│ Label it with atom map number = 1                           │
└─────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────┐
│                      ATOM MAPPING                           │
├─────────────────────────────────────────────────────────────┤
│ Try RXNMapper first                                         │
│   ↓ (if token limit exceeded)                               │
│ Fallback to localmapper                                     │
│   ↓                                                         │
│ Result: Reaction with all atoms numbered                    │
└─────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────┐
│                     ATOM TRACKING                           │
├─────────────────────────────────────────────────────────────┤
│ 1. Find atom #1 in reactants                                │
│ 2. Look for same number in products                         │
│ 3. Identify which product molecule contains it              │
│ 4. Convert back to compound ID                              │
└─────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────┐
│                    OUTPUT RESULTS                           │
├─────────────────────────────────────────────────────────────┤
│ Success → Print tracking results to console                 │
│ Failure → Save to "Error Reactions.csv"                     │
│ Optional: Generate reaction visualization image             │
└─────────────────────────────────────────────────────────────┘
```

---

## Step-by-Step Workflow

### Step 1: Initialize and Load Data

**What happens:**
```python
# Command-line arguments are parsed
args = parser.parse_args()

# Required files are located
reactions_file = 'CBRdb_R.csv'             # Contains reaction equations
molecule_reference_file = 'CBRdb_C.csv'    # Contains compound structures
```

**Example data:**
```csv
# CBRdb_R.csv
id,reaction
R00001,C00031 + C00002 <=> C00668 + C00008 + C00009

# CBRdb_C.csv
compound_id,smiles
C00031,C(C(=O)O)O[C@H]1[C@@H](C(=O)[O-])O[C@@H](C(=O)[O-])[C@@H]1O
C00002,C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-])O)O)N
```

**Result:** DataFrames `df_input_reactions` and `df_molecules` are populated.

---

### Step 2: Select the Specific Reaction

**What happens:**
```python
# User specifies which row to process
selected_row = df_input_reactions.iloc[args.row_index]
selected_id = selected_row['id']           # e.g., "R00001"
raw_reaction = selected_row['reaction']    # e.g., "C00031 + C00002 <=> C00668 + ..."
```

**The reason:** The system processes ONE reaction per execution. For batch processing, you'd run the code multiple times (e.g., in parallel using a job scheduler).

---

### Step 3: Remove Stoichiometric Coefficients

**What happens:**
```python
reaction_no_stoich = remove_coefficients(raw_reaction)
```

**Example transformation:**
```
Before: 2 C00031 + n C00002 <=> m+1 C00668 + C00008
After:  C00031 + C00002 <=> C00668 + C00008
```

**The reason:** Coefficients like `2`, `n`, `m+1` confuse the atom mapping algorithms. We only care about which molecules participate, not their quantities.

---

### Step 4: Parse Reaction Sides

**What happens:**
```python
substrates_raw, products_raw = raw_reaction.split('<=>')
substrate_ids = ['C00031', 'C00002']  # Extract compound IDs
product_ids = ['C00668', 'C00008', 'C00009']
```

**Result:** Two lists of compound IDs representing the left (substrates) and right (products) sides of the reaction.

---

### Step 5: Handle Cofactors (Optional)

This is where the "helper" molecules that don't participate in the main transformation are filtered out.

#### Option A: Static Cofactor List (`--static_cofactors`)
```python
# Predefined list of common cofactors
cofactor_ids = ['C00002', 'C00003', 'C00004', ...]  # ATP, NAD+, NADP+, H2O, etc.
```

**Example:**
```
Before: Glucose + ATP ⇄ Glucose-6-P + ADP
After:  Glucose ⇄ Glucose-6-P
```

#### Option B: Dynamic Cofactor Detection (`--dynamic_cofactors`)
```python
# Remove any compound that appears on BOTH sides
shared = set(substrate_ids) & set(product_ids)
```

**Example:**
```
Before: 2-phosphoglycerate + H2O ⇄ phosphoenolpyruvate + H2O
After:  2-phosphoglycerate ⇄ phosphoenolpyruvate
```

**Safety check:** If removing cofactors would empty either side, the filtering is skipped and marked as `dynamic_skipped`.

---

### Step 6: Convert Compound IDs to SMILES

**What happens:**
```python
# Look up each compound's chemical structure
compound_smiles_map = {
    'C00031': 'C(C(=O)O)O[C@H]1[C@@H](C(=O)[O-])...',
    'C00002': 'C1=NC(=C2C(=N1)N(C=N2)...',
    ...
}

substrate_smiles = '.'.join([compound_smiles_map['C00031'], ...])
# Result: "SMILES1.SMILES2.SMILES3"  (dot-separated)
```

**SMILES notation:** A text-based way to represent molecular structures. Example:
- `C` = carbon atom
- `O` = oxygen atom
- `C(=O)O` = carboxylic acid group (-COOH)
- `@` = stereochemistry marker (3D orientation)

---

### Step 7: Build Canonical Index (First Run Only)

**What happens:**
```python
df_canonical_index = get_or_build_canonical_index(df_molecules, canonical_index_file)
```

**Purpose:** Create a lookup table that maps canonical SMILES → compound IDs. This speeds up reverse lookups when we find the product.

**Example entry:**
```csv
compound_id,original_smiles,fragment_smiles,canonical_smiles
C00031,"CCO.OCC","CCO","CCO"
C00031,"CCO.OCC","OCC","CCO"
```

**Note:** If a molecule has multiple fragments (separated by `.`), each fragment gets its own entry.

---

### Step 8: Select the Substrate Molecule

Multiple substrates might exist. We need to pick ONE to track.

**Selection strategies:**

#### Default: Highest Heavy Atom Count
```python
mol_idx = max(range(len(substrate_mols)), 
              key=lambda i: substrate_mols[i].GetNumHeavyAtoms())
```
Selects the biggest molecule (most carbon, nitrogen, oxygen, etc.).

#### `--highest_CCC`: Highest Chiral Center Count
```python
mol_idx = max(range(len(substrate_mols)), key=selection_key_highest)
```
Selects the molecule with the most stereogenic centers (important for tracking chirality).

#### `--lowest_CCC`: Lowest Chiral Center Count
```python
mol_idx = min(range(len(substrate_mols)), key=selection_key_lowest)
```
Selects the simplest molecule in terms of stereochemistry.

**Example:**
```
Substrates: [Glucose (C6H12O6), ATP (C10H16N5O13P3)]
Heavy atoms: [24, 31]
Selected: ATP (has more heavy atoms)
```

---

### Step 9: Select the Specific Atom

**What happens:**
```python
mol = substrate_mols[mol_idx]
atom_idx = max(range(mol.GetNumAtoms()), 
               key=lambda i: mol.GetAtomWithIdx(i).GetDegree())
```

**Atom degree:** Number of bonds an atom has.
- Degree 1: Terminal atom (like -OH, -CH3)
- Degree 2: Chain atom
- Degree 3: Branching point
- Degree 4: Highly connected (like central carbons in rings)

**Why pick the highest degree?** These atoms are typically:
- More central to the molecule's structure
- More likely to be preserved in products
- More chemically significant

**Example:**
```
Glucose: C-C-C-C-C-C-O
         ↑ (this carbon has degree 3, selected)
```

**Action:** The selected atom is labeled with map number `1`:
```python
substrate_mols[mol_idx] = label_atom(substrate_mols[mol_idx], atom_idx, 1)
```

---

### Step 10: Perform Atom Mapping

This is the **core algorithmic step** where atoms across reactants and products are matched.

#### Input to Mapper:
```
Reactant: [CH3:1]CH2OH (atom map numbers only on selected atom)
Products: CH3CHO.H2O
```

#### Process:
```python
def perform_atom_mapping(reaction_smiles, reaction_id):
    try:
        # Try RXNMapper first (transformer-based AI model)
        mapped_rxn = rxn_mapper.get_attention_guided_atom_maps([reaction_smiles])[0]['mapped_rxn']
        return mapped_rxn, "RXNMapper", None
    except Exception as e:
        # Check if error is due to token limit (>512 tokens)
        if is_token_limit_error:
            # Fallback to localmapper
            mapped_rxn = local_mapper.get_atom_map(reaction_smiles)
            return mapped_rxn, "localmapper_fallback", None
```

#### Output from Mapper:
```
[CH3:1][CH2:2][OH:3]>>[CH3:1][CH:2]=[O:3].[H:4][OH:3]
```

Every atom now has a number! Same numbers = same atom across the reaction.

**Why two mappers?**
- **RXNMapper**: State-of-the-art, but limited to reactions with <512 tokens (about 30-50 atoms)
- **localmapper**: Handles larger reactions but may be less accurate

---

### Step 11: Track the Atom Through Products

**What happens:**
```python
# 1. Find which atom was labeled as #1 in reactants
target_map_num = 1

# 2. Look for atom #1 in products
matched = trace_atom(target_map_num, product_mols)

# 3. Identify which product molecule contains it
if matched:
    clean_smiles = remove_atom_map_numbers(matched)
    product_cid = lookup_compound_id(clean_smiles)
```

**Matching strategy (in priority order):**

1. **Exact match:** Find atom with map number `1` in products
2. **Element + degree match:** If #1 missing, find same element with same connectivity
3. **Element match:** Last resort - match by element type only

**Fallback mechanism:**
If tracking fails with the primary atom and localmapper was used, the system tries the **second-most connected atom**:
```python
if not tracking_successful and mapper_used == "localmapper_fallback":
    # Try alternative atom
    alt_atom_idx = degrees[1][0]  # Second highest degree
    alt_product_cid, alt_atom_mapped_rxn, ... = perform_atom_tracking(...)
```

**Example:**
```
Mapped reaction: [C:1]H3CH2OH >> [C:1]H3CHO.H2O
Atom #1 is in: CH3CHO
Lookup: CH3CHO → C00084 (acetaldehyde)
Result: tracking = "C00031 => C00084"
```

---

### Step 12: Handle Edge Cases

#### Case 1: Variable Coefficient Reactions
```
Reaction: n C00001 <=> m C00001 + C00002
```
If a compound appears on both sides AND has variable coefficients (`n`, `m`, `x`), it's treated as a **self-reaction**:
```python
if is_variable_coeff_reaction and shared:
    main_cid = list(shared)[0]
    substrate_cid = main_cid
    product_cid = main_cid  # Same compound
    reason = 'special_variable_coeff'
```

#### Case 2: Invalid SMILES
If no valid substrate molecules can be created:
```python
if not substrate_mols:
    log_and_exit_early("No valid substrate SMILES found", ...)
```
This saves error details to `Error Reactions.csv`.

---

### Step 13: Generate Output

**Console output:**
```
reaction_id: R00001
tracking: C00031 => C00668
mapped_rxns: [C:1]CC(O)...>>[C:1]CC(=O)...
reaction_no_stoich: C00031 + C00002 <=> C00668 + C00008
reaction_after_cofactor: C00031 <=> C00668
selection_method: highest heavy atom count (default)
mapper_used: RXNMapper
cofactor_handling: static
mapping_error: 
```

**Fields explained:**
- `reaction_id`: Unique identifier for the reaction
- `tracking`: **Main result** - shows substrate compound → product compound transformation
- `mapped_rxns`: Full atom-mapped reaction SMILES (for verification)
- `reaction_no_stoich`: Cleaned reaction without coefficients
- `reaction_after_cofactor`: Reaction after removing cofactors
- `selection_method`: How the substrate was chosen
- `mapper_used`: Which algorithm performed the mapping
- `cofactor_handling`: Which cofactor strategy was applied
- `mapping_error`: Any errors encountered (empty if successful)

---

### Step 14: Visualization (Optional)

**If `--visualize` flag is set:**
```python
visualize_mapped_rxn(mapped_rxns, selected_id, save_images=True)
```

**Result:** Saves a PNG image to `Reactions Visualizations/R00001.png` showing:
- All reactant molecules with atoms numbered
- Arrow
- All product molecules with atoms numbered
- Color-coded atom correspondences

---

## Key Components

### 1. RDKit (Molecular Processing)
**What it does:**
- Parses SMILES strings into molecular objects
- Manipulates atom properties (adding map numbers)
- Canonicalizes structures for comparison
- Counts chiral centers and heavy atoms

**Key functions:**
- `Chem.MolFromSmiles()`: Text → Molecule
- `Chem.MolToSmiles()`: Molecule → Text
- `atom.SetAtomMapNum()`: Label atoms
- `Chem.FindMolChiralCenters()`: Find stereogenic centers

### 2. RXNMapper (AI Atom Mapping)
**Technology:** Transformer-based neural network trained on millions of reactions

**How it works:**
1. Tokenizes reaction SMILES
2. Uses attention mechanism to identify atom correspondences
3. Outputs reaction with all atoms numbered

**Limitation:** Maximum 512 tokens (~30-50 heavy atoms)

### 3. localmapper (Fallback Mapping)
**Technology:** Graph-based algorithm

**How it works:**
1. Represents molecules as graphs (atoms = nodes, bonds = edges)
2. Uses maximum common substructure to find atom matches
3. Handles larger reactions than RXNMapper

**Trade-off:** More robust to size, potentially less accurate on complex rearrangements

### 4. FileLock (Thread Safety)
**Purpose:** Allows multiple instances of the code to run simultaneously without corrupting output files

**How it works:**
```python
with FileLock(lock_file, timeout=60):
    # Only one process can execute this block at a time
    df_combined.to_csv(output_file, index=False)
```

### 5. Canonical Index
**Purpose:** Fast reverse lookup from SMILES → compound ID

**Without index:**
- Need to canonicalize EVERY molecule in database for each lookup
- O(N) complexity per lookup

**With index:**
- Pre-computed canonical forms
- O(1) hash lookup
- Stored in `CBRdb_CanonicalIndex.csv`

---

## Input Requirements

### Required Files

#### 1. `CBRdb_R.csv`
**Columns:**
- `id`: Reaction identifier (e.g., "R00001")
- `reaction`: Equation with compound IDs (e.g., "C00031 + C00002 <=> C00668 + C00008")

**Format rules:**
- Must use `<=>` as separator between reactants and products
- Can include coefficients (2, n, m+1, etc.)
- Use `+` to separate multiple compounds

#### 2. `CBRdb_C.csv`
**Columns:**
- `compound_id`: Unique identifier (e.g., "C00031")
- `smiles`: Chemical structure in SMILES notation

**Requirements:**
- Every compound in `CBRdb_R.csv` must have an entry
- SMILES can contain multiple fragments (dot-separated)
- SMILES should be valid RDKit-parseable format

### File Location
The code tries to find the files in:
1. Same directory as the script
2. One directory up (`../`)
3. Two directories up (`../../`)

---

## Command-Line Options

### Required Arguments

#### `--row_index` (integer)
**What it does:** Specifies which row (0-indexed) from `CBRdb_R.csv` to process

**Example:**
```bash
python script.py --row_index 0      # Process first reaction
python script.py --row_index 42     # Process 43rd reaction
```

### Cofactor Handling Flags

#### `--static_cofactors`
**What it does:** Removes 51 predefined common cofactors

**Cofactor list includes:**
- Energy carriers: ATP, ADP, AMP (C00002, C00008, C00020)
- Electron carriers: NAD+, NADH, NADP+, NADPH (C00003, C00004, C00005, C00006)
- Small molecules: H2O, O2, CO2, NH3, Pi (C00001, C00007, C00011, C00014, C00009)
- Coenzyme A derivatives
- And 36 more...

**When to use:** 
- Studying core metabolic transformations
- Want to focus on carbon skeleton changes
- Tracking main substrate-product relationships

#### `--dynamic_cofactors`
**What it does:** Removes any compound appearing on BOTH sides of the reaction

**Example:**
```
Before: A + B + H2O ⇄ C + D + H2O
After:  A + B ⇄ C + D
```

**When to use:**
- Reactions with solvents or catalysts
- Want to auto-detect "spectator" molecules
- Unknown cofactor set

**Can combine:** `--static_cofactors --dynamic_cofactors` uses both strategies

### Substrate Selection Flags

#### `--highest_CCC`
**What it does:** Selects substrate with **most chiral centers**

**Example:**
```
Substrates:
- Glucose (5 chiral centers)
- Pyruvate (0 chiral centers)

Selected: Glucose
```

**When to use:**
- Studying stereochemical transformations
- Tracking chirality through reactions
- Focus on complex molecules

#### `--lowest_CCC`
**What it does:** Selects substrate with **fewest chiral centers**

**When to use:**
- Studying simple molecules
- Avoiding complex stereochemistry
- Focus on achiral transformations

**Default (no flag):** Selects substrate with **most heavy atoms** (largest molecule)

### Visualization Flag

#### `--visualize`
**What it does:** Generates PNG image of the mapped reaction

**Output location:** `Reactions Visualizations/[reaction_id].png`

**Image content:**
- All molecules drawn as 2D structures
- Atom map numbers shown on atoms
- Color-coded to show correspondences
- High resolution (300 DPI)

**When to use:**
- Verifying atom mapping accuracy
- Creating figures for papers/presentations
- Debugging mapping errors

---

## Output Results

### Success Output

**Console format:**
```
reaction_id: R00001
tracking: C00031 => C00668
mapped_rxns: [C:1](C(=O)[O-])O[C@H]1...>>[C:1](C(=O)[O-])O...
reaction_no_stoich: C00031 + C00002 <=> C00668 + C00008 + C00009
reaction_after_cofactor: C00031 <=> C00668
selection_method: highest heavy atom count (default)
mapper_used: RXNMapper
cofactor_handling: static
mapping_error: 
```
**Important Note:**
The code makes a text-based file for each mapped reaction. For a single .csv file containing all the mapped reactions, we suggest
using combine_out_files.py. (can be found in the https://github.com/ELIFE-ASU/CBRdb/tree/main/atom_tracking)
It is worth noting that the folder for the output and error files can be modified through the shell script.

### Error Output

**If tracking fails, saved to:** `Error Reactions.csv`

**Columns:**
- `reaction_id`: Which reaction failed
- `tracking`: Empty (tracking failed)
- `mapped_rxns`: May be present if mapping succeeded but tracking failed
- `reaction_no_stoich`: Cleaned reaction
- `reaction_after_cofactor`: Reaction after cofactor removal
- `selection_method`: How substrate was chosen
- `mapper_used`: Which mapper was attempted
- `cofactor_handling`: Cofactor strategy used
- `mapping_error`: Detailed error message
- `error_message`: Same as mapping_error
- `reaction_smiles`: Full reaction in SMILES format (added during post-processing)

**Common error types:**

1. **Token limit errors:**
```
mapping_error: RXNMapper: The size of tensor a (512) must match the size of tensor b (634)
mapper_used: rxnmapper_failed_no_fallback
```

2. **Atom not found:**
```
mapping_error: Atom tracking failed - atom not found in products
mapper_used: RXNMapper
```

3. **Invalid SMILES:**
```
error_message: No valid substrate SMILES found for reaction R12345
```

4. **Both mappers failed:**
```
mapping_error: RXNMapper: token limit; localmapper: graph error
mapper_used: both_failed
```

### Cofactor Handling Codes

**In output field `cofactor_handling`:**
- `none`: No cofactor filtering
- `static`: Only static cofactor list used
- `dynamic`: Only dynamic filtering used
- `both`: Both strategies combined
- `static_skipped`: Dynamic filtering would empty a side, reverted to static only
- `dynamic_skipped`: Dynamic filtering would empty a side, reverted to none
- `both_skipped`: Both would cause issues, no filtering applied

---

## Error Handling

### Automatic Error Recovery

#### 1. Mapper Fallback
```
RXNMapper fails (token limit)
    ↓
Try localmapper
    ↓
If successful, continue with tracking
```

#### 2. Alternative Atom Selection
```
Tracking fails with highest-degree atom (when using localmapper)
    ↓
Try second-highest degree atom
    ↓
If successful, use this tracking result
```

#### 3. Dynamic Cofactor Safety
```
Remove shared compounds
    ↓
Check if any side is now empty
    ↓
If yes, restore original compounds and mark as "skipped"
```

### Error File Processing

**Post-execution enhancement:**
```python
process_error_reactions_smiles(error_output_file, compound_smiles_map)
```

**What it does:**
- Reads `Error Reactions.csv`
- Adds `reaction_smiles` column with full SMILES representation
- Allows manual inspection of failed reactions in chemical structure format

---

## Technical Details

### SMILES Notation Basics

**Examples:**
```
Methane:       C
Ethanol:       CCO
Acetic acid:   CC(=O)O
Benzene:       c1ccccc1
Glucose:       C(C1C(C(C(C(O1)O)O)O)O)O
```

**Special characters:**
- `C`: Aliphatic carbon
- `c`: Aromatic carbon
- `()`: Branching
- `=`: Double bond
- `#`: Triple bond
- `[]`: Atom specification with properties
- `@`: Stereochemistry
- `.`: Fragment separator

### Atom Mapping Example

**Input reaction:**
```
CC(O)C >> CC(=O)C
(2-propanol >> acetone)
```

**After atom mapping:**
```
[C:1][C:2]([O:3])[C:4]>>[C:1][C:2](=[O:3])[C:4]
```

**Interpretation:**
- Carbon 1 stays as carbon 1 (methyl group)
- Carbon 2 stays as carbon 2 (central carbon)
- Oxygen 3 stays as oxygen 3 (but changes from single to double bond)
- Carbon 4 stays as carbon 4 (other methyl group)

### Canonicalization

**Purpose:** Create a unique representation for each molecule

**Example:**
```
Input SMILES:
- "CCO"     (ethanol written left-to-right)
- "OCC"     (ethanol written right-to-left)

Canonical form:
- "CCO"     (both convert to the same)
```

**Why it matters:** Allows reliable comparison and lookup

### Chiral Centers

**Definition:** An atom (usually carbon) bonded to four different groups

**Example:**
```
Alanine: C[C@H](N)C(=O)O
         ↑ This carbon is chiral (attached to: CH3, NH2, COOH, H)
```

**Notation:**
- `@`: Counterclockwise configuration (S)
- `@@`: Clockwise configuration (R)

---

## Usage Examples

### Example 1: Basic Atom Tracking

**Command:**
```bash
python AtomTracking_Enhanced.py --row_index 0
```

**Input (CBRdb_R.csv, row 0):**
```
id,reaction
R00001,C00031 + C00002 <=> C00668 + C00008 + C00009
```

**Process:**
1. Loads reaction R00001
2. Finds SMILES for C00031 (D-glucose 1-phosphate) and C00002 (ATP)
3. Selects ATP (has more heavy atoms)
4. Picks most connected phosphorus atom
5. Maps reaction with RXNMapper
6. Tracks phosphorus → finds it in C00009 (Pi, inorganic phosphate)

**Output:**
```
reaction_id: R00001
tracking: C00002 => C00009
mapped_rxns: [complex SMILES with atom numbers]
selection_method: highest heavy atom count (default)
mapper_used: RXNMapper
```

---

### Example 2: Focusing on Core Transformation

**Command:**
```bash
python AtomTracking_Enhanced.py --row_index 5 --static_cofactors
```

**Input:**
```
id,reaction
R00006,C00031 + C00002 <=> C00668 + C00008 + C00009
```

**Effect of `--static_cofactors`:**
- Removes C00002 (ATP), C00008 (ADP), C00009 (Pi)
- Leaves only: C00031 ⇄ C00668

**Output:**
```
tracking: C00031 => C00668
reaction_after_cofactor: C00031 <=> C00668
cofactor_handling: static
```

**Result:** Tracks the glucose phosphate transformation, ignoring energy transfer

---

### Example 3: Tracking Stereochemistry

**Command:**
```bash
python AtomTracking_Enhanced.py --row_index 10 --highest_CCC
```

**Input:**
```
id,reaction
R00011,C00022 + C00010 <=> C00036 + C00011
```
(Pyruvate + CoA ⇄ Acetyl-CoA + CO2)

**Effect of `--highest_CCC`:**
- C00022 (pyruvate): 0 chiral centers
- C00010 (CoA): 5+ chiral centers
- **Selects:** CoA

**Output:**
```
tracking: C00010 => C00036
selection_method: highest chiral center count
```

---

### Example 4: Large Reaction with Fallback

**Command:**
```bash
python AtomTracking_Enhanced.py --row_index 100
```

**Scenario:** Reaction with 60+ heavy atoms (exceeds RXNMapper limit)

**Process:**
1. RXNMapper attempts mapping
2. Fails with: "should be at most 512 tokens"
3. Automatically switches to localmapper
4. Successfully maps and tracks

**Output:**
```
mapper_used: localmapper_fallback
mapping_error: 
```

---

### Example 5: Generating Visualizations

**Command:**
```bash
python AtomTracking_Enhanced.py --row_index 25 --visualize
```

**Additional output:**
- File created: `Reactions Visualizations/R00025.png`
- Contains: Color-coded 2D molecular structures with atom numbers

**Use case:** Verify that atom mapping is correct by visual inspection

---

### Example 6: Batch Processing (Shell Script)

**Windows PowerShell:**
```powershell
# Process reactions 0-99 with static cofactor removal
for ($i=0; $i -lt 100; $i++) {
    python AtomTracking_Enhanced.py --row_index $i --static_cofactors
}
```

**Linux/Mac Bash:**
```bash
# Process all reactions in parallel (4 at a time)
parallel -j 4 python AtomTracking_Enhanced.py --row_index {} --static_cofactors ::: {0..99}
```

---

### Example 7: Combined Cofactor Strategies

**Command:**
```bash
python AtomTracking_Enhanced.py --row_index 50 --static_cofactors --dynamic_cofactors
```

**Effect:**
1. First removes 51 static cofactors
2. Then removes any remaining shared compounds
3. Maximum filtering for core transformations

**Output:**
```
cofactor_handling: both
```

---

## Expected Results Summary

### What You Should See After Running

#### Successful Tracking
✅ **Console output** with 9 fields showing complete tracking information  
✅ **Substrate → Product transformation** clearly identified  
✅ **Mapper used** (RXNMapper or localmapper_fallback)  
✅ **Optional PNG image** (if --visualize used)  

#### Failed Tracking
❌ **Entry added** to `Error Reactions.csv`  
❌ **Detailed error message** explaining what went wrong  
❌ **Script exits** with code 1  
---

## Troubleshooting

### Common Issues

#### Issue 1: "File not found"
**Error:** `FileNotFoundError: CBRdb_R.csv not found`

**Solution:**
- Ensure CSV files are in the same directory as script
- Or in parent directories (../ or ../../)
- Check file names are exact (case-sensitive on Linux/Mac)

---

#### Issue 2: "Row index out of bounds"
**Error:** `Row index 500 out of bounds.`

**Solution:**
```python
# Check how many reactions exist
import pandas as pd
df = pd.read_csv('CBRdb_R.csv')
print(f"Total reactions: {len(df)}")
print(f"Valid indices: 0 to {len(df)-1}")
```
---

#### Issue 3: "localmapper not available"
**Error:** `mapper_used: rxnmapper_failed_no_fallback`

**Solution:**
```bash
pip install localmapper
```
---

#### Issue 4: All tracking results are "NotFound"
**Possible causes:**
- Atom map numbers not preserved by mapper
- Product molecules not in canonical index
- Reaction is intramolecular (atom changes within same molecule)

**Debug steps:**
1. Use `--visualize` to inspect atom mapping
2. Check if `mapped_rxns` output has numbers
3. Verify products are in `CBRdb_C.csv`

---

#### Issue 5: "Dynamic filtering skipped"
**Message:** `cofactor_handling: dynamic_skipped`

**Explanation:** Removing shared compounds would empty reactants or products

**Example:**
```
Reaction: A + B <=> B + C
Shared: B
After removal: A <=> C  ✓ (OK)

Reaction: A <=> A + B
Shared: A
After removal: (empty) <=> B  ✗ (SKIPPED)
```
---

## Citation and Credits

**Core dependencies:**
- **RDKit:** Open-source cheminformatics toolkit
- **RXNMapper:** Schwaller, P., et al. (2021). "Extraction of organic chemistry grammar from unsupervised learning of chemical reactions"
- **localmapper:** Shuan, Chen, et al. (2024). "Precise atom-to-atom mapping for organic reactions via human-in-the-loop machine learning"
