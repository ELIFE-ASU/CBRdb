import os
import re

import pandas as pd
from rdkit import Chem as Chem
from rdkit import RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

# suppress SettingWithCopyWarning
pd.options.mode.chained_assignment = None


def scan_for_kegg_pattern(s, prefix):
    if prefix in s:
        return [i for i in s.replace(';', ' ').replace('[', ' ').replace(']', ' ').split()
                if i.startswith(prefix) and len(i) == len(prefix) + 5]
    return []


def make_lists_printable(df):
    list_like_cols = df.loc[:, df.iloc[0].apply(lambda x: isinstance(x, (set, list)))]
    printable_df = df.copy(deep=True)
    printable_df[list_like_cols.columns] = list_like_cols.map(' '.join)
    return printable_df


def scan_for_terms(s):
    for j in OK:
        if j in s or s.startswith(('R1', 'R0')):
            return s
    return None


def remove_terms(s):
    for j in OK:
        if j in s or s.startswith(('R1', 'R0')):
            return None
    return s


def flag_missing_reactions(x):
    missing_reactions = [i for i in x if i not in reactions.index]
    return float('NaN') if not missing_reactions else missing_reactions


def load_reactions_data(target_dir):
    a = os.listdir(os.path.abspath(target_dir))
    reactions = {}
    for i in a:
        try:
            with open(f"{target_dir}{i}/{i}.data", 'r') as f:
                entries = {}
                for line in f.readlines()[:-1]:
                    if not line.startswith(' '):
                        field = line.split(' ')[0]
                        entries[field] = line.lstrip(field).lstrip().rstrip('\n')
                    else:
                        entries[field] += '; ' + line.lstrip().rstrip('\n')
                reactions[i] = entries
        except:
            pass
    return pd.DataFrame(reactions).fillna('').T


def get_reaction_ids_substr(reactions, substr="incomplete reaction"):
    incomplete_reaction_ids = reactions[
        reactions['COMMENT'].str.contains(substr, case=False, na=False)].index.tolist()
    return incomplete_reaction_ids


def main(target_dir='../data/kegg_data_R/', data_dir="data/"):
    print("Input: An unzipped, locally-downloaded folder of raw KEGG reaction data files.", flush=True)
    print("Output: CSVs with relevant reaction info, ties to other databases and flags for removal.", flush=True)
    global reactions, OK

    f_save_files = False

    prefixes = {
        'RCLASS': r'RC\d{5}',
        'BRITE': r'br\d{5}',
        'PATHWAY': r'rn\d{5}',
        'MODULE': r'M\d{5}',
        'ORTHOLOGY': r'K\d{5}'
    }

    OK = ['STEP', 'REACTION', 'SIMILAR', 'TO', 'SEE', '+', r'R\d{5}']

    # This needs to point to the directory where the KEGG reaction data files are stored.
    reactions = load_reactions_data(target_dir)

    # Parse each reaction's ties to other databases.
    # Flag reactions with glycan participants and "overall" reactions
    # See [this table](https://www.kegg.jp/kegg/tables/br08210.html) for info
    reactions = reactions.assign(
        OVERALL_FLAG=reactions['ENTRY'].str.contains('Overall', case=False),
        GLYCAN_FLAG=reactions['EQUATION'].str.contains('G'),
        EQUATION_ENTRIES=reactions['EQUATION'].str.findall(r'C\d{5}|G\d{5}'),
        ENZYME_ENTRIES=reactions['ENZYME'].str.replace(';', ' ').str.split(),
        DBLINKS_RHEA_ENTRIES=reactions['DBLINKS'].str.lstrip('RHEA: ').str.split())

    for col, pattern in prefixes.items():
        reactions[f'{col}_ENTRIES'] = reactions[col].str.findall(pattern)

    # KEGG no longer maintains the RCLASS database, BUT some edges in KEGG pathways
    # are still tagged with RCLASS entries instead of their constituent reactions.
    # Let's reconstruct RCLASS from REACTION (since we get that for free).
    r_rclass = reactions['RCLASS'].str.split('; ').explode().to_frame()
    r_rclass = r_rclass.assign(RC=r_rclass['RCLASS'].str.extract(r'(RC\d{5})'),
                               CPAIR=r_rclass['RCLASS'].str.findall(r'(C\d{5}_C\d{5})')).drop('RCLASS', axis=1).explode(
        'CPAIR')

    if f_save_files:
        rclass_db = r_rclass.reset_index(names='REACTION').groupby(by=['RC', 'CPAIR'])['REACTION'].apply(set).apply(
            list).to_frame()
        make_lists_printable(rclass_db).to_csv(os.path.join(data_dir, 'rclass_db.csv'))

    # Use comments to identify multi-step reactions
    # artificial shortcuts compressing multiple existing reaction IDs into a single net reaction.
    to_replace = {'MULTISTEP': 'MULTI-STEP',
                  '3STEP': '3-STEP',
                  '; SEE R': ' (SEE R',
                  'REACTION;': 'REACTION '}
    dm = reactions['COMMENT'].replace('', float('NaN')).dropna().str.upper().reset_index()
    dm['ORIGINAL_COMMENT'] = dm['COMMENT'].copy(deep=True)
    for k, v in to_replace.items():
        dm['COMMENT'] = dm['COMMENT'].str.replace(k, v)
    dm['COMMENT'] = dm['COMMENT'].str.split(';')
    s1 = r'COMMENT.str.contains("STEP REACTION")'
    s2 = r'COMMENT.str.contains(r"R\d{5}\+|R\d{5} \+|R\d{5} \,")'
    s3 = r'COMMENT.str.contains("STEP OF|PART OF|THE LAST|THE FIRST|THE FORMER|THE LATTER|POSSIBLY|PROBABLY IDENTICAL")'
    dm = dm.explode('COMMENT').reset_index(drop=True).query(f'{s1} and {s2} and not {s3}')

    to_replace = {
        'TWO': '2',
        'THREE': '3',
        'FOUR': '4',
        'MULTI': 'N',
        ' + ': '+'
    }

    for k, v in to_replace.items():
        dm['COMMENT'] = dm['COMMENT'].str.replace(k, v).str.strip()

    dm = dm.sort_values('COMMENT')

    # Simplified version of the selected code
    format_to_remove = r'SIMILAR TO R\d{5}, R\d{5}\+R\d{5}'
    for _ in range(2, 15):
        to_remove = dm['COMMENT'].str.findall(format_to_remove + r'\)').explode().dropna()
        for i, str_to_remove in to_remove.items():
            dm.loc[i, 'COMMENT'] = dm.loc[i, 'COMMENT'].replace(str_to_remove, '')
        format_to_remove += r'\+R\d{5}'

    format_to_remove = r'SIMILAR TO R\d{5}\+R\d{5}'
    for _ in range(2, 15):
        to_remove = dm['COMMENT'].str.findall(format_to_remove + r'\)').explode().dropna()
        for i, str_to_remove in to_remove.items():
            dm.loc[i, 'COMMENT'] = dm.loc[i, 'COMMENT'].replace(str_to_remove, '')
        format_to_remove += r'\+R\d{5}'

    reaction_string_format = r'R\d{5}'
    dm = dm[dm['COMMENT'].str.contains(reaction_string_format)]

    dm['COMMENT_FILTERED'] = (
        dm['COMMENT']
        .str.replace(' OR ', ', ')
        .str.split()
        .explode()
        .apply(scan_for_terms)
        .dropna()
        .groupby(level=0)
        .apply(' '.join)
    )
    dm['OTHER_COMMENTS'] = (
        dm['COMMENT']
        .str.split()
        .explode()
        .apply(remove_terms)
        .dropna()
        .groupby(level=0)
        .apply(' '.join)
    )
    dm = dm.query('~COMMENT_FILTERED.str.contains("SIMILAR TO")')

    to_remove = ['SEE ', '(', ')', 'REACTION ', '-STEP']
    dm['info'] = dm['COMMENT_FILTERED'].copy()
    for i in to_remove:
        dm['info'] = dm['info'].str.replace(i, '').str.strip(', ')

    dm['n_steps_cited'] = dm['info'].str[0].replace('N', float('nan')).astype(float)

    removed = {}
    for i in dm.index:
        noself = (
            dm.loc[i, 'info'][1:]
            .replace(dm.loc[i, 'index'] + '+', '')
            .replace('+' + dm.loc[i, 'index'], '')
            .replace(dm.loc[i, 'index'], '')
        )
        if '+' in noself:
            removed[i] = [j.strip() for j in noself.split(',')]

    dm['steps_shown'] = removed
    dm = (
        dm.explode('steps_shown')
        .dropna(subset=['steps_shown'])
        .sort_index()
        .reset_index(drop=True)
    )
    dm['steps_shown'] = dm['steps_shown'].str.split('+')
    dm['n_steps_shown'] = dm['steps_shown'].apply(len)

    # Identify reactions with incomplete multistep shortcuts
    has_rns_missing = dm.loc[dm['steps_shown'].apply(flag_missing_reactions).dropna().index, 'index'].to_list()
    incomplete_multistep = dm.loc[dm['index'].isin(has_rns_missing), 'index'].unique()
    print('Reactions with incomplete multistep shortcuts:', len(incomplete_multistep), incomplete_multistep, flush=True)

    # Filter out incomplete multistep reactions
    dm = dm.query('index not in @incomplete_multistep')
    dm['step_group'] = dm['steps_shown'].apply('+'.join)

    # Flag multistep shortcuts and generate data files
    reactions['MULTISTEP_FLAG'] = reactions.index.isin(dm['index'].unique())
    reactions['STEP_ENTRIES'] = dm.groupby('index')['step_group'].apply(list).fillna('')
    reactions_printable = make_lists_printable(reactions)
    if f_save_files:
        make_lists_printable(dm).to_csv('../reactions_multistep_intermediate_processing.csv.zip',
                                        compression='zip',
                                        encoding='utf-8')
        reactions_printable.to_csv(os.path.join(data_dir, 'reactions_processed_full.csv.zip'),
                                   compression='zip',
                                   encoding='utf-8')

    reactions_processed = pd.concat(
        [reactions_printable.filter(like='FLAG'),
         reactions_printable.filter(like='ENTRIES')],
        axis=1)
    if f_save_files:
        reactions_processed.to_csv(os.path.join(data_dir, 'reactions_processed_basic.csv.zip'),
                                   compression='zip',
                                   encoding='utf-8')

    # Select only the reactions that are overall_flagged or glycan_flagged
    reactions_shortcut = reactions_printable.query('OVERALL_FLAG == True')
    print('Reactions that are shortcuts:', len(reactions_shortcut), flush=True)
    reactions_shortcut = reactions_shortcut['ENTRY'].to_list()
    reactions_shortcut = list(set([i.split()[0].strip() for i in reactions_shortcut]))

    reactions_glycan = reactions_printable.query('GLYCAN_FLAG == True')
    print('Reactions with glycan participants:', len(reactions_glycan), flush=True)
    reactions_glycan = reactions_glycan['ENTRY'].to_list()
    reactions_glycan = list(set([i.split()[0].strip() for i in reactions_glycan]))

    reactions_incomplete = get_reaction_ids_substr(reactions, substr="incomplete reaction")
    print('Reactions with incomplete data:', len(reactions_incomplete), flush=True)

    reactions_general = get_reaction_ids_substr(reactions, substr="general reaction")
    print('General reactions:', len(reactions_general), flush=True)

    data = {
        'Reaction': reactions_shortcut + reactions_glycan + reactions_incomplete + reactions_general,
        'Type': ['shortcut'] * len(reactions_shortcut)
                + ['glycan'] * len(reactions_glycan)
                + ['incomplete'] * len(reactions_incomplete)
                + ['general'] * len(reactions_general)
    }
    data = pd.DataFrame(data)
    # Sort the data
    data = data.sort_values(by=['Reaction'])
    data.to_csv(os.path.join(data_dir, 'R_IDs_bad.dat'), index=False)


def load_bad_entries(bad_file, target_str="molless"):
    with open(bad_file, 'r') as file:
        return [line.split(',')[0].strip() for line in file if target_str in line]


def get_reactions_with_substring(reactions_df, substring):
    return reactions_df[reactions_df['reaction'].str.contains(substring, case=False, na=False)]


def standardize_mol(mol):
    # Standardize the molecule
    mol.UpdatePropertyCache(strict=False)
    Chem.SetConjugation(mol)
    Chem.SetHybridization(mol)
    # Normalize the molecule
    Chem.SanitizeMol(mol, sanitizeOps=(Chem.SANITIZE_ALL ^ Chem.SANITIZE_CLEANUP ^ Chem.SANITIZE_PROPERTIES))
    rdMolStandardize.NormalizeInPlace(mol)
    # kekulize the molecule
    # Chem.Kekulize(mol)
    # Update the properties
    mol.UpdatePropertyCache(strict=False)
    return mol


def side_to_dict(side):
    result = {}
    for component in map(str.strip, side.split('+')):
        match = re.match(r'(\d*)\s*(C\d+)', component)
        if match:
            count = int(match.group(1)) if match.group(1) else 1
            molecule = match.group(2)
            result[molecule] = count
    return result


def eq_to_dict(eq):
    return map(side_to_dict, eq.split('<=>'))


def fix_halogen_cid(data):
    target_dir = r"..\data\kegg_data_C"
    target_dir = os.path.abspath(target_dir)
    hal_exp = ['F', 'Cl', 'Br', 'I']
    # remove C13373
    # data.remove("C13373")  # Xe is not a halogen
    n_data = len(data)
    n_hal = len(hal_exp)
    n_comb = n_data * n_hal
    num_range = range(n_comb)
    # Reshape num_range to a 2D array
    num_range = [num_range[i:i + n_hal] for i in range(0, n_comb, n_hal)]
    cids = []
    smis = []
    smis_dict = {}
    cids_dict = {}
    for i in range(n_data):
        tmp_file = os.path.join(target_dir, data[i], data[i] + ".mol")
        with open(tmp_file, 'r') as f:
            file_data = f.read()
        # Initialize the list for the current data[i]
        cids_dict[data[i]] = []
        smis_dict[data[i]] = []
        # replace the X with the halogen
        for j, hal in enumerate(hal_exp):
            idx = num_range[i][j]
            mol = Chem.MolFromMolBlock(file_data.replace("X", hal))
            mol = standardize_mol(mol)
            smi = Chem.MolToSmiles(mol, allHsExplicit=True)
            # Determine the compound id
            cid = {
                "C00462": {"F": "C16487", "Cl": "C01327", "Br": "C13645", "I": "C05590"},
                "C01812": {"F": "C06108", "Cl": "C06755"}
            }.get(data[i], {}).get(hal)

            # Generate the compound id if not found
            if cid is None:
                cid = f"C{int(99000 + idx):05d}"
            # print(f"Compound {data[i]} with halogen {hal}, idx {idx} -> {cid}")
            cids_dict[data[i]].append(cid)
            smis_dict[data[i]].append(smi)
    return cids_dict, smis_dict


def fix_halogen_reactions():
    c_id_file = "data/kegg_data_C.csv.zip"
    r_id_file = "data/kegg_data_R.csv.zip"
    r_id_file = "data/atlas_data_kegg_R.csv.zip"
    # r_id_file = "data/atlas_data_R.csv.zip"
    hal_exp = ['F', 'Cl', 'Br', 'I']

    c_id_bad_file = "data/C_IDs_bad.dat"
    # Load the data
    data = load_bad_entries(c_id_bad_file, target_str="X group")
    # remove C13373
    data.remove("C13373")
    print(f"bad files {data}")
    cids_dict, smis_dict = fix_halogen_cid(data)

    # load the compounds data
    compounds = pd.read_csv(c_id_file, compression='zip')
    # # Get the compounds with the halogens
    # compounds = compounds[compounds['compound_id'].isin(data)]
    # print(compounds['compound_id'].values)
    target_dir = r"..\data\kegg_data_C"
    target_dir = os.path.abspath(target_dir)
    # Load the data
    reactions = pd.read_csv(r_id_file, compression='zip')
    re_set = set()
    for i in range(len(data)):
        # Get the equations
        equations = get_reactions_with_substring(reactions, data[i])
        re_set.update(equations['id'].values)
    re_set = list(re_set)
    print(f"Reactions with the halogens {re_set}")
    print(cids_dict)
    # loop over the reactions
    for i in range(len(re_set)):
        # Get the reaction
        reaction = reactions[reactions['id'] == re_set[i]]
        eq = reaction['reaction'].values
        # for each key value in the dictionary of cid
        print(eq)
        # break the eq
        lhs, rhs = eq_to_dict(eq[0])
        print(lhs, rhs)
        keys_list = list(lhs.keys())
        values_list = list(lhs.values())

        for i, key in enumerate(keys_list):
            if key in data:
                for j, val in enumerate(cids_dict[key]):
                    # # Get the index of the key
                    # idx = data.index(key)
                    # # Get the new value
                    new_value = cids_dict[key][j]
                    print(f"Replacing {key} with {new_value}")
                    # construct the new equation
                    lhs[new_value] = lhs.pop(key)
        #         print(f"Replacing {key} with {new_value}")
        #         # Update the lhs
        #         lhs[new_value] = lhs.pop(key)

        # # loop over the lhs
        # for key, value in lhs.items():
        #     for j, hal in enumerate(hal_exp):
        #         if key in data:
        #             # Get the index of the key
        #             idx = data.index(key)
        #             # Get the new value
        #             new_value = cids_dict[key][j]
        #             print(f"Replacing {key} with {new_value}")
        #             # Update the lhs
        #             lhs[new_value] = lhs.pop(key)
        #             # # Update the lhs
        #             # lhs[new_value] = lhs.pop(key)
        #     # if key in data:
        #     #     # Get the index of the key
        #     #     idx = data.index(key)
        #     #     # Get the new value
        #     #     new_value = cids_dict[key][hal_exp.index(value)]
        #     #     # Update the lhs
        #     #     lhs[new_value] = lhs.pop(key)


if __name__ == "__main__":
    print("Program started", flush=True)
    # main()
    fix_halogen_reactions()
    print("Program finished", flush=True)
