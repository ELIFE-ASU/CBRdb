import os

import pandas as pd
from rdkit import RDLogger

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

# suppress SettingWithCopyWarning
pd.options.mode.chained_assignment = None


def scan_for_kegg_pattern(s, prefix):
    """
    Scans a string for KEGG patterns with a specific prefix.

    Parameters:
    s (str): The string to scan for KEGG patterns.
    prefix (str): The prefix to look for in the KEGG patterns.

    Returns:
    list: A list of KEGG patterns that match the specified prefix.
    """
    if prefix in s:
        return [i for i in s.replace(';', ' ').replace('[', ' ').replace(']', ' ').split()
                if i.startswith(prefix) and len(i) == len(prefix) + 5]
    return []


def make_lists_printable(df):
    """
    Converts list-like columns in a DataFrame to printable strings.

    Parameters:
    df (DataFrame): The pandas DataFrame to be processed.

    Returns:
    DataFrame: A copy of the DataFrame with list-like columns converted to strings.
    """
    list_like_cols = df.loc[:, df.iloc[0].apply(lambda x: isinstance(x, (set, list)))]
    printable_df = df.copy(deep=True)
    printable_df[list_like_cols.columns] = list_like_cols.map(' '.join)
    return printable_df


def scan_for_terms(s):
    """
    Scans a string for specific terms defined in the global variable OK.

    Parameters:
    s (str): The string to scan for terms.

    Returns:
    str or None: The original string if any term from OK is found or if the string starts with 'R1' or 'R0', otherwise None.
    """
    for j in OK:
        if j in s or s.startswith(('R1', 'R0')):
            return s
    return None


def remove_terms(s):
    """
    Removes specific terms defined in the global variable OK from the string.

    Parameters:
    s (str): The string to remove terms from.

    Returns:
    str or None: None if any term from OK is found or if the string starts with 'R1' or 'R0', otherwise the original string.
    """
    for j in OK:
        if j in s or s.startswith(('R1', 'R0')):
            return None
    return s


def flag_missing_reactions(x):
    """
    Identifies missing reactions from a given list.

    Parameters:
    x (list): A list of reaction identifiers to check.

    Returns:
    float or list: NaN if no reactions are missing, otherwise a list of missing reactions.
    """
    missing_reactions = [i for i in x if i not in reactions.index]
    return float('NaN') if not missing_reactions else missing_reactions


def load_reactions_data(target_dir):
    """
    Loads reaction data from a specified directory.

    Parameters:
    target_dir (str): The directory containing the reaction data files.

    Returns:
    DataFrame: A pandas DataFrame containing the loaded reaction data.
    """
    a = os.listdir(target_dir)
    reactions = {}
    for i in a:
        file = os.path.abspath(f"{target_dir}/{i}/{i}.data")
        try:
            with open(file, 'r') as f:
                entries = {}
                for line in f.readlines()[:-1]:
                    if not line.startswith(' '):
                        field = line.split(' ')[0]
                        entries[field] = line.lstrip(field).lstrip().rstrip('\n')
                    else:
                        entries[field] += '; ' + line.lstrip().rstrip('\n')
                reactions[i] = entries
        except:
            print(f"Error reading {i}.data", flush=True)
            raise FileExistsError
    return pd.DataFrame(reactions).fillna('').T


def get_reaction_ids_substr(reactions, substr="incomplete reaction"):
    """
    Retrieves reaction IDs from a DataFrame where the 'COMMENT' column contains a specific substring.

    Parameters:
    reactions (DataFrame): The pandas DataFrame containing reaction data.
    substr (str): The substring to search for within the 'COMMENT' column. Default is "incomplete reaction".

    Returns:
    list: A list of reaction IDs where the 'COMMENT' column contains the specified substring.
    """
    incomplete_reaction_ids = reactions[
        reactions['COMMENT'].str.contains(substr, case=False, na=False)].index.tolist()
    return incomplete_reaction_ids


def clean_reaction_shortcuts(r_file='../data/kegg_data_R.csv.zip', data_dir='../data/'):
    """
    Cleans and processes KEGG reaction data, identifying and flagging multistep shortcuts, glycan participants,
    and reactions with incomplete or general data.

    Parameters:
    r_file (str): The path to the directory containing the preprocessed KEGG reaction csv.
    data_dir (str): The path to the directory where processed data files will be saved.

    Returns:
    None
    """
    # Prepare the full path of the files
    r_file = os.path.abspath(r_file)
    global reactions, OK

    f_save_files = False

    OK = ['STEP', 'REACTION', 'SIMILAR', 'TO', 'SEE', '+', r'R\d{5}']
    # Import reaction metadata. Note that "overall" column tags top-level reactions in 
    # [this table](https://www.kegg.jp/kegg/tables/br08210.html)
    reactions = pd.read_csv(r_file, header=0, index_col=0)
    reactions.columns = reactions.columns.str.upper()

    # Use comments to identify multi-step reactions
    # artificial shortcuts compressing multiple existing reaction IDs into a single net reaction.
    
    dm = reactions['COMMENT'].replace('', float('NaN')).dropna().str.upper().reset_index()
    dm['ORIGINAL_COMMENT'] = dm['COMMENT'].copy(deep=True)

    f_save_files = False

    OK = ['STEP', 'REACTION', 'SIMILAR', 'TO', 'SEE', '+', r'R\d{5}']
    
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
            .replace(dm.loc[i, 'id'] + '+', '')
            .replace('+' + dm.loc[i, 'id'], '')
            .replace(dm.loc[i, 'id'], '')
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
    has_rns_missing = dm.loc[dm['steps_shown'].apply(flag_missing_reactions).dropna().index, 'id'].to_list()
    incomplete_multistep = dm.loc[dm['id'].isin(has_rns_missing), 'id'].unique()
    print('Reactions with incomplete multistep shortcuts:', len(incomplete_multistep), incomplete_multistep, flush=True)

    # Filter out incomplete multistep reactions
    dm = dm.query('index not in @incomplete_multistep')
    dm['step_group'] = dm['steps_shown'].apply('+'.join)

    # Flag multistep shortcuts and generate data files
    reactions['MULTISTEP_FLAG'] = reactions.index.isin(dm['id'].unique())
    reactions['STEP_ENTRIES'] = dm.groupby('id')['step_group'].apply(list).fillna('')
    reactions_printable = make_lists_printable(reactions)
    if f_save_files:
        make_lists_printable(dm).to_csv(os.path.join(data_dir, 'reactions_multistep_intermediate_processing.csv.zip'),
                                        compression='zip',
                                        encoding='utf-8', index=False)
        reactions_printable.to_csv(os.path.join(data_dir, 'reactions_processed_full.csv.zip'),
                                   compression='zip',
                                   encoding='utf-8', index=False)

    reactions_processed = pd.concat(
        [reactions_printable.filter(like='FLAG'),
         reactions_printable.filter(like='ENTRIES')],
        axis=1)
    if f_save_files:
        reactions_processed.to_csv(os.path.join(data_dir, 'reactions_processed_basic.csv.zip'),
                                   compression='zip',
                                   encoding='utf-8', index=False)

    # Select only the reactions that are overall_flagged or glycan_flagged
    reactions_shortcut = reactions_printable.query('OVERALL.notna()')
    print('Reactions that are shortcuts:', len(reactions_shortcut), flush=True)
    reactions_shortcut = reactions_shortcut['ID'].to_list()
    reactions_shortcut = list(set([i.split()[0].strip() for i in reactions_shortcut]))

    reactions_incomplete = get_reaction_ids_substr(reactions, substr="incomplete reaction")
    print('Reactions with incomplete data:', len(reactions_incomplete), flush=True)

    reactions_general = get_reaction_ids_substr(reactions, substr="general reaction")
    print('General reactions:', len(reactions_general), flush=True)

    data = {
        'Reaction': reactions_shortcut + reactions_incomplete + reactions_general,
        'Type': ['shortcut'] * len(reactions_shortcut)
                + ['incomplete'] * len(reactions_incomplete)
                + ['general'] * len(reactions_general)
    }
    data = pd.DataFrame(data)
    # Sort the data
    data = data.sort_values(by=['Reaction'])
    data.to_csv(os.path.join(data_dir, 'R_IDs_bad.dat'), index=False)
    return data
